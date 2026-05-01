#!/usr/bin/env python3
"""
Extract TSS coordinates from FANTOM5 CAGE peak file.

The FANTOM CAGE file contains:
- Genomic coordinates in first column (chr:start..end,strand)
- TPM expression values across ~1,800 samples
- Already normalized (RLE normalization)

This script:
1. Parses CAGE peak coordinates
2. Computes aggregate expression (mean/max TPM across samples)
3. Outputs BED6 format with weight
4. Optionally lifts hg19 → hg38

Usage:
    python parse_fantom_cage.py \\
        --input TSS_db/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz \\
        --output data/ref/fantom_cage_tss.bed \\
        --genome-build hg19 \\
        --liftover hg19ToHg38.over.chain
"""

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description="Parse FANTOM5 CAGE peaks")

    parser.add_argument("--input", required=True,
                        help="FANTOM CAGE file (hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz)")
    parser.add_argument("--output", required=True,
                        help="Output BED file")
    parser.add_argument("--genome-build", default="hg19",
                        choices=["hg19", "hg38"],
                        help="Input genome build (default: hg19)")
    parser.add_argument("--liftover",
                        help="Chain file for liftover (hg19ToHg38.over.chain)")
    parser.add_argument("--aggregation",
                        choices=["mean", "max", "median"],
                        default="mean",
                        help="How to aggregate TPM across samples (default: mean)")
    parser.add_argument("--min-tpm", type=float, default=0.5,
                        help="Minimum TPM to keep peak (default: 0.5)")

    return parser.parse_args()


def parse_cage_coordinate(coord_str):
    """
    Parse FANTOM CAGE coordinate format.

    Format: chr10:100013403..100013414,-

    Returns: (chr, start, end, strand)
    """
    if coord_str.startswith("0") or "STAT" in coord_str:
        return None, None, None, None

    try:
        # Split chr and coords
        chrom, rest = coord_str.rsplit(":", 1)

        # Split coords and strand
        coords, strand = rest.rsplit(",", 1)

        # Parse start..end
        start, end = coords.split("..")

        return chrom, int(start), int(end), strand

    except Exception as e:
        print(f"Warning: Failed to parse '{coord_str}': {e}", file=sys.stderr)
        return None, None, None, None


def aggregate_tpm(tpm_values, method="mean"):
    """Aggregate TPM across samples."""
    if method == "mean":
        return np.mean(tpm_values)
    elif method == "max":
        return np.max(tpm_values)
    elif method == "median":
        return np.median(tpm_values)
    else:
        raise ValueError(f"Unknown aggregation method: {method}")


def parse_fantom_cage(input_file, aggregation="mean", min_tpm=0.5):
    """
    Parse FANTOM CAGE file and extract TSS regions.

    Returns DataFrame with columns: chr, start, end, name, score, strand
    """

    print(f"Parsing FANTOM CAGE file: {input_file}")
    print(f"Aggregation method: {aggregation}")
    print(f"Minimum TPM: {min_tpm}")

    # Read file
    opener = gzip.open if input_file.endswith('.gz') else open

    peaks = []
    header_cols = None

    with opener(input_file, 'rt') as f:
        for line_num, line in enumerate(f):
            # Skip metadata lines
            if line.startswith("##"):
                continue

            # Parse header (column names)
            if line_num == 0 or header_cols is None:
                if not line.startswith("##"):
                    header_cols = line.strip().split('\t')
                    print(f"Found {len(header_cols)} columns")
                    print(f"Sample columns: {len(header_cols) - 7} (excluding annotation columns)")
                    continue

            fields = line.strip().split('\t')

            # Parse coordinate
            coord_str = fields[0]
            chrom, start, end, strand = parse_cage_coordinate(coord_str)

            if chrom is None:
                continue  # Skip stats rows

            # Get TPM values (columns 8 onwards are sample TPM values)
            tpm_values = []
            for val in fields[7:]:  # Skip first 7 annotation columns
                try:
                    tpm_values.append(float(val))
                except ValueError:
                    continue

            if not tpm_values:
                continue

            # Aggregate TPM
            agg_tpm = aggregate_tpm(tpm_values, aggregation)

            # Filter by minimum TPM
            if agg_tpm < min_tpm:
                continue

            # Create peak entry
            peaks.append({
                'chr': chrom,
                'start': start,
                'end': end,
                'name': coord_str,
                'score': int(agg_tpm * 100),  # Scale to 0-1000+ range
                'strand': strand,
                'tpm': agg_tpm
            })

            if (line_num + 1) % 10000 == 0:
                print(f"Processed {line_num + 1} lines, kept {len(peaks)} peaks...")

    print(f"\nTotal peaks kept: {len(peaks)}")

    df = pd.DataFrame(peaks)

    # Sort by position
    df = df.sort_values(['chr', 'start', 'end'])

    return df


def run_liftover(bed_file, chain_file, output_file):
    """
    Run liftOver to convert hg19 → hg38.

    Requires liftOver binary in PATH or TTS_db/ directory.
    """
    import subprocess

    print(f"\nLifting over from hg19 to hg38...")
    print(f"Chain file: {chain_file}")

    # Check for liftOver binary
    liftover_binary = None
    for location in ["liftOver", "TTS_db/liftOver", "./liftOver"]:
        try:
            result = subprocess.run([location, "--version"],
                                   capture_output=True, text=True)
            if result.returncode == 0 or "liftOver" in result.stderr:
                liftover_binary = location
                break
        except FileNotFoundError:
            continue

    if not liftover_binary:
        print("ERROR: liftOver not found!", file=sys.stderr)
        print("Download from: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver")
        print("Or use existing: TTS_db/liftOver")
        sys.exit(1)

    print(f"Using liftOver: {liftover_binary}")

    # Run liftOver
    unmapped_file = bed_file + ".unmapped"

    cmd = [
        liftover_binary,
        bed_file,
        chain_file,
        output_file,
        unmapped_file
    ]

    print(f"Running: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"liftOver failed: {result.stderr}", file=sys.stderr)
        sys.exit(1)

    # Count results
    with open(output_file) as f:
        lifted = sum(1 for line in f)

    with open(unmapped_file) as f:
        unmapped = sum(1 for line in f if not line.startswith('#'))

    print(f"Liftover complete:")
    print(f"  Lifted: {lifted} peaks")
    print(f"  Unmapped: {unmapped} peaks")
    print(f"  Success rate: {lifted/(lifted+unmapped)*100:.1f}%")

    return output_file


def main():
    args = parse_args()

    print("=" * 80)
    print("FANTOM5 CAGE Peak Parser")
    print("=" * 80)
    print()

    # Parse CAGE file
    df = parse_fantom_cage(
        args.input,
        aggregation=args.aggregation,
        min_tpm=args.min_tpm
    )

    # Select BED6 columns
    bed_df = df[['chr', 'start', 'end', 'name', 'score', 'strand']]

    # Write temporary file
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    temp_file = str(output_path) + ".temp"
    bed_df.to_csv(temp_file, sep='\t', header=False, index=False)

    print(f"\nWrote {len(bed_df)} peaks to {temp_file}")

    # Liftover if needed
    if args.genome_build == "hg19" and args.liftover:
        final_file = run_liftover(temp_file, args.liftover, args.output)

        # Remove temp file
        Path(temp_file).unlink()

    elif args.genome_build == "hg19":
        print("\n" + "=" * 80)
        print("WARNING: Input is hg19, but no liftover chain provided!")
        print("=" * 80)
        print()
        print("Your data is in hg19, but flAIr needs hg38.")
        print()
        print("To convert to hg38:")
        print(f"  python {__file__} \\")
        print(f"      --input {args.input} \\")
        print(f"      --output {args.output} \\")
        print("      --liftover TTS_db/hg19ToHg38.over.chain")
        print()
        print("For now, saving hg19 coordinates to:", temp_file)
        print()

        # Rename temp to output
        Path(temp_file).rename(args.output)
        final_file = args.output

    else:
        # Already hg38
        Path(temp_file).rename(args.output)
        final_file = args.output

    # Summary statistics
    final_df = pd.read_csv(final_file, sep='\t', header=None,
                            names=['chr', 'start', 'end', 'name', 'score', 'strand'])

    print("\n" + "=" * 80)
    print("Summary Statistics")
    print("=" * 80)
    print(f"Total TSS peaks: {len(final_df):,}")
    print(f"Chromosomes: {final_df['chr'].nunique()}")
    print(f"Strand distribution:")
    print(f"  +: {(final_df['strand'] == '+').sum():,}")
    print(f"  -: {(final_df['strand'] == '-').sum():,}")
    print()
    print("Score (expression) distribution:")
    print(final_df['score'].describe())
    print()
    print("=" * 80)
    print("COMPLETE!")
    print("=" * 80)
    print(f"\nOutput: {final_file}")
    print()
    print("Next steps:")
    print("1. Merge with refTSS and FANTOM_TSS_human.bed:")
    print(f"   python src/extract/standardize_priors.py \\")
    print(f"       --tss-cage {final_file} \\")
    print("       --tss-reftss TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt \\")
    print("       --tss-fantom TSS_db/FANTOM_TSS_human.bed \\")
    print("       --output data/ref/")
    print()
    print("2. Or skip this CAGE file and just use the simpler BED files you have:")
    print("   (FANTOM_TSS_human.bed + refTSS are already comprehensive)")


if __name__ == "__main__":
    main()
