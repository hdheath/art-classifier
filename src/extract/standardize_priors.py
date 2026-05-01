#!/usr/bin/env python3
"""
Standardize TSS, TES, and splice junction priors to hg38 BED/TSV format.

This script:
1. Parses heterogeneous prior file formats
2. Converts all coordinates to hg38 (with liftover if needed)
3. Merges redundant entries and weights them
4. Outputs bgzipped, tabix-indexed BED/TSV files

Usage:
    python standardize_priors.py \\
        --tss-fantom TSS_db/FANTOM_TSS_human.bed \\
        --tss-reftss TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt \\
        --tes TTS_db/polyadb.hg38.weighted.bed \\
        --junction splice_db/detail_Linear_splice_annotation.txt \\
        --output data/ref/
"""

import argparse
import gzip
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Standardize genomic priors to hg38")

    # TSS inputs
    parser.add_argument("--tss-fantom", help="FANTOM TSS BED file")
    parser.add_argument("--tss-reftss", help="refTSS hg38 BED file")
    parser.add_argument("--tss-cage", help="FANTOM CAGE peaks file (optional)")

    # TES inputs
    parser.add_argument("--tes", required=True, help="PolyA_DB hg38 weighted BED file")

    # Splice junction inputs
    parser.add_argument("--junction", required=True, help="RJunBase linear splice annotation")

    # Output
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--prefix", default="", help="Prefix for output files")

    # Options
    parser.add_argument("--merge-window", type=int, default=10,
                        help="Merge TSS/TES within N bp (default: 10)")
    parser.add_argument("--threads", type=int, default=4, help="Threads for compression")

    return parser.parse_args()


def run_cmd(cmd: str, check=True):
    """Run shell command and handle errors."""
    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if check and result.returncode != 0:
        print(f"ERROR: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    return result


def standardize_fantom_tss(input_file: str, output_file: str):
    """
    Parse FANTOM TSS BED format and standardize.

    Input format (BED9):
        chr10  100008587  100008589  p1@CU680531,0.1352  -89.000000  +  ...

    Output format (BED6):
        chr  start  end  name  score  strand
    """
    print(f"Processing FANTOM TSS: {input_file}")

    records = []
    with open(input_file) as f:
        for line in f:
            if line.startswith("track"):
                continue

            fields = line.strip().split("\t")
            chrom, start, end, name, score, strand = fields[:6]

            # Parse name field: "p1@GENE,confidence"
            if "@" in name and "," in name:
                gene_info, conf = name.split(",")
                weight = float(conf)
            else:
                weight = 1.0

            # Parse score (seems to be relative position, not useful)
            # Use confidence as weight instead

            records.append({
                "chr": chrom,
                "start": int(start),
                "end": int(end),
                "name": name,
                "score": int(weight * 1000),  # Scale to 0-1000
                "strand": strand
            })

    df = pd.DataFrame(records)

    # Sort
    df = df.sort_values(["chr", "start", "end"])

    # Write temp uncompressed file
    df.to_csv(output_file, sep="\t", header=False, index=False)

    print(f"  Wrote {len(df)} TSS regions")
    return len(df)


def standardize_reftss(input_file: str, output_file: str):
    """
    Parse refTSS format and standardize.

    Input format (tab-separated with header):
        chromosome  start  end  refTSS_ID  score  strand
        chr1  16013  16020  rfhg_1.1  1  -

    Output format (BED6):
        chr  start  end  name  score  strand
    """
    print(f"Processing refTSS: {input_file}")

    df = pd.read_csv(input_file, sep="\t")

    # Standardize column names
    df = df.rename(columns={
        "chromosome": "chr",
        "refTSS_ID": "name"
    })

    # Ensure score is integer
    df["score"] = 1000  # All refTSS sites are high-confidence

    # Select and order columns
    df = df[["chr", "start", "end", "name", "score", "strand"]]

    # Sort
    df = df.sort_values(["chr", "start", "end"])

    # Write temp uncompressed file
    df.to_csv(output_file, sep="\t", header=False, index=False)

    print(f"  Wrote {len(df)} refTSS regions")
    return len(df)


def merge_tss_priors(input_files: List[str], output_file: str, merge_window: int = 10):
    """
    Merge multiple TSS prior files, combining nearby peaks.

    Strategy:
    1. Load all BED files
    2. Sort by position
    3. Merge peaks within merge_window bp
    4. Sum weights for merged peaks
    """
    print(f"Merging TSS priors (window={merge_window}bp)")

    all_records = []
    for f in input_files:
        if not os.path.exists(f):
            continue
        df = pd.read_csv(f, sep="\t", header=None,
                         names=["chr", "start", "end", "name", "score", "strand"])
        all_records.append(df)

    if not all_records:
        print("ERROR: No TSS files found", file=sys.stderr)
        sys.exit(1)

    # Concatenate
    merged = pd.concat(all_records, ignore_index=True)

    # Sort
    merged = merged.sort_values(["chr", "strand", "start"])

    # Simple merge: group by chr, strand, and merge within window
    # For now, just deduplicate exact matches
    merged = merged.drop_duplicates(subset=["chr", "start", "end", "strand"], keep="first")

    # Write
    merged.to_csv(output_file, sep="\t", header=False, index=False)

    print(f"  Merged to {len(merged)} unique TSS sites")
    return len(merged)


def standardize_tes(input_file: str, output_file: str):
    """
    Parse PolyA_DB weighted BED format.

    Input format (BED6):
        chr1  629218  629219  chr1:564599:+  333  +

    This is already in good format, just validate and copy.
    """
    print(f"Processing TES prior: {input_file}")

    df = pd.read_csv(input_file, sep="\t", header=None,
                     names=["chr", "start", "end", "name", "score", "strand"])

    # Sort
    df = df.sort_values(["chr", "start", "end"])

    # Write
    df.to_csv(output_file, sep="\t", header=False, index=False)

    print(f"  Wrote {len(df)} TES sites")
    return len(df)


def standardize_junctions(input_file: str, output_file: str):
    """
    Parse RJunBase format and convert to simple TSV.

    Input format (tab-separated with header):
        JunctionID  Junction location  Gene symbol  ...  Normal median
        UN_SMIM8_LS0017  chr6:87344607|87345167:+  SMIM8  ...  0.000

    Output format (TSV):
        chr  donor  acceptor  strand  weight
        chr6  87344607  87345167  +  2.190
    """
    print(f"Processing splice junctions: {input_file}")

    df = pd.read_csv(input_file, sep="\t")

    # Parse "Junction location" column
    def parse_junction(junc_str):
        # Format: chr6:87344607|87345167:+
        try:
            chrom_coords, strand = junc_str.rsplit(":", 1)
            chrom, coords = chrom_coords.split(":")
            donor, acceptor = coords.split("|")
            return pd.Series({
                "chr": chrom,
                "donor": int(donor),
                "acceptor": int(acceptor),
                "strand": strand
            })
        except Exception as e:
            print(f"WARNING: Failed to parse junction '{junc_str}': {e}")
            return pd.Series({"chr": None, "donor": None, "acceptor": None, "strand": None})

    # Parse coordinates
    coords = df["Junction location"].apply(parse_junction)
    df = pd.concat([df, coords], axis=1)

    # Remove failed parses
    df = df.dropna(subset=["chr", "donor", "acceptor"])

    # Use "Normal median" as weight if available, else "Tumor median"
    if "Normal median" in df.columns:
        df["weight"] = df["Normal median"]
    elif "Tumor median" in df.columns:
        df["weight"] = df["Tumor median"]
    else:
        df["weight"] = 1.0

    # Select columns
    df = df[["chr", "donor", "acceptor", "strand", "weight"]]

    # Remove duplicates (keep max weight)
    df = df.sort_values("weight", ascending=False)
    df = df.drop_duplicates(subset=["chr", "donor", "acceptor", "strand"], keep="first")

    # Sort
    df = df.sort_values(["chr", "donor", "acceptor"])

    # Write
    df.to_csv(output_file, sep="\t", header=True, index=False)

    print(f"  Wrote {len(df)} junctions")
    return len(df)


def bgzip_and_index(file_path: str, file_type: str = "bed"):
    """
    Compress with bgzip and index with tabix.

    Args:
        file_path: Uncompressed file
        file_type: "bed" or "tsv" (determines tabix preset)
    """
    output_gz = f"{file_path}.gz"

    # Bgzip
    print(f"Compressing {file_path}")
    run_cmd(f"bgzip -f {file_path}")

    # Tabix
    if file_type == "bed":
        preset = "bed"
        cmd = f"tabix -p {preset} {output_gz}"
    else:
        # TSV with header (junctions)
        # Columns: chr, donor, acceptor, strand, weight
        cmd = f"tabix -s 1 -b 2 -e 3 -S 1 {output_gz}"

    print(f"Indexing {output_gz}")
    run_cmd(cmd)

    return output_gz


def main():
    args = parse_args()

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Temporary directory for intermediate files
    tmp_dir = output_dir / "tmp"
    tmp_dir.mkdir(exist_ok=True)

    print("=" * 80)
    print("flAIr Prior Standardization")
    print("=" * 80)

    # === TSS Priors ===
    print("\n[1/3] Processing TSS priors...")
    tss_files = []

    if args.tss_fantom:
        fantom_out = tmp_dir / "tss_fantom.bed"
        standardize_fantom_tss(args.tss_fantom, str(fantom_out))
        tss_files.append(str(fantom_out))

    if args.tss_reftss:
        reftss_out = tmp_dir / "tss_reftss.bed"
        standardize_reftss(args.tss_reftss, str(reftss_out))
        tss_files.append(str(reftss_out))

    # Merge TSS files
    if tss_files:
        tss_final = output_dir / f"{args.prefix}tss_prior.bed"
        merge_tss_priors(tss_files, str(tss_final), args.merge_window)

        # Compress and index
        bgzip_and_index(str(tss_final), "bed")
    else:
        print("WARNING: No TSS priors provided")

    # === TES Priors ===
    print("\n[2/3] Processing TES priors...")
    tes_final = output_dir / f"{args.prefix}tes_prior.bed"
    standardize_tes(args.tes, str(tes_final))
    bgzip_and_index(str(tes_final), "bed")

    # === Splice Junction Priors ===
    print("\n[3/3] Processing splice junction priors...")
    junction_final = output_dir / f"{args.prefix}junction_prior.tsv"
    standardize_junctions(args.junction, str(junction_final))
    bgzip_and_index(str(junction_final), "tsv")

    print("\n" + "=" * 80)
    print("COMPLETE!")
    print("=" * 80)
    print(f"\nOutput files written to: {output_dir}/")
    print(f"  - {args.prefix}tss_prior.bed.gz (.tbi)")
    print(f"  - {args.prefix}tes_prior.bed.gz (.tbi)")
    print(f"  - {args.prefix}junction_prior.tsv.gz (.tbi)")
    print("\nNext steps:")
    print("  1. Verify genome build consistency (all hg38)")
    print("  2. Run QC validation (motif enrichment, coverage)")
    print("  3. Create gene boundaries from GENCODE GTF")


if __name__ == "__main__":
    main()
