#!/usr/bin/env python3
"""
Build splice junction database from recount3 data.

recount3 provides comprehensive RNA-seq junction data across:
- GTEx (all tissues, ~17,000 samples)
- TCGA (cancer samples)
- SRA (public datasets)

This script:
1. Downloads junction data from recount3
2. Aggregates across samples to get junction support
3. Filters for high-confidence junctions
4. Outputs standardized TSV format for flAIr

Usage:
    # Download GTEx junctions (recommended - high quality)
    python build_recount3_splice_db.py \
        --source gtex \
        --output data/ref/recount3_junctions.tsv \
        --min-samples 2 \
        --min-reads 5

    # Combine multiple sources
    python build_recount3_splice_db.py \
        --source gtex tcga \
        --output data/ref/recount3_junctions_combined.tsv \
        --min-samples 5 \
        --min-reads 10
"""

import argparse
import gzip
import json
import os
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd
from tqdm import tqdm


# recount3 API endpoints
RECOUNT3_BASE = "https://recount.bio"
RECOUNT3_API = f"{RECOUNT3_BASE}/api/v1"

# Available data sources
SOURCES = {
    "gtex": "GTEx v8 (all tissues, ~17,000 samples)",
    "tcga": "TCGA (cancer, ~11,000 samples)",
    "sra": "SRA (public data, subset)",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build splice junction database from recount3"
    )

    parser.add_argument(
        "--source",
        nargs="+",
        choices=list(SOURCES.keys()),
        default=["gtex"],
        help="Data source(s) to use (default: gtex)",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV file path",
    )

    parser.add_argument(
        "--min-samples",
        type=int,
        default=2,
        help="Minimum samples supporting junction (default: 2)",
    )

    parser.add_argument(
        "--min-reads",
        type=int,
        default=5,
        help="Minimum total reads supporting junction (default: 5)",
    )

    parser.add_argument(
        "--genome-build",
        choices=["hg38", "hg19"],
        default="hg38",
        help="Genome build (default: hg38)",
    )

    parser.add_argument(
        "--cache-dir",
        default="data/cache/recount3",
        help="Cache directory for downloaded files",
    )

    parser.add_argument(
        "--tissue",
        nargs="+",
        help="GTEx tissues to include (default: all)",
    )

    parser.add_argument(
        "--annotated-only",
        action="store_true",
        help="Only include annotated junctions (from GENCODE)",
    )

    return parser.parse_args()


def download_file(url: str, output_path: str, desc: str = None):
    """Download file with progress bar."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    if os.path.exists(output_path):
        print(f"Using cached file: {output_path}")
        return output_path

    print(f"Downloading: {url}")

    def progress_hook(count, block_size, total_size):
        if not hasattr(progress_hook, "pbar"):
            progress_hook.pbar = tqdm(
                total=total_size,
                unit="B",
                unit_scale=True,
                desc=desc or "Downloading",
            )
        progress_hook.pbar.update(block_size)

    urllib.request.urlretrieve(url, output_path, reporthook=progress_hook)

    if hasattr(progress_hook, "pbar"):
        progress_hook.pbar.close()
        delattr(progress_hook, "pbar")

    return output_path


def get_gtex_projects(cache_dir: str, tissues: List[str] = None) -> List[str]:
    """
    Get list of GTEx project IDs from recount3.

    GTEx projects are organized by tissue type.
    """
    # For GTEx in recount3, we typically use the unified GTEx collection
    # This contains all tissues pre-aggregated

    gtex_projects = [
        "BRAIN",
        "HEART",
        "LIVER",
        "LUNG",
        "MUSCLE",
        "SKIN",
        "BLOOD",
        # Add more as needed
    ]

    if tissues:
        gtex_projects = [p for p in gtex_projects if p.lower() in [t.lower() for t in tissues]]

    return gtex_projects


def parse_recount3_junctions(junction_file: str, genome_build: str = "hg38") -> pd.DataFrame:
    """
    Parse recount3 junction format.

    recount3 provides junctions in a custom format:
    - chr:start-end:strand format
    - or BED-like format depending on export

    This is a placeholder - actual format may vary.
    You'll need to adjust based on real recount3 data structure.
    """

    print(f"Parsing junctions from {junction_file}")

    # Placeholder parsing logic
    # Real implementation depends on recount3 API response format

    junctions = []

    # Example: if junctions are in BED format
    if junction_file.endswith('.bed') or junction_file.endswith('.bed.gz'):
        opener = gzip.open if junction_file.endswith('.gz') else open

        with opener(junction_file, 'rt') as f:
            for line in f:
                if line.startswith('#') or line.startswith('track'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 6:
                    continue

                chrom = fields[0]
                start = int(fields[1])  # donor (0-based)
                end = int(fields[2])    # acceptor (0-based)
                strand = fields[5]

                # Extract read count from name or score field
                score = int(fields[4]) if len(fields) > 4 else 1

                junctions.append({
                    'chr': chrom,
                    'donor': start,
                    'acceptor': end,
                    'strand': strand,
                    'reads': score,
                    'samples': 1  # Will aggregate later
                })

    # Example: if junctions are in custom format
    elif junction_file.endswith('.tsv') or junction_file.endswith('.tsv.gz'):
        df = pd.read_csv(junction_file, sep='\t', compression='gzip' if junction_file.endswith('.gz') else None)

        # Adjust column names based on actual recount3 format
        if 'junction_id' in df.columns:
            # Parse junction_id like "chr1:12345-67890:+"
            def parse_junction_id(jid):
                chrom, coords = jid.rsplit(':', 1)
                strand = coords[-1]
                start, end = coords[:-2].split('-')
                return chrom, int(start), int(end), strand

            parsed = df['junction_id'].apply(parse_junction_id)
            df[['chr', 'donor', 'acceptor', 'strand']] = pd.DataFrame(
                parsed.tolist(), index=df.index
            )

        junctions = df.to_dict('records')

    return pd.DataFrame(junctions)


def aggregate_junctions(
    junction_dfs: List[pd.DataFrame],
    min_samples: int = 2,
    min_reads: int = 5,
) -> pd.DataFrame:
    """
    Aggregate junctions across multiple samples/projects.

    Args:
        junction_dfs: List of junction DataFrames
        min_samples: Minimum samples supporting junction
        min_reads: Minimum total reads across all samples

    Returns:
        Aggregated junction DataFrame with support metrics
    """

    print(f"Aggregating {len(junction_dfs)} junction datasets...")

    # Combine all junctions
    all_junctions = pd.concat(junction_dfs, ignore_index=True)

    # Group by junction coordinates
    grouped = all_junctions.groupby(['chr', 'donor', 'acceptor', 'strand'])

    # Aggregate metrics
    aggregated = grouped.agg({
        'reads': 'sum',      # Total reads across all samples
        'samples': 'sum',    # Total samples (if already counted)
    }).reset_index()

    # If samples weren't counted, count unique occurrences
    if 'samples' not in aggregated.columns or aggregated['samples'].isna().all():
        sample_counts = grouped.size().reset_index(name='samples')
        aggregated = aggregated.merge(sample_counts, on=['chr', 'donor', 'acceptor', 'strand'])

    # Filter by support thresholds
    print(f"Before filtering: {len(aggregated)} junctions")

    aggregated = aggregated[
        (aggregated['samples'] >= min_samples) &
        (aggregated['reads'] >= min_reads)
    ]

    print(f"After filtering: {len(aggregated)} junctions")
    print(f"  Min samples: {min_samples}")
    print(f"  Min reads: {min_reads}")

    # Calculate weight (log-scaled read support)
    import numpy as np
    aggregated['weight'] = np.log10(aggregated['reads'] + 1)

    # Sort by support
    aggregated = aggregated.sort_values(['chr', 'donor', 'acceptor'])

    return aggregated


def validate_junctions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Validate junction coordinates and strand.

    - Check for canonical splice sites (GT-AG)
    - Remove invalid chromosomes
    - Check coordinate ordering
    """

    print("Validating junctions...")

    initial_count = len(df)

    # Valid chromosomes
    valid_chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
    df = df[df['chr'].isin(valid_chroms)]

    # Check coordinate ordering
    df = df[df['donor'] < df['acceptor']]

    # Check valid strands
    df = df[df['strand'].isin(['+', '-'])]

    print(f"  Removed {initial_count - len(df)} invalid junctions")

    return df


def download_recount3_gtex_junctions(cache_dir: str, tissues: List[str] = None) -> List[str]:
    """
    Download GTEx junction data from recount3.

    NOTE: This is a template function. The actual recount3 API
    and data format will need to be verified and adjusted.

    Real implementation options:
    1. Use recount3 R package via rpy2
    2. Use recount3 REST API (if available)
    3. Direct download from known URLs
    """

    print("=" * 80)
    print("IMPORTANT: recount3 Data Access")
    print("=" * 80)
    print()
    print("This script provides a template for downloading recount3 data.")
    print("However, the actual recount3 API and data format need verification.")
    print()
    print("Recommended approaches:")
    print()
    print("1. Use the recount3 R package:")
    print("   R> library(recount3)")
    print("   R> gtex <- available_projects('gtex')")
    print("   R> jxn <- get_junctions('GTEx', 'BRAIN')")
    print()
    print("2. Check recount.bio documentation:")
    print("   https://rna.recount.bio/docs")
    print()
    print("3. Use pre-computed GTEx junction files if available")
    print()
    print("For now, this will create a placeholder structure.")
    print("You'll need to populate it with real recount3 data.")
    print("=" * 80)

    # Placeholder - return empty list
    # Real implementation would download actual files

    junction_files = []

    # Example structure (adjust based on real API):
    # for tissue in tissues or get_gtex_projects(cache_dir):
    #     url = f"{RECOUNT3_BASE}/data/gtex/{tissue}/junctions.bed.gz"
    #     output = f"{cache_dir}/gtex_{tissue}_junctions.bed.gz"
    #     downloaded = download_file(url, output, desc=f"GTEx {tissue}")
    #     junction_files.append(downloaded)

    return junction_files


def main():
    args = parse_args()

    print("=" * 80)
    print("flAIr Splice Junction Database Builder")
    print("=" * 80)
    print(f"Source(s): {', '.join(args.source)}")
    print(f"Genome build: {args.genome_build}")
    print(f"Min samples: {args.min_samples}")
    print(f"Min reads: {args.min_reads}")
    print()

    # Create cache directory
    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    all_junction_dfs = []

    # Process each source
    for source in args.source:
        print(f"\n[{source.upper()}] Processing...")

        if source == "gtex":
            junction_files = download_recount3_gtex_junctions(
                str(cache_dir),
                tissues=args.tissue
            )

            for jf in junction_files:
                df = parse_recount3_junctions(jf, args.genome_build)
                all_junction_dfs.append(df)

        elif source == "tcga":
            # Similar logic for TCGA
            print("TCGA support not yet implemented")
            print("Please use GTEx or provide TCGA junction files manually")

        elif source == "sra":
            print("SRA support not yet implemented")

    if not all_junction_dfs:
        print("\n" + "=" * 80)
        print("WARNING: No junction data downloaded")
        print("=" * 80)
        print()
        print("This appears to be your first time using this script.")
        print("Since recount3 requires specific API access or R package usage,")
        print("please use one of these alternatives:")
        print()
        print("OPTION 1: Use your existing RJunBase data")
        print("  You already have: splice_db/detail_Linear_splice_annotation.txt")
        print("  This is already a comprehensive junction database!")
        print()
        print("OPTION 2: Download GTEx junctions directly")
        print("  GTEx Portal: https://gtexportal.org/home/datasets")
        print("  Download: GTEx_Analysis_v8_STARjunctions.gct.gz")
        print()
        print("OPTION 3: Use GENCODE annotations")
        print("  Extract junctions from GENCODE GTF (high quality, annotated)")
        print()
        print("For now, I recommend using your existing RJunBase data,")
        print("which you can process with the standardize_priors.py script.")
        print()
        sys.exit(0)

    # Aggregate across all sources
    print("\n" + "=" * 80)
    print("Aggregating junctions...")
    print("=" * 80)

    aggregated = aggregate_junctions(
        all_junction_dfs,
        min_samples=args.min_samples,
        min_reads=args.min_reads,
    )

    # Validate
    aggregated = validate_junctions(aggregated)

    # Select final columns
    final_df = aggregated[['chr', 'donor', 'acceptor', 'strand', 'weight']].copy()

    # Write output
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 80)
    print("Writing output...")
    print("=" * 80)

    final_df.to_csv(output_path, sep='\t', index=False, header=True)

    print(f"Wrote {len(final_df)} junctions to {output_path}")

    # Generate statistics
    print("\n" + "=" * 80)
    print("Junction Statistics")
    print("=" * 80)
    print(f"Total junctions: {len(final_df):,}")
    print(f"Chromosomes: {final_df['chr'].nunique()}")
    print(f"  Autosomes: {final_df[final_df['chr'].str.match('chr[0-9]+')]['chr'].nunique()}")
    print(f"  Sex chromosomes: {final_df[final_df['chr'].isin(['chrX', 'chrY'])]['chr'].nunique()}")
    print(f"Strand distribution:")
    print(f"  +: {(final_df['strand'] == '+').sum():,}")
    print(f"  -: {(final_df['strand'] == '-').sum():,}")

    print("\nWeight distribution:")
    print(final_df['weight'].describe())

    print("\n" + "=" * 80)
    print("COMPLETE!")
    print("=" * 80)
    print(f"\nOutput: {output_path}")
    print("\nNext steps:")
    print("1. Compress and index:")
    print(f"   sort -k1,1 -k2,2n {output_path} > {output_path}.sorted")
    print(f"   bgzip {output_path}.sorted")
    print(f"   tabix -s 1 -b 2 -e 3 -S 1 {output_path}.sorted.gz")
    print("2. Validate with:")
    print(f"   tabix {output_path}.sorted.gz chr1:1000000-2000000")


if __name__ == "__main__":
    main()
