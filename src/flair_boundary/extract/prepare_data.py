#!/usr/bin/env python3
"""
Prepare training datasets for TSS and TTS boundary classifiers.

This script:
1. Loads positive sites from reference databases (FANTOM, refTSS, polyAdb)
2. Merges and deduplicates within a distance window
3. Generates matched negative (non-boundary) sites
4. Extracts sequence features from GRCh38 genome
5. Saves feature matrices as compressed .npz for fast training

Usage:
    # Prepare TSS training data
    python prepare_data.py \\
        --type tss \\
        --genome /path/to/GRCh38.primary_assembly.genome.fa \\
        --fantom /path/to/TSS_db/FANTOM_TSS_human.bed \\
        --reftss /path/to/TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt \\
        --output /path/to/data/derived/tss_features.npz \\
        --neg-ratio 3 \\
        --max-sites 200000

    # Prepare TTS training data
    python prepare_data.py \\
        --type tts \\
        --genome /path/to/GRCh38.primary_assembly.genome.fa \\
        --polyadb /path/to/TTS_db/polyadb.hg38.weighted.bed \\
        --polyadb-pas /path/to/TTS_db/polyAdb_human.PAS.txt \\
        --output /path/to/data/derived/tts_features.npz \\
        --neg-ratio 3
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# Allow running from this directory
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from flair_boundary.extract.features import (
    SequenceFeatureExtractor,
    generate_hard_negatives,
    generate_negative_sites,
    load_positive_sites,
)


def parse_args():
    p = argparse.ArgumentParser(description="Prepare boundary classifier training data")

    p.add_argument("--type", required=True, choices=["tss", "tts"],
                   help="Boundary type to prepare")
    p.add_argument("--genome", required=True,
                   help="GRCh38 genome FASTA (must be indexed with samtools faidx)")
    p.add_argument("--output", required=True,
                   help="Output .npz file for features and labels")

    # TSS inputs
    p.add_argument("--fantom",
                   help="FANTOM permissive CAGE peaks BED9 (hg38). Used as the "
                        "primary positives source. Must be verified hg38 — see "
                        "qc/database_qc_report.md for the build audit.")
    p.add_argument("--fantom-fair",
                   help="FANTOM5 hg38 'fair' robust CAGE peaks BED9 (hg38-native, "
                        "from Lizio et al. 2017). Used as the width target via "
                        "spatial overlap with FANTOM positives — peaks that fall "
                        "inside a fair cluster inherit the cluster's width "
                        "(end-start) as the regression target.")
    p.add_argument("--fantom-fair-max-width", type=int, default=5000,
                   help="Drop fair-peak clusters wider than this (bp) — guards "
                        "against rare liftOver pathologies.")

    # TTS inputs
    p.add_argument("--polyadb",
                   help="PolyA_DB BED6 file (hg38)")
    p.add_argument("--polyadb-pas",
                   help="PolyA_DB PAS annotation file (optional, for score filtering)")

    # Options
    p.add_argument("--neg-ratio", type=int, default=3,
                   help="Negatives per positive (default: 3)")
    p.add_argument("--max-sites", type=int, default=0,
                   help="Max positive sites to use (0 = all)")
    p.add_argument("--min-score", type=float, default=0,
                   help="Minimum BED score for positive sites")
    p.add_argument("--merge-window", type=int, default=25,
                   help="Merge sites within this distance (bp, default: 25)")
    p.add_argument("--neg-min-dist", type=int, default=500,
                   help="Min distance of negatives from positives (bp)")
    p.add_argument("--hard-neg-frac", type=float, default=0.5,
                   help="Fraction of negatives drawn from inside annotated gene "
                        "bodies (hard negatives). 0 = no hard negatives, 1 = all "
                        "hard negatives. Requires --gtf. Default: 0.5")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for negative sampling")
    p.add_argument("--threads", type=int, default=1,
                   help="Parallel workers for feature extraction")
    p.add_argument("--chrom-filter", default="",
                   help="Comma-separated list of chromosomes to EXCLUDE (e.g. chrM,chrUn)")
    p.add_argument("--gtf",
                   help="GENCODE GTF for gene-distance context features (TSS only; optional)")

    return p.parse_args()


# ---------------------------------------------------------------------------
# Database-specific loaders
# ---------------------------------------------------------------------------

def load_fantom_tss(bed_file: str, min_score: float = 0) -> pd.DataFrame:
    """Load FANTOM permissive CAGE peaks BED9 (hg38).

    Used as the primary positives source for the boundary classifier.
    Width labels are NOT assigned here — they come from the
    FANTOM5-hg38 fair-peaks spatial join in :func:`assign_width_labels`.

    IMPORTANT: this loader assumes the BED is hg38. The FANTOM5 permissive
    peak file as distributed (FANTOM_TSS_human.bed) is **hg19**; the build
    audit (qc/database_qc_report.md) confirmed every gene-anchored peak
    in the original file landed in hg19 windows. The hg38 version we use
    is produced by ``liftOver FANTOM_TSS_human.bed hg19ToHg38.over.chain
    FANTOM_TSS_human.hg38.bed``.

    Boundary coordinate: BED `start` for `+` strand, `end-1` for `-`
    strand. (`thickStart == start` for every row in this file, so there
    is no separate summit annotation to use.)
    """
    print(f"  Loading FANTOM TSS: {bed_file}")
    df = pd.read_csv(
        bed_file, sep="\t", header=None, skiprows=0,
        usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8],
        names=["chrom", "start", "end", "name", "score", "strand",
               "thickStart", "thickEnd", "rgb"],
        dtype={"chrom": str, "start": int, "end": int, "name": str,
               "score": str, "strand": str,
               "thickStart": int, "thickEnd": int, "rgb": str},
    )

    # Drop liftOver-pathology rows (peaks that ballooned to >500bp during
    # the hg19->hg38 lift; the 99th percentile of legitimate FANTOM peaks
    # is 59bp).
    n_before = len(df)
    df = df[(df["end"] - df["start"]).clip(lower=1) <= 500].reset_index(drop=True)
    if len(df) < n_before:
        print(f"    Dropped {n_before - len(df)} rows >500bp (liftOver pathologies)")

    # Parse confidence from name field "p1@GENE,0.1352"
    def parse_conf(name):
        if "," in str(name):
            try:
                return float(str(name).rsplit(",", 1)[1])
            except Exception:
                return 0.0
        return 0.0

    df["conf"] = df["name"].apply(parse_conf)
    if min_score > 0:
        df = df[df["conf"] >= min_score]

    # Strand-aware boundary coordinate
    df["pos"] = np.where(df["strand"] == "+", df["start"], df["end"] - 1).astype(int)
    df["score"] = (df["conf"] * 1000).astype(int)
    df["label"] = 1
    print(f"    → {len(df):,} sites (strand-aware BED edge)")
    return df[["chrom", "pos", "strand", "score", "label"]]


def load_cage_cluster_bed(
    bed_file: str,
    max_width: int = 5000,
) -> pd.DataFrame:
    """Load a generic merged CAGE cluster BED (hg38) for use as the
    width-target source.

    Currently used with the Lizio et al. 2017 ``hg38_fair_CAGE_peaks_phase1and2.bed``
    file (the FANTOM5 hg38-native robust set). Each row is a merged CAGE
    cluster with strand. The per-cluster width
    (`end - start`) is the regression target for the TSS width head: it
    represents the empirically-observed breadth of the initiation zone for
    each promoter.

    Filters out clusters wider than ``max_width`` bp to drop liftOver
    pathologies (one cluster lifted to 2.45 Mb in QC).

    Returns
    -------
    pd.DataFrame
        Columns: chrom, start, end, strand, width.
    """
    print(f"  Loading CAGE cluster BED: {bed_file}")
    df = pd.read_csv(
        bed_file, sep="\t", header=None,
        usecols=[0, 1, 2, 3, 4, 5],
        names=["chrom", "start", "end", "name", "score", "strand"],
        dtype={"chrom": str, "start": int, "end": int,
               "name": str, "score": str, "strand": str},
    )
    n_before = len(df)
    df["width"] = (df["end"].astype(int) - df["start"].astype(int)).clip(lower=1)
    df = df[df["width"] <= max_width].reset_index(drop=True)
    n_dropped = n_before - len(df)
    print(
        f"    → {len(df):,} clusters "
        f"({n_dropped} dropped as >{max_width}bp liftOver outliers); "
        f"width median={int(df['width'].median())}bp, "
        f"p90={int(df['width'].quantile(0.9))}bp, "
        f"p99={int(df['width'].quantile(0.99))}bp, "
        f"max={int(df['width'].max())}bp"
    )
    return df[["chrom", "start", "end", "strand", "width"]]


def assign_width_labels(
    positives: pd.DataFrame,
    cat_clusters: pd.DataFrame,
) -> pd.DataFrame:
    """Assign FANTOM-CAT cluster widths to TSS positives via spatial overlap.

    For each positive (chrom, pos, strand), find the FANTOM-CAT cluster on
    the same chrom+strand whose interval ``[start, end)`` contains pos.
    If multiple clusters overlap, take the narrowest (most specific). If
    none overlap, the width label is NaN — those rows stay positives for
    the classifier but are excluded from width-regressor training.

    The implementation uses ``np.searchsorted`` against per-(chrom, strand)
    sorted start/end arrays — O(N log M) total instead of the O(N·M)
    pairwise comparison.
    """
    print(f"  Assigning width labels via FANTOM-CAT spatial overlap...")
    pos = positives.reset_index(drop=True).copy()
    pos["width"] = np.nan

    # Build per-(chrom, strand) sorted arrays
    cat_index = {}
    for (chrom, strand), grp in cat_clusters.groupby(["chrom", "strand"]):
        order = np.argsort(grp["start"].values)
        cat_index[(chrom, strand)] = {
            "start": grp["start"].values[order],
            "end":   grp["end"].values[order],
            "width": grp["width"].values[order],
        }

    n_assigned = 0
    for (chrom, strand), grp in pos.groupby(["chrom", "strand"]):
        idx = cat_index.get((chrom, strand))
        if idx is None:
            continue
        positions = grp["pos"].values
        # candidate cluster i is one whose start <= pos
        cand = np.searchsorted(idx["start"], positions, side="right") - 1
        in_range = (cand >= 0) & (cand < len(idx["start"])) & (idx["end"][np.clip(cand, 0, len(idx["end"])-1)] > positions)
        # NB: this finds the rightmost cluster whose start<=pos.
        # If clusters are non-overlapping (the FANTOM-CAT case), this is
        # also the unique containing cluster.
        widths = np.where(in_range, idx["width"][np.clip(cand, 0, len(idx["width"])-1)], np.nan)
        pos.loc[grp.index, "width"] = widths
        n_assigned += int(np.isfinite(widths).sum())

    n_total = len(pos)
    print(
        f"    → {n_assigned:,} / {n_total:,} positives assigned a width label "
        f"({100*n_assigned/max(n_total,1):.1f}%); "
        f"the rest stay positives but are excluded from the width regressor."
    )
    return pos


def load_polyadb(bed_file: str, pas_file: Optional[str] = None) -> pd.DataFrame:
    """
    Load PolyA_DB BED6 file (hg38).

    If pas_file is provided, merge to get PAS signal type as quality weight.
    Prioritizes canonical AATAAA sites (score=1000) over NoPAS (score=200).
    """
    print(f"  Loading PolyA_DB: {bed_file}")
    df = pd.read_csv(bed_file, sep="\t", header=None,
                     names=["chrom", "start", "end", "name", "score", "strand"])

    # TTS position: end-1 for +, start for -
    df["pos"] = df.apply(
        lambda r: r["end"] - 1 if r["strand"] == "+" else r["start"], axis=1
    )
    df["label"] = 1

    if pas_file:
        print(f"  Merging PAS annotations from: {pas_file}")
        pas_df = pd.read_csv(pas_file, sep="\t")
        # PAS_ID matches BED name column
        pas_df = pas_df.rename(columns={
            "PAS_ID": "name",
            "PAS Signal": "pas_signal",
            "Intron/exon location": "intron_exon_loc",
            "PSE": "pse",
        })
        # PolyA_DB writes signal names with the RNA letter U (AAUAAA);
        # map to DNA spelling so the lookup actually matches.
        pas_df["pas_signal"] = (
            pas_df["pas_signal"].astype(str).str.replace("U", "T", regex=False)
        )
        signal_weights = {
            "AATAAA": 1000,
            "ATTAAA":  800,
            "Arich":   400,
            "OtherPAS": 300,
            "NoPAS":   200,
        }
        pas_df["pas_score"] = pas_df["pas_signal"].map(signal_weights).fillna(300)
        merge_cols = ["name", "pas_score"]
        for col in ["intron_exon_loc", "pse"]:
            if col in pas_df.columns:
                merge_cols.append(col)
        # Dedup PAS_ID before merging — duplicate IDs would multiply BED rows
        # and silently inflate the positive count.
        n_before = len(pas_df)
        pas_df = pas_df.drop_duplicates(subset="name", keep="first")
        if len(pas_df) < n_before:
            print(f"    Dropped {n_before - len(pas_df)} duplicate PAS_IDs before merge")
        df = df.merge(pas_df[merge_cols], on="name", how="left", validate="many_to_one")
        df["score"] = df["pas_score"].fillna(df["score"]).astype(int)

    print(f"    → {len(df):,} sites (1bp resolution; no width target)")
    # NOTE: intron_exon_loc and pse columns from polyAdb are NOT included
    # in the training feature set — they cause label leakage because they
    # are only non-zero for positive sites (from the database).
    # PolyA_DB is single-bp resolution so no `width` column is returned —
    # the width regressor only trains on TSS positives.
    return df[["chrom", "pos", "strand", "score", "label"]]


# ---------------------------------------------------------------------------
# Merging / deduplication
# ---------------------------------------------------------------------------

def merge_sites(dfs: list, window: int = 25) -> pd.DataFrame:
    """
    Merge sites from multiple sources, deduplicating within `window` bp.

    Sites within `window` bp on the same strand are merged (keep highest score).
    """
    if not dfs:
        raise ValueError("No DataFrames to merge")

    combined = pd.concat(dfs, ignore_index=True)
    combined = combined.sort_values(["chrom", "strand", "pos"])

    # Greedy merge: group into clusters within window
    merged = []
    for (chrom, strand), grp in combined.groupby(["chrom", "strand"]):
        grp = grp.sort_values("pos").reset_index(drop=True)
        cluster_start = None
        cluster_best = None
        for _, row in grp.iterrows():
            if cluster_start is None:
                cluster_start = row["pos"]
                cluster_best = row.copy()
            elif row["pos"] - cluster_start <= window:
                if row["score"] > cluster_best["score"]:
                    cluster_best = row.copy()
            else:
                merged.append(cluster_best)
                cluster_start = row["pos"]
                cluster_best = row.copy()
        if cluster_best is not None:
            merged.append(cluster_best)

    result = pd.DataFrame(merged).reset_index(drop=True)
    result["label"] = 1
    cols = ["chrom", "pos", "strand", "score", "label"]
    if "width" in result.columns:
        cols.append("width")
    return result[cols]


# ---------------------------------------------------------------------------
# Genomic context features (database-level, not sequence-level)
# ---------------------------------------------------------------------------

def compute_proximity_features(
    all_sites: pd.DataFrame,
    positives: pd.DataFrame,
) -> np.ndarray:
    """
    Compute inter-boundary proximity features for every site in all_sites.

    Uses the full merged positive reference set to characterize the local
    density of known boundaries around each candidate position.

    Features (3 total):
      - nn_dist_bp:        Log1p distance (bp) to nearest reference boundary
                           on the same strand. 0 for reference sites themselves;
                           large for isolated negatives.
      - n_within_500bp:    Count of reference boundaries within 500bp
                           (same strand). Captures promoter cluster density.
      - is_clustered:      1 if any reference boundary within 100bp, else 0.
                           Reflects whether this is part of a TSS/PAS cluster.

    These are computed from the reference database, so they are available
    both for positives (where they reflect true biology) and negatives
    (where they serve as a distance-to-nearest-real-boundary signal).

    Parameters
    ----------
    all_sites : pd.DataFrame
        Full site table (positives + negatives) with chrom, pos, strand.
    positives : pd.DataFrame
        Reference positive sites (chrom, pos, strand) — used as the
        ground-truth boundary set for distance computation.

    Returns
    -------
    np.ndarray, shape (n_sites, 1)
        Column: [nn_dist_log1p]
        Note: n_within_500bp and is_clustered are computed internally for
        logging but excluded from the returned array to prevent data leakage.
    """
    print("  Computing proximity features from reference boundaries...")

    # Build per-(chrom, strand) sorted position arrays from positives
    pos_index: dict = {}
    for (chrom, strand), grp in positives.groupby(["chrom", "strand"]):
        pos_index[(chrom, strand)] = np.sort(grp["pos"].values)

    nn_dists = np.zeros(len(all_sites), dtype=np.float32)
    n_within = np.zeros(len(all_sites), dtype=np.float32)
    is_clust = np.zeros(len(all_sites), dtype=np.float32)

    for i, row in enumerate(all_sites.itertuples(index=False)):
        key = (row.chrom, row.strand)
        ref_positions = pos_index.get(key)

        if ref_positions is None or len(ref_positions) == 0:
            nn_dists[i] = np.log1p(1_000_000)  # no reference on this chrom/strand
            continue

        dists = np.abs(ref_positions - row.pos)

        # Nearest neighbor (exclude self — distance 0 means it IS a reference site)
        nonzero = dists[dists > 0]
        nn_dists[i] = np.log1p(nonzero.min()) if len(nonzero) > 0 else 0.0
        n_within[i] = float((dists <= 500).sum())
        is_clust[i] = 1.0 if (dists <= 100).any() else 0.0

    # Return only the continuous distance signal — the binary cluster features
    # (n_within_500bp, is_clustered) are excluded from training features to
    # avoid data leakage: negatives are placed away from positives, so those
    # features would be trivially 0 for negatives regardless of site count.
    # At inference time, cluster density comes from the HDBSCAN read clusters.
    print(f"    nn_dist_log1p: mean={nn_dists.mean():.2f}  "
          f"(n_within_500bp={n_within.mean():.1f}, is_clustered={is_clust.mean()*100:.1f}% — logged only)")
    return nn_dists.reshape(-1, 1)


PROXIMITY_FEATURE_NAMES = [
    "ctx_nn_dist_log1p",     # log1p(bp to nearest reference boundary, same strand)
    # ctx_n_within_500bp and ctx_is_clustered are computed but excluded from
    # training to avoid leakage; used only at inference from HDBSCAN clusters.
]


def compute_genomic_location_features(
    all_sites: pd.DataFrame,
    boundary_type: str,
    gtf_file: Optional[str] = None,
) -> np.ndarray:
    """
    Compute gene-body context features: distance to nearest annotated gene.

    Both TSS and TTS use the same GTF-derived features:
      - gene_dist_log1p:   log1p(bp to nearest gene body start on same strand)
      - is_intergenic:     1 if >1kb from any annotated gene
      - is_antisense:      1 if nearest gene is on opposite strand
      - gene_body_frac:    relative position within nearest gene (0=start, 1=end)

    NOTE: For TTS, we previously used polyAdb annotation columns
    (intron_exon_loc, pse), but these caused *label leakage* because they
    are only non-zero for positive sites from the database and are trivially
    0 for randomly generated negatives, enabling perfect class separation.

    For sites where annotation is unavailable, defaults to 0 (neutral).

    Parameters
    ----------
    all_sites : pd.DataFrame
        Full site table with chrom, pos, strand.
    boundary_type : str
        'tss' or 'tts'
    gtf_file : str or None
        Path to GENCODE GTF (used for gene-distance features).
        If None, returns zeros for all gene-distance features.

    Returns
    -------
    np.ndarray, shape (n_sites, 4)
    """
    n = len(all_sites)

    # Both TSS and TTS use GTF-derived gene-distance features.
    # NOTE: For TTS, we previously used polyAdb annotation columns
    # (intron_exon_loc, pse) but these caused label leakage since they
    # are only available for positive sites from the database and are
    # trivially 0 for negatives, creating perfect class separation.
    gene_dist = np.full(n, np.log1p(1_000_000), dtype=np.float32)
    is_intergenic = np.ones(n, dtype=np.float32)
    is_antisense = np.zeros(n, dtype=np.float32)
    gene_body_frac = np.full(n, 0.5, dtype=np.float32)

    if gtf_file and Path(gtf_file).exists():
        print(f"  Loading gene coordinates from GTF for {boundary_type.upper()} context features...")
        # Parse gene-level entries only (fast)
        genes = []
        with open(gtf_file) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip().split("\t")
                if len(parts) < 9 or parts[2] != "gene":
                    continue
                genes.append({
                    "chrom": parts[0],
                    "start": int(parts[3]),
                    "end":   int(parts[4]),
                    "strand": parts[6],
                })
        gene_df = pd.DataFrame(genes)
        print(f"    Loaded {len(gene_df):,} gene records")

        # For each site, find nearest gene on same chromosome.
        # NB: use a scalar groupby key (not a 1-element list) so the loop
        # variable receives the chrom directly on both pandas 1.x (py3.10)
        # and pandas 2.x (py3.11+). The list form returns a 1-tuple key
        # only on the newer version.
        for chrom, grp_sites in all_sites.groupby("chrom"):
            chrom_genes = gene_df[gene_df["chrom"] == chrom]
            if chrom_genes.empty:
                continue
            g_starts = chrom_genes["start"].values
            g_ends = chrom_genes["end"].values
            g_strands = chrom_genes["strand"].values

            for idx in grp_sites.index:
                pos = all_sites.at[idx, "pos"]
                site_strand = all_sites.at[idx, "strand"]

                # Distance: 0 if inside gene, else gap to nearest end
                inside = (g_starts <= pos) & (pos <= g_ends)
                if inside.any():
                    dist_arr = np.zeros(inside.sum())
                    nearest_idx = np.where(inside)[0][0]
                else:
                    dist_arr = np.minimum(
                        np.abs(g_starts - pos), np.abs(g_ends - pos)
                    )
                    nearest_idx = np.argmin(dist_arr)
                    dist_arr = dist_arr[[nearest_idx]]

                nearest_dist = dist_arr[0] if len(dist_arr) > 0 else 1_000_000
                gene_dist[idx] = np.log1p(nearest_dist)
                is_intergenic[idx] = 1.0 if nearest_dist > 1000 else 0.0
                is_antisense[idx] = 1.0 if g_strands[nearest_idx] != site_strand else 0.0
                # Position within gene body (0 = TSS, 1 = TES for same-strand gene)
                g_len = max(g_ends[nearest_idx] - g_starts[nearest_idx], 1)
                if site_strand == "+":
                    frac = (pos - g_starts[nearest_idx]) / g_len
                else:
                    frac = (g_ends[nearest_idx] - pos) / g_len
                gene_body_frac[idx] = float(np.clip(frac, 0, 1))
    else:
        if gtf_file:
            print(f"  WARNING: GTF not found at {gtf_file}, skipping gene-distance features")
        else:
            print("  No GTF provided, gene-distance features will be 0")

    return np.column_stack([gene_dist, is_intergenic, is_antisense, gene_body_frac])


GENOMIC_LOCATION_FEATURE_NAMES = {
    "tss": [
        "ctx_gene_dist_log1p",      # log1p(bp to nearest gene body)
        "ctx_is_intergenic",        # 1 if >1kb from any gene
        "ctx_is_antisense",         # 1 if nearest gene is on opposite strand
        "ctx_gene_body_frac",       # relative position within nearest gene
    ],
    "tts": [
        "ctx_gene_dist_log1p",      # log1p(bp to nearest gene body)
        "ctx_is_intergenic",        # 1 if >1kb from any gene
        "ctx_is_antisense",         # 1 if nearest gene is on opposite strand
        "ctx_gene_body_frac",       # relative position within nearest gene
    ],
}


# ---------------------------------------------------------------------------
# Feature extraction (with optional parallelism)
# ---------------------------------------------------------------------------

def _extract_worker(args):
    """Module-level worker for multiprocessing (must be picklable)."""
    chunk, genome_fasta, boundary_type = args
    ext = SequenceFeatureExtractor(genome_fasta, boundary_type)
    return ext.extract_batch(chunk)


def extract_features_parallel(
    df: pd.DataFrame,
    genome_fasta: str,
    boundary_type: str,
    n_workers: int = 1,
) -> np.ndarray:
    """Extract features, optionally splitting across workers."""
    if n_workers <= 1:
        extractor = SequenceFeatureExtractor(genome_fasta, boundary_type)
        print(f"  Extracting features for {len(df):,} sites (single-threaded)...")
        X = extractor.extract_batch(df)
    else:
        from multiprocessing import Pool

        chunk_size = max(1, len(df) // n_workers)
        chunks = [df.iloc[i:i+chunk_size] for i in range(0, len(df), chunk_size)]
        # Pass genome_fasta and boundary_type as part of the argument tuple
        # so the worker function is a plain module-level callable (picklable)
        work_items = [(chunk, genome_fasta, boundary_type) for chunk in chunks]

        print(f"  Extracting features across {n_workers} workers ({len(df):,} sites)...")
        with Pool(n_workers) as pool:
            parts = pool.map(_extract_worker, work_items)
        X = np.vstack(parts)

    return X


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    print("=" * 70)
    print(f"flAIr: Preparing {args.type.upper()} training data")
    print("=" * 70)
    print()

    # Excluded chromosomes
    excluded_chroms = set(args.chrom_filter.split(",")) if args.chrom_filter else {"chrM"}
    print(f"Excluding chromosomes: {excluded_chroms}")

    # --- Load positive sites ---
    print("\n[1/4] Loading positive sites from reference databases...")
    positive_dfs = []

    if args.type == "tss":
        if args.fantom:
            df_f = load_fantom_tss(args.fantom, args.min_score)
            positive_dfs.append(df_f)
        if not positive_dfs:
            print("ERROR: No TSS input files provided (--fantom required)")
            sys.exit(1)
    else:  # tts
        if args.polyadb:
            df_p = load_polyadb(args.polyadb, args.polyadb_pas)
            positive_dfs.append(df_p)
        if not positive_dfs:
            print("ERROR: No TTS input files provided (--polyadb required)")
            sys.exit(1)

    # Merge databases
    if len(positive_dfs) > 1:
        print(f"\n  Merging {len(positive_dfs)} databases (window={args.merge_window}bp)...")
        positives = merge_sites(positive_dfs, window=args.merge_window)
    else:
        positives = positive_dfs[0]

    # Filter chromosomes
    positives = positives[~positives["chrom"].isin(excluded_chroms)]
    print(f"  Total positive sites: {len(positives):,}")

    # Subsample if requested
    if args.max_sites and len(positives) > args.max_sites:
        positives = positives.sample(args.max_sites, random_state=args.seed)
        print(f"  Subsampled to {len(positives):,} sites")

    # --- Assign width labels via FANTOM5-hg38 fair peaks (TSS only) ---
    if args.type == "tss" and getattr(args, "fantom_fair", None):
        if not Path(args.fantom_fair).exists():
            print(f"  WARNING: --fantom-fair path does not exist: {args.fantom_fair}; "
                  f"width regressor will have no labels")
            positives["width"] = np.nan
        else:
            fair_clusters = load_cage_cluster_bed(
                args.fantom_fair, max_width=args.fantom_fair_max_width,
            )
            positives = assign_width_labels(positives, fair_clusters)
    else:
        # TTS has no width target (PolyA_DB is single-bp)
        positives["width"] = np.nan

    # NOTE on negative generation: negatives are placed ≥neg_min_dist from any
    # positive. This means ctx_is_clustered and ctx_n_within_500bp will be
    # trivially 0 for all negatives if neg_min_dist >= 500bp (the default).
    # To prevent this from becoming a perfect-separation leak, negatives need
    # neg_min_dist < 500bp OR these features should be excluded from sequence
    # training and only used at inference. We handle this by setting
    # neg_min_dist to 150bp when proximity features are requested, so some
    # negatives fall within 500bp of a real boundary.

    # --- Generate negative sites ---
    # Use a reduced min_distance so proximity features (n_within_500bp) are
    # not trivially 0 for all negatives. 150bp keeps negatives clearly
    # separated from positives while allowing some to land near real boundaries.
    effective_neg_min_dist = min(args.neg_min_dist, 150)
    n_total_neg = args.neg_ratio * len(positives)

    use_hard = args.hard_neg_frac > 0 and getattr(args, "gtf", None) and Path(args.gtf).exists()
    if use_hard:
        n_hard = int(round(n_total_neg * args.hard_neg_frac))
        n_easy = n_total_neg - n_hard
    else:
        if args.hard_neg_frac > 0:
            print(
                f"  NOTE: --hard-neg-frac={args.hard_neg_frac} requires --gtf; "
                f"falling back to all easy negatives"
            )
        n_hard, n_easy = 0, n_total_neg

    print(f"\n[2/4] Generating {n_total_neg:,} negatives "
          f"({n_easy:,} random + {n_hard:,} hard, "
          f"min_dist={effective_neg_min_dist}bp)...")

    neg_dfs = []
    if n_easy > 0:
        df_easy = generate_negative_sites(
            positive_df=positives,
            genome_fasta=args.genome,
            n_negatives=n_easy,
            min_distance=effective_neg_min_dist,
            seed=args.seed,
        )
        df_easy = df_easy[~df_easy["chrom"].isin(excluded_chroms)]
        neg_dfs.append(df_easy)
        print(f"  Random negatives:  {len(df_easy):,}")

    if n_hard > 0:
        df_hard = generate_hard_negatives(
            positive_df=positives,
            gtf_file=args.gtf,
            n_negatives=n_hard,
            boundary_type=args.type,
            min_distance=effective_neg_min_dist,
            seed=args.seed + 1,
        )
        df_hard = df_hard[~df_hard["chrom"].isin(excluded_chroms)]
        neg_dfs.append(df_hard)
        print(f"  Hard (gene-body) negatives: {len(df_hard):,}")

    negatives = pd.concat(neg_dfs, ignore_index=True) if neg_dfs else pd.DataFrame()
    print(f"  Total negatives: {len(negatives):,}")

    # Combine
    all_sites = pd.concat([positives, negatives], ignore_index=True)
    all_sites = all_sites.sample(frac=1, random_state=args.seed).reset_index(drop=True)
    print(f"  Combined: {len(all_sites):,} total sites ({(all_sites['label']==1).sum():,} pos, {(all_sites['label']==0).sum():,} neg)")

    # --- Extract sequence features ---
    print(f"\n[3/4] Extracting sequence features (window: {args.type.upper()})...")
    X_seq = extract_features_parallel(
        all_sites, args.genome, args.type, n_workers=args.threads
    )

    # Remove rows where sequence extraction failed (out-of-bounds / missing chrom)
    valid_mask = ~np.isnan(X_seq).any(axis=1)
    X_seq = X_seq[valid_mask]
    all_sites_valid = all_sites.iloc[valid_mask].reset_index(drop=True)

    failed = (~valid_mask).sum()
    if failed > 0:
        print(f"  WARNING: {failed} sites failed extraction (out-of-bounds or missing chrom)")

    # --- Compute genomic context features ---
    print(f"\n  Computing genomic context features...")

    # 1. Proximity features (database-level nearest-neighbor)
    X_prox = compute_proximity_features(all_sites_valid, positives)

    # 2. Genomic location features (gene body / exon-intron / PAS annotation)
    X_loc = compute_genomic_location_features(
        all_sites_valid,
        boundary_type=args.type,
        gtf_file=getattr(args, "gtf", None),
    )

    # Combine all features
    X = np.hstack([X_seq, X_prox, X_loc])
    print(f"  Feature matrix: {X.shape} "
          f"(seq={X_seq.shape[1]}, proximity={X_prox.shape[1]}, location={X_loc.shape[1]})")

    y = all_sites_valid["label"].values
    scores = all_sites_valid["score"].values
    chroms = all_sites_valid["chrom"].values
    positions = all_sites_valid["pos"].values
    strands = all_sites_valid["strand"].values
    # Width is the regression target for the width head. NaN for negatives
    # (they have no defined initiation/termination zone width). Stored
    # separately from y so the regressor can be trained on positives only.
    if "width" in all_sites_valid.columns:
        widths = all_sites_valid["width"].astype(float).fillna(-1).values
    else:
        widths = np.full(len(all_sites_valid), -1.0, dtype=float)

    # --- Get feature names ---
    extractor = SequenceFeatureExtractor(args.genome, args.type)
    feature_names = (
        extractor.get_feature_names()
        + PROXIMITY_FEATURE_NAMES
        + GENOMIC_LOCATION_FEATURE_NAMES[args.type]
    )

    # --- Save ---
    print(f"\n[4/4] Saving to {args.output}...")
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        str(out_path),
        X=X,
        y=y,
        scores=scores,
        chroms=chroms,
        positions=positions,
        strands=strands,
        widths=widths,
        feature_names=np.array(feature_names),
        boundary_type=np.array([args.type]),
    )

    print(f"  Saved {X.shape[0]} samples × {X.shape[1]} features")
    print(f"  Class balance: {(y==1).sum()} positives, {(y==0).sum()} negatives")
    print("\nDone!")


if __name__ == "__main__":
    main()
 