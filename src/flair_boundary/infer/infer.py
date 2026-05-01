#!/usr/bin/env python3
"""
Hybrid boundary scoring and selection for FLAIR transcriptome isoforms.

This module implements the inference side of the hybrid boundary model:
  1. Sequence classifier score  — from trained XGBoost model
  2. Read-level evidence score  — from BAM file at candidate boundaries
  3. Composite score            — product of the above with artifact penalty
  4. Boundary selection         — highest-scoring candidate per splice chain

The composite score formula (per candidate boundary b):
  score(b) = seq_conf(b) * read_support_norm(b) * (1 - artifact_penalty(b))

Where:
  seq_conf:         XGBoost probability output P(boundary | sequence)
  read_support_norm: log1p(read_count) normalized within cluster, scaled by density rank
  artifact_penalty: soft-clipping fraction (TSS) or A-richness penalty (TTS)

Usage:
    python infer.py \\
        --candidates candidates.bed \\
        --bam aligned_reads.bam \\
        --genome GRCh38.primary_assembly.genome.fa \\
        --tss-model model_weights/tss/xgboost.pkl \\
        --tts-model model_weights/tts/xgboost.pkl \\
        --output scored_boundaries.tsv \\
        --type tss

    # As a library:
    from flair_boundary.infer.infer import HybridBoundaryScorer
    scorer = HybridBoundaryScorer(
        genome_fasta="GRCh38.fa",
        tss_model_path="model_weights/tss/xgboost.pkl",
        tts_model_path="model_weights/tts/xgboost.pkl",
    )
    scored = scorer.score_candidates(candidates_df, bam_path="aligned.bam")
    selected = scorer.select_boundaries(scored)
"""

import argparse
import pickle
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from flair_boundary.extract.features import SequenceFeatureExtractor


# ---------------------------------------------------------------------------
# Read-level feature computation
# ---------------------------------------------------------------------------

class ReadLevelFeatures:
    """
    Compute read-level evidence features at candidate boundary positions.

    Features extracted:
      - n_reads:          Total reads supporting this boundary
      - soft_clip_frac:   Fraction of reads with soft-clip at boundary (5' end artifact indicator)
      - a_rich_frac:      Fraction of reads with A-rich tails downstream (3' artifact indicator)
      - density_rank:     Normalized rank of this boundary within its cluster (0=lowest, 1=highest)
      - dispersion:       Boundary position dispersion (std of positions within cluster, normalized)
    """

    def __init__(self, bam_path: str, min_mapq: int = 10):
        self.bam = pysam.AlignmentFile(bam_path, "rb")
        self.min_mapq = min_mapq

    def close(self):
        self.bam.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def compute_at_position(
        self,
        chrom: str,
        pos: int,
        strand: str,
        boundary_type: str,
        window: int = 50,
    ) -> Dict[str, float]:
        """
        Compute read-level features at a single boundary position.

        Parameters
        ----------
        chrom : str
        pos : int
            0-based boundary coordinate.
        strand : str
            '+' or '-'
        boundary_type : str
            'tss' or 'tts'
        window : int
            Window around position to fetch reads.
        """
        start = max(0, pos - window)
        end = pos + window

        reads = list(self.bam.fetch(chrom, start, end))
        reads = [r for r in reads if not r.is_unmapped and r.mapping_quality >= self.min_mapq]

        # Strand filter
        if strand == "+":
            reads = [r for r in reads if not r.is_reverse]
        else:
            reads = [r for r in reads if r.is_reverse]

        if not reads:
            return {
                "n_reads": 0,
                "soft_clip_frac": 0.0,
                "a_rich_frac": 0.0,
                "boundary_dispersion": 0.0,
            }

        n_reads = len(reads)

        # Soft-clipping fraction at the boundary end
        sc_count = 0
        boundary_positions = []
        a_rich_count = 0

        for read in reads:
            cigar = read.cigartuples
            if cigar is None:
                continue

            if boundary_type == "tss":
                # For TSS: check soft-clipping at 5' end of read
                # 5' end = left end for + strand, right end for - strand
                if strand == "+":
                    if cigar[0][0] == 4:  # soft clip at start
                        sc_count += 1
                    boundary_positions.append(read.reference_start)
                else:
                    if cigar[-1][0] == 4:  # soft clip at end
                        sc_count += 1
                    boundary_positions.append(read.reference_end)

            else:  # tts
                # For TTS: check soft-clipping at 3' end and A-richness
                if strand == "+":
                    if cigar[-1][0] == 4:  # soft clip at end
                        sc_count += 1
                    boundary_positions.append(read.reference_end)
                    # Check A-richness in soft-clipped sequence
                    seq = read.query_sequence
                    if seq and cigar[-1][0] == 4:
                        clip_len = cigar[-1][1]
                        clip_seq = seq[-clip_len:]
                        if len(clip_seq) >= 5:
                            a_frac = clip_seq.count("A") / len(clip_seq)
                            if a_frac >= 0.7:
                                a_rich_count += 1
                else:
                    if cigar[0][0] == 4:  # soft clip at start
                        sc_count += 1
                    boundary_positions.append(read.reference_start)
                    seq = read.query_sequence
                    if seq and cigar[0][0] == 4:
                        clip_len = cigar[0][1]
                        clip_seq = seq[:clip_len]
                        if len(clip_seq) >= 5:
                            # Complement A-richness for minus strand = T in clip
                            t_frac = clip_seq.count("T") / len(clip_seq)
                            if t_frac >= 0.7:
                                a_rich_count += 1

        sc_frac = sc_count / max(n_reads, 1)
        a_rich_frac = a_rich_count / max(n_reads, 1)

        # Boundary position dispersion (spread of read endpoints near this boundary)
        if len(boundary_positions) > 1:
            bp_arr = np.array(boundary_positions)
            # Keep only those within 2*window of our target
            bp_arr = bp_arr[np.abs(bp_arr - pos) <= window]
            dispersion = float(np.std(bp_arr)) if len(bp_arr) > 1 else 0.0
        else:
            dispersion = 0.0

        return {
            "n_reads": n_reads,
            "soft_clip_frac": sc_frac,
            "a_rich_frac": a_rich_frac,
            "boundary_dispersion": dispersion,
        }

    def compute_batch(
        self,
        df: pd.DataFrame,
        boundary_type: str,
        window: int = 50,
    ) -> pd.DataFrame:
        """
        Compute read-level features for all positions in df.

        df must have columns: chrom, pos, strand
        Returns df with additional feature columns.
        """
        results = []
        for _, row in df.iterrows():
            feats = self.compute_at_position(
                row["chrom"], int(row["pos"]), row["strand"],
                boundary_type, window
            )
            results.append(feats)

        feat_df = pd.DataFrame(results, index=df.index)
        return pd.concat([df, feat_df], axis=1)


# ---------------------------------------------------------------------------
# Composite scoring
# ---------------------------------------------------------------------------

def compute_read_support_score(
    n_reads: np.ndarray,
    density_ranks: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Normalize read counts to [0, 1] using log1p scaling and optional density rank.

    score = log1p(n_reads) / max(log1p(n_reads))
    If density_ranks provided: score = 0.7*log_norm + 0.3*density_rank_norm
    """
    log_counts = np.log1p(n_reads.astype(float))
    max_log = log_counts.max()
    if max_log == 0:
        norm = np.zeros_like(log_counts)
    else:
        norm = log_counts / max_log

    if density_ranks is not None:
        max_rank = density_ranks.max()
        rank_norm = density_ranks / max(max_rank, 1)
        norm = 0.7 * norm + 0.3 * rank_norm

    return norm


def compute_artifact_penalty(
    df: pd.DataFrame,
    boundary_type: str,
) -> np.ndarray:
    """
    Compute artifact penalty score for each candidate.

    TSS: soft-clipping fraction (high soft-clip → likely internal fragment)
    TTS: A-rich fraction (high A-richness → likely internal priming artifact)

    Penalty is in [0, 1]; higher = more likely artifact.
    Final factor applied: (1 - penalty)
    """
    if boundary_type == "tss":
        # Soft-clip at TSS is actually normal (reads start exactly at TSS)
        # High soft-clip > 0.7 of reads suggests reads don't actually start there
        # but are clipped internal fragments — penalize
        sc = df["soft_clip_frac"].values if "soft_clip_frac" in df.columns else np.zeros(len(df))
        # Moderate clipping (0.2-0.5) is normal at TSS; penalize if >0.7
        penalty = np.clip((sc - 0.5) / 0.5, 0, 1)
    else:  # tts
        # A-richness in clipped reads strongly suggests internal priming
        a_rich = df["a_rich_frac"].values if "a_rich_frac" in df.columns else np.zeros(len(df))
        # Dispersion also matters: high dispersion at TTS suggests fuzzy termination
        disp = df["boundary_dispersion"].values if "boundary_dispersion" in df.columns else np.zeros(len(df))
        disp_penalty = np.clip(disp / 100.0, 0, 1)  # normalize: 100bp dispersion = max penalty
        penalty = 0.7 * a_rich + 0.3 * disp_penalty

    return np.clip(penalty, 0, 1)


def composite_score(
    seq_conf: np.ndarray,
    read_support: np.ndarray,
    artifact_penalty: np.ndarray,
    seq_weight: float = 0.5,
    read_weight: float = 0.5,
) -> np.ndarray:
    """
    Compute composite boundary score.

    score = (seq_weight * seq_conf + read_weight * read_support) * (1 - artifact_penalty)

    Parameters are tunable; defaults weight sequence and read support equally,
    then apply artifact penalty multiplicatively.
    """
    base = seq_weight * seq_conf + read_weight * read_support
    return base * (1.0 - artifact_penalty)


# ---------------------------------------------------------------------------
# Main scorer class
# ---------------------------------------------------------------------------

class HybridBoundaryScorer:
    """
    Complete hybrid boundary scoring pipeline.

    Loads trained sequence classifiers for TSS and TTS, then scores
    candidate boundary positions using both sequence and read-level evidence.

    Parameters
    ----------
    genome_fasta : str
        Path to GRCh38 genome FASTA.
    tss_model_path : str or None
        Path to trained TSS classifier .pkl.
    tts_model_path : str or None
        Path to trained TTS classifier .pkl.
    seq_weight : float
        Weight for sequence classifier score (0-1).
    read_weight : float
        Weight for read-level evidence (0 = sequence only).
    """

    def __init__(
        self,
        genome_fasta: str,
        tss_model_path: Optional[str] = None,
        tts_model_path: Optional[str] = None,
        seq_weight: float = 0.5,
        read_weight: float = 0.5,
    ):
        self.genome_fasta = genome_fasta
        self.seq_weight = seq_weight
        self.read_weight = read_weight

        # Load models
        self.models: Dict[str, object] = {}
        self.extractors: Dict[str, SequenceFeatureExtractor] = {}

        for btype, path in [("tss", tss_model_path), ("tts", tts_model_path)]:
            if path and Path(path).exists():
                print(f"Loading {btype.upper()} model from {path}...")
                with open(path, "rb") as f:
                    self.models[btype] = pickle.load(f)
                self.extractors[btype] = SequenceFeatureExtractor(genome_fasta, btype)
                print(f"  Loaded.")
            else:
                if path:
                    print(f"WARNING: {btype.upper()} model not found at {path}")

    def score_candidates(
        self,
        candidates: pd.DataFrame,
        boundary_type: str,
        bam_path: Optional[str] = None,
        density_ranks: Optional[np.ndarray] = None,
    ) -> pd.DataFrame:
        """
        Score a set of candidate boundaries.

        Parameters
        ----------
        candidates : pd.DataFrame
            Must have columns: chrom, pos, strand.
            Optionally: cluster_id (for grouping candidates).
        boundary_type : str
            'tss' or 'tts'
        bam_path : str or None
            BAM file for read-level features. If None, read score = 0.
        density_ranks : np.ndarray or None
            Pre-computed density rank per candidate (from HDBSCAN cluster).
            If None, uses only read count for read support.

        Returns
        -------
        pd.DataFrame
            candidates + columns: seq_conf, read_support, artifact_penalty, composite_score
        """
        df = candidates.copy().reset_index(drop=True)
        n = len(df)

        if n == 0:
            return df

        # 1. Sequence classifier score
        if boundary_type in self.models:
            extractor = self.extractors[boundary_type]
            print(f"  Extracting sequence features for {n} candidates...")
            X = extractor.extract_batch(df)
            X = np.nan_to_num(X, nan=0.0)
            clf = self.models[boundary_type]
            seq_conf = clf.predict_proba(X)[:, 1]
        else:
            print(f"  WARNING: No {boundary_type} model loaded, using seq_conf=0.5")
            seq_conf = np.full(n, 0.5)

        df["seq_conf"] = seq_conf

        # 2. Read-level features
        if bam_path and self.read_weight > 0:
            print(f"  Computing read-level features from {bam_path}...")
            with ReadLevelFeatures(bam_path) as rfe:
                df = rfe.compute_batch(df, boundary_type)
            n_reads = df["n_reads"].values.astype(float)
            artifact_pen = compute_artifact_penalty(df, boundary_type)
        else:
            df["n_reads"] = 0
            df["soft_clip_frac"] = 0.0
            df["a_rich_frac"] = 0.0
            df["boundary_dispersion"] = 0.0
            n_reads = np.zeros(n)
            artifact_pen = np.zeros(n)

        # 3. Read support normalization (within-cluster if cluster_id present)
        if "cluster_id" in df.columns:
            read_support = np.zeros(n)
            for cid, grp in df.groupby("cluster_id"):
                idx = grp.index
                dr = density_ranks[idx] if density_ranks is not None else None
                read_support[idx] = compute_read_support_score(
                    n_reads[idx], dr
                )
        else:
            read_support = compute_read_support_score(n_reads, density_ranks)

        df["read_support"] = read_support
        df["artifact_penalty"] = artifact_pen

        # 4. Composite score
        df["composite_score"] = composite_score(
            seq_conf=seq_conf,
            read_support=read_support,
            artifact_penalty=artifact_pen,
            seq_weight=self.seq_weight,
            read_weight=self.read_weight,
        )

        return df

    def select_boundaries(
        self,
        scored_df: pd.DataFrame,
        group_col: str = "splice_chain",
        min_score: float = 0.3,
        multi_boundary_threshold: float = 0.5,
        multi_boundary_min_dist: int = 50,
    ) -> pd.DataFrame:
        """
        Select the best boundary per splice junction chain.

        For each splice chain, selects the highest-scoring candidate. If multiple
        candidates exceed `multi_boundary_threshold` and are >min_dist apart,
        they are both retained (representing distinct isoforms).

        Parameters
        ----------
        scored_df : pd.DataFrame
            Output from score_candidates().
        group_col : str
            Column to group candidates by (e.g. 'splice_chain').
        min_score : float
            Minimum composite score to consider a candidate at all.
        multi_boundary_threshold : float
            If a secondary candidate exceeds this and is far enough from primary,
            report it as a distinct isoform.
        multi_boundary_min_dist : int
            Minimum distance (bp) for two boundaries to be called distinct isoforms.

        Returns
        -------
        pd.DataFrame
            Selected boundaries with is_primary and is_alternative columns.
        """
        if group_col not in scored_df.columns:
            # No grouping — just filter by score and return top hit
            passing = scored_df[scored_df["composite_score"] >= min_score].copy()
            passing["is_primary"] = True
            passing["is_alternative"] = False
            return passing.sort_values("composite_score", ascending=False)

        selected = []
        for chain_id, grp in scored_df.groupby(group_col):
            grp = grp.sort_values("composite_score", ascending=False).reset_index(drop=True)

            # Primary: highest score
            primary = grp.iloc[0].copy()
            if primary["composite_score"] < min_score:
                continue

            primary["is_primary"] = True
            primary["is_alternative"] = False
            primary["splice_chain_id"] = chain_id
            selected.append(primary)

            # Secondary candidates (alternative isoforms)
            for _, candidate in grp.iloc[1:].iterrows():
                if candidate["composite_score"] < multi_boundary_threshold:
                    break  # sorted descending, no point checking further
                dist = abs(candidate["pos"] - primary["pos"])
                if dist >= multi_boundary_min_dist:
                    alt = candidate.copy()
                    alt["is_primary"] = False
                    alt["is_alternative"] = True
                    alt["splice_chain_id"] = chain_id
                    selected.append(alt)

        if not selected:
            return pd.DataFrame()

        result = pd.DataFrame(selected)
        return result.sort_values(
            ["splice_chain_id", "is_primary", "composite_score"],
            ascending=[True, False, False]
        )


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Score and select hybrid boundaries for FLAIR transcriptome",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--candidates", required=True,
                   help="BED file of candidate boundary positions "
                        "(chrom, start, end, name, score, strand, [chain_id])")
    p.add_argument("--genome", required=True,
                   help="GRCh38 genome FASTA")
    p.add_argument("--output", required=True,
                   help="Output TSV file for scored/selected boundaries")
    p.add_argument("--type", required=True, choices=["tss", "tts"],
                   help="Boundary type to score")

    p.add_argument("--tss-model",
                   help="Path to trained TSS model .pkl")
    p.add_argument("--tts-model",
                   help="Path to trained TTS model .pkl")
    p.add_argument("--bam",
                   help="BAM file for read-level features (optional)")

    p.add_argument("--seq-weight", type=float, default=0.5,
                   help="Weight for sequence classifier score")
    p.add_argument("--read-weight", type=float, default=0.5,
                   help="Weight for read-level evidence score")
    p.add_argument("--min-score", type=float, default=0.3,
                   help="Minimum composite score to report")
    p.add_argument("--multi-threshold", type=float, default=0.5,
                   help="Threshold for calling alternative isoform")
    p.add_argument("--multi-min-dist", type=int, default=50,
                   help="Min distance (bp) for alternative isoform call")
    p.add_argument("--select", action="store_true",
                   help="Run boundary selection (requires splice_chain column)")

    return p.parse_args()


def main():
    args = parse_args()

    print("=" * 70)
    print(f"flAIr: Hybrid Boundary Scoring ({args.type.upper()})")
    print("=" * 70)
    print()

    # Load candidates
    print(f"Loading candidates from {args.candidates}...")
    col_names = ["chrom", "start", "end", "name", "score", "strand"]
    try:
        cands = pd.read_csv(args.candidates, sep="\t", header=None)
        if cands.shape[1] >= 7:
            cands.columns = col_names + [f"col{i}" for i in range(7, cands.shape[1]+1)]
            cands = cands.rename(columns={"col7": "splice_chain"})
        else:
            cands.columns = col_names[:cands.shape[1]]
    except Exception as e:
        print(f"ERROR reading candidates: {e}")
        sys.exit(1)

    # Set pos as the boundary coordinate
    if args.type == "tss":
        cands["pos"] = cands.apply(
            lambda r: r["start"] if r["strand"] == "+" else r["end"] - 1, axis=1
        )
    else:
        cands["pos"] = cands.apply(
            lambda r: r["end"] - 1 if r["strand"] == "+" else r["start"], axis=1
        )

    print(f"  Loaded {len(cands):,} candidates")

    # Build scorer
    scorer = HybridBoundaryScorer(
        genome_fasta=args.genome,
        tss_model_path=args.tss_model,
        tts_model_path=args.tts_model,
        seq_weight=args.seq_weight,
        read_weight=args.read_weight,
    )

    # Score
    print(f"\nScoring candidates...")
    scored = scorer.score_candidates(
        candidates=cands,
        boundary_type=args.type,
        bam_path=args.bam,
    )

    # Select
    if args.select:
        print(f"\nSelecting boundaries (min_score={args.min_score})...")
        result = scorer.select_boundaries(
            scored,
            group_col="splice_chain" if "splice_chain" in scored.columns else "chrom",
            min_score=args.min_score,
            multi_boundary_threshold=args.multi_threshold,
            multi_boundary_min_dist=args.multi_min_dist,
        )
        print(f"  Selected {len(result):,} boundaries "
              f"({result.get('is_primary', pd.Series(dtype=bool)).sum()} primary, "
              f"{result.get('is_alternative', pd.Series(dtype=bool)).sum()} alternative)")
    else:
        result = scored[scored["composite_score"] >= args.min_score].copy()
        print(f"  Passing score threshold: {len(result):,} / {len(scored):,}")

    # Save
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out_path, sep="\t", index=False)
    print(f"\nSaved scored boundaries: {out_path}")

    # Summary stats
    if len(result) > 0:
        print(f"\nScore distribution:")
        print(f"  seq_conf:       mean={result['seq_conf'].mean():.3f}  "
              f"median={result['seq_conf'].median():.3f}")
        print(f"  read_support:   mean={result['read_support'].mean():.3f}  "
              f"median={result['read_support'].median():.3f}")
        print(f"  composite:      mean={result['composite_score'].mean():.3f}  "
              f"median={result['composite_score'].median():.3f}")

    print("\nDone!")


if __name__ == "__main__":
    main()
