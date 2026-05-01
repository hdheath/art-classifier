#!/usr/bin/env python3
"""
Sequence feature extraction for TSS and TTS boundary classifiers.

Extracts k-mer composition, positional k-mer profiles, and biological motif
indicators from genomic sequence windows centered on candidate boundaries.

Supports two boundary types:
  - TSS (5' end): Captures promoter-associated sequence signals
  - TTS (3' end): Captures polyadenylation-associated sequence signals

Usage:
    from flair_boundary.extract.features import SequenceFeatureExtractor
    extractor = SequenceFeatureExtractor(genome_fa, boundary_type='tss')
    features = extractor.extract_batch(positions_df)
"""

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pysam

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Window sizes around boundary (upstream, downstream in bp)
WINDOW_CONFIG = {
    "tss": {"upstream": 200, "downstream": 100},   # promoter-proximal region
    "tts": {"upstream": 100, "downstream": 200},   # PAS-proximal region
}

# K-mer sizes to compute
KMER_SIZES = [2, 3, 4]

# PAS (polyadenylation signal) hexamers for TTS
PAS_HEXAMERS = [
    "AATAAA",   # canonical AAUAAA (most common, ~39%)
    "ATTAAA",   # AUUAAA variant (~13.7%)
    "AGTAAA",   #
    "TATAAA",   #
    "CATAAA",   #
    "GATAAA",   #
    "AATATA",   #
    "AATACA",   #
    "AATAGA",   #
    "ACTAAA",   #
    "AACAAA",   #
    "AAGAAA",   #
]

# TATA-box and Initiator motifs for TSS
TSS_MOTIFS = {
    "TATA_box":   r"TATA[AT]A[AT]",   # TATA-box in promoter
    "INR":        r"[TC][TC]A[ACGT][AT][TC][TC]",  # Initiator element
    "DPE":        r"[AG]G[AT][CT][GT]",  # Downstream promoter element
    "BRE":        r"[GC][GC][GA]CGCC",   # TFIIB recognition element
    "CpG_island": r"CG",               # CpG dinucleotide (proxy for island)
}

# Nucleotide alphabet
BASES = ["A", "C", "G", "T"]


# ---------------------------------------------------------------------------
# Core extractor class
# ---------------------------------------------------------------------------

class SequenceFeatureExtractor:
    """
    Extracts sequence-based features from a genome FASTA for boundary training.

    Parameters
    ----------
    genome_fasta : str or Path
        Path to indexed genome FASTA (must have .fai index).
    boundary_type : str
        'tss' or 'tts' — controls window size and motif sets.
    kmer_sizes : list[int]
        Which k-mer sizes to include (default: [2, 3, 4]).
    """

    def __init__(
        self,
        genome_fasta: str,
        boundary_type: str = "tss",
        kmer_sizes: Optional[List[int]] = None,
    ):
        if boundary_type not in ("tss", "tts"):
            raise ValueError(f"boundary_type must be 'tss' or 'tts', got '{boundary_type}'")

        self.boundary_type = boundary_type
        self.kmer_sizes = kmer_sizes or KMER_SIZES
        self.genome = pysam.FastaFile(str(genome_fasta))

        wc = WINDOW_CONFIG[boundary_type]
        self.upstream = wc["upstream"]
        self.downstream = wc["downstream"]
        self.window_size = self.upstream + self.downstream

        # Build k-mer vocabulary
        self._kmer_vocab = self._build_kmer_vocab()

        # Positional k-mer bins (divide window into N bins of 3-mers)
        self._n_bins = 10
        self._bin_size = self.window_size // self._n_bins

    # -----------------------------------------------------------------------
    # Public API
    # -----------------------------------------------------------------------

    def extract_batch(self, df: pd.DataFrame) -> np.ndarray:
        """
        Extract features for a batch of positions.

        Parameters
        ----------
        df : pd.DataFrame
            Must have columns: chrom, pos, strand.
            pos is 0-based coordinate of the boundary.

        Returns
        -------
        np.ndarray
            Shape (n_samples, n_features). NaN rows for failed extractions.
        """
        rows = []
        for _, row in df.iterrows():
            feats = self._extract_one(row["chrom"], int(row["pos"]), row["strand"])
            rows.append(feats)
        return np.array(rows, dtype=np.float32)

    def get_feature_names(self) -> List[str]:
        """Return ordered list of feature names (matches extract_batch columns)."""
        names = []

        # 1. K-mer composition (global window)
        for k in self.kmer_sizes:
            for kmer in sorted(_all_kmers(k)):
                names.append(f"kmer_{kmer}")

        # 2. Positional k-mer profile (3-mers in N bins)
        for b in range(self._n_bins):
            for kmer in sorted(_all_kmers(3)):
                names.append(f"pos_bin{b}_{kmer}")

        # 3. Nucleotide composition in sub-windows
        for region in ["upstream", "mid", "downstream"]:
            for base in BASES:
                names.append(f"base_{region}_{base}")

        # 4. Biological motifs
        if self.boundary_type == "tss":
            for motif_name in TSS_MOTIFS:
                names.append(f"motif_{motif_name}_count")
            for motif_name in TSS_MOTIFS:
                names.append(f"motif_{motif_name}_pos")
            names.append("GC_content")
            names.append("CpG_obs_over_exp")
        else:  # tts
            for hexamer in PAS_HEXAMERS:
                names.append(f"pas_{hexamer}_count")
            for hexamer in PAS_HEXAMERS:
                names.append(f"pas_{hexamer}_min_dist")
            names.append("Arich_frac_downstream50")
            names.append("GC_content")
            names.append("AT_skew")

        return names

    # -----------------------------------------------------------------------
    # Internal helpers
    # -----------------------------------------------------------------------

    def _extract_one(self, chrom: str, pos: int, strand: str) -> List[float]:
        """Extract features for a single boundary position."""
        seq = self._fetch_sequence(chrom, pos, strand)
        if seq is None:
            n_feats = len(self.get_feature_names())
            return [np.nan] * n_feats

        feats = []

        # 1. Global k-mer composition
        for k in self.kmer_sizes:
            feats.extend(self._kmer_composition(seq, k))

        # 2. Positional k-mer profile (3-mers per bin)
        feats.extend(self._positional_kmer_profile(seq))

        # 3. Nucleotide composition per sub-region
        feats.extend(self._regional_composition(seq))

        # 4. Boundary-type specific biological motifs
        if self.boundary_type == "tss":
            feats.extend(self._tss_motif_features(seq))
        else:
            feats.extend(self._tts_motif_features(seq))

        return feats

    def _fetch_sequence(
        self, chrom: str, pos: int, strand: str
    ) -> Optional[str]:
        """
        Fetch genomic sequence window around boundary.

        For TSS/TTS, 'pos' is the boundary coordinate (0-based).
        Window is strand-aware: flipped to 5'→3' orientation for minus strand.
        """
        try:
            chrom_len = self.genome.get_reference_length(chrom)
        except KeyError:
            # Try without/with 'chr' prefix
            alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
            try:
                chrom_len = self.genome.get_reference_length(alt)
                chrom = alt
            except KeyError:
                return None

        if strand == "+":
            start = pos - self.upstream
            end = pos + self.downstream
        else:
            # For minus strand: flip upstream/downstream
            start = pos - self.downstream
            end = pos + self.upstream

        # Clamp to chromosome bounds
        if start < 0 or end > chrom_len:
            return None

        seq = self.genome.fetch(chrom, start, end).upper()

        if strand == "-":
            seq = _reverse_complement(seq)

        if len(seq) != self.window_size:
            return None

        return seq

    def _kmer_composition(self, seq: str, k: int) -> List[float]:
        """Normalized k-mer frequency vector."""
        vocab = sorted(_all_kmers(k))
        counts = {kmer: 0 for kmer in vocab}
        total = 0
        for i in range(len(seq) - k + 1):
            sub = seq[i : i + k]
            if sub in counts:
                counts[sub] += 1
                total += 1
        denom = max(total, 1)
        return [counts[kmer] / denom for kmer in vocab]

    def _positional_kmer_profile(self, seq: str) -> List[float]:
        """3-mer composition in each of N equal-width bins."""
        vocab = sorted(_all_kmers(3))
        feats = []
        for b in range(self._n_bins):
            start = b * self._bin_size
            end = min(start + self._bin_size, len(seq))
            sub = seq[start:end]
            feats.extend(self._kmer_composition(sub, 3))
        return feats

    def _regional_composition(self, seq: str) -> List[float]:
        """Nucleotide composition in upstream, mid, downstream thirds."""
        n = len(seq)
        regions = {
            "upstream":   seq[: n // 3],
            "mid":        seq[n // 3 : 2 * n // 3],
            "downstream": seq[2 * n // 3 :],
        }
        feats = []
        for region_name in ("upstream", "mid", "downstream"):
            sub = regions[region_name]
            denom = max(len(sub), 1)
            for base in BASES:
                feats.append(sub.count(base) / denom)
        return feats

    def _tss_motif_features(self, seq: str) -> List[float]:
        """TSS-specific motif features: TATA-box, INR, CpG density."""
        feats = []

        # Motif counts
        for motif_name, pattern in TSS_MOTIFS.items():
            matches = list(re.finditer(pattern, seq))
            feats.append(float(len(matches)))

        # Position of first motif occurrence (normalized, -1 if absent)
        for motif_name, pattern in TSS_MOTIFS.items():
            m = re.search(pattern, seq)
            pos = (m.start() / len(seq)) if m else -1.0
            feats.append(pos)

        # Global GC content
        gc = (seq.count("G") + seq.count("C")) / max(len(seq), 1)
        feats.append(gc)

        # CpG observed/expected ratio
        cpg_obs = seq.count("CG")
        c_count = seq.count("C")
        g_count = seq.count("G")
        n = len(seq)
        cpg_exp = (c_count / n) * (g_count / n) * n if n > 0 else 0
        cpg_oe = cpg_obs / cpg_exp if cpg_exp > 0 else 0.0
        feats.append(cpg_oe)

        return feats

    def _tts_motif_features(self, seq: str) -> List[float]:
        """TTS-specific motif features: PAS hexamers, A-richness."""
        feats = []

        # PAS hexamer counts in full window
        for hexamer in PAS_HEXAMERS:
            count = seq.count(hexamer)
            feats.append(float(count))

        # Min distance from boundary (position self.upstream) to each hexamer
        boundary_idx = self.upstream
        for hexamer in PAS_HEXAMERS:
            positions = [m.start() for m in re.finditer(hexamer, seq)]
            if positions:
                min_dist = min(abs(p - boundary_idx) for p in positions)
                feats.append(float(min_dist) / self.window_size)
            else:
                feats.append(-1.0)

        # A-richness in 50nt downstream of boundary
        downstream50 = seq[boundary_idx : boundary_idx + 50]
        a_rich = downstream50.count("A") / max(len(downstream50), 1)
        feats.append(a_rich)

        # GC content
        gc = (seq.count("G") + seq.count("C")) / max(len(seq), 1)
        feats.append(gc)

        # AT skew: (A-T)/(A+T)
        a_count = seq.count("A")
        t_count = seq.count("T")
        at_skew = (a_count - t_count) / max(a_count + t_count, 1)
        feats.append(at_skew)

        return feats

    def _build_kmer_vocab(self) -> Dict[int, List[str]]:
        return {k: sorted(_all_kmers(k)) for k in self.kmer_sizes}


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def _all_kmers(k: int) -> List[str]:
    """All DNA k-mers of length k (A/C/G/T only)."""
    if k == 1:
        return list(BASES)
    sub = _all_kmers(k - 1)
    return [b + s for b in BASES for s in sub]


def _reverse_complement(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def load_positive_sites(
    bed_file: str,
    boundary_type: str,
    min_score: float = 0,
    max_sites: Optional[int] = None,
) -> pd.DataFrame:
    """
    Load positive boundary sites from a BED file.

    Parameters
    ----------
    bed_file : str
        BED6 file with cols: chrom, start, end, name, score, strand.
    boundary_type : str
        'tss' or 'tts' — used to determine which coordinate is the boundary.
    min_score : float
        Minimum BED score to include (filters low-confidence sites).
    max_sites : int, optional
        Subsample to this many sites (for memory/speed).

    Returns
    -------
    pd.DataFrame with columns: chrom, pos, strand, score, label=1
    """
    df = pd.read_csv(
        bed_file,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name", "score", "strand"],
        comment="#",
    )

    # Filter low-confidence
    if min_score > 0:
        df = df[df["score"] >= min_score]

    # Choose boundary coordinate
    # TSS: transcription start = start for + strand, end-1 for - strand
    # TTS: transcription termination = end-1 for + strand, start for - strand
    if boundary_type == "tss":
        df["pos"] = df.apply(
            lambda r: r["start"] if r["strand"] == "+" else r["end"] - 1, axis=1
        )
    else:  # tts
        df["pos"] = df.apply(
            lambda r: r["end"] - 1 if r["strand"] == "+" else r["start"], axis=1
        )

    df["label"] = 1

    # Subsample if requested
    if max_sites and len(df) > max_sites:
        df = df.sample(max_sites, random_state=42)

    return df[["chrom", "pos", "strand", "score", "label"]].reset_index(drop=True)


def generate_hard_negatives(
    positive_df: pd.DataFrame,
    gtf_file: str,
    n_negatives: int,
    boundary_type: str,
    min_distance: int = 500,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Generate "hard" negative sites drawn from inside annotated gene bodies.

    The default `generate_negative_sites` samples anywhere on the genome, which
    means negatives are dominated by intergenic positions and the model can
    learn "is this in a gene-rich region?" rather than "is this a TSS/TTS?".

    Hard negatives are sampled from inside gene bodies on the same strand,
    excluding a `min_distance` window around any positive. For TSS, this means
    "internal exon/intron positions on real genes" — the kind of false-positive
    candidate ted.py actually has to reject. Same for TTS.

    Parameters
    ----------
    positive_df : pd.DataFrame
        Reference positive sites (chrom, pos, strand).
    gtf_file : str
        GENCODE GTF for gene coordinates.
    n_negatives : int
        Number of hard negatives to draw.
    boundary_type : str
        'tss' or 'tts' — controls how to avoid the actual TSS/TTS region of
        each gene (we don't want hard negatives that *are* TSSs/TTSs).
    min_distance : int
        Minimum distance from any positive site (bp).
    seed : int
    """
    rng = np.random.default_rng(seed)

    # Parse gene records from GTF
    genes = []
    with open(gtf_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            genes.append((parts[0], int(parts[3]) - 1, int(parts[4]), parts[6]))
    if not genes:
        print("    WARNING: no genes parsed from GTF; returning empty hard negatives")
        return pd.DataFrame(columns=["chrom", "pos", "strand", "score", "label"])

    # Per-(chrom, strand) positive position arrays for distance check
    pos_index: Dict[Tuple[str, str], np.ndarray] = {}
    for (chrom, strand), grp in positive_df.groupby(["chrom", "strand"]):
        pos_index[(chrom, strand)] = np.sort(grp["pos"].values)

    # For each gene, compute the "interior" region we can sample from.
    # Skip the first/last 200bp of the gene to avoid actual TSS/TTS regions.
    SKIP = 200
    gene_intervals = []  # (chrom, lo, hi, strand)
    for chrom, gstart, gend, strand in genes:
        lo = gstart + SKIP
        hi = gend - SKIP
        if hi - lo < 100:
            continue
        gene_intervals.append((chrom, lo, hi, strand))
    if not gene_intervals:
        return pd.DataFrame(columns=["chrom", "pos", "strand", "score", "label"])

    # Length-weight sampling
    lengths = np.array([hi - lo for _, lo, hi, _ in gene_intervals], dtype=np.float64)
    probs = lengths / lengths.sum()

    negatives = []
    attempts = 0
    max_attempts = n_negatives * 30
    while len(negatives) < n_negatives and attempts < max_attempts:
        attempts += 1
        gi = int(rng.choice(len(gene_intervals), p=probs))
        chrom, lo, hi, strand = gene_intervals[gi]
        pos = int(rng.integers(lo, hi))
        ref_positions = pos_index.get((chrom, strand))
        if ref_positions is not None and len(ref_positions) > 0:
            if np.abs(ref_positions - pos).min() < min_distance:
                continue
        negatives.append(
            {"chrom": chrom, "pos": pos, "strand": strand, "score": 0, "label": 0}
        )

    if len(negatives) < n_negatives:
        print(
            f"    WARNING: only {len(negatives)}/{n_negatives} hard negatives "
            f"generated after {attempts} attempts (positives may be too dense)"
        )

    return pd.DataFrame(negatives)


def generate_negative_sites(
    positive_df: pd.DataFrame,
    genome_fasta: str,
    n_negatives: int,
    min_distance: int = 500,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Generate negative (non-boundary) sites for training.

    Strategy: Random genomic positions on valid chromosomes, filtered to
    be at least min_distance away from any positive site.

    Parameters
    ----------
    positive_df : pd.DataFrame
        Positive sites (chrom, pos, strand).
    genome_fasta : str
        Path to indexed genome FASTA.
    n_negatives : int
        Number of negative sites to generate.
    min_distance : int
        Minimum distance from any positive site (bp).
    seed : int
        Random seed.

    Returns
    -------
    pd.DataFrame with columns: chrom, pos, strand, score=0, label=0
    """
    rng = np.random.default_rng(seed)
    fa = pysam.FastaFile(str(genome_fasta))

    # Get chromosome lengths for main chromosomes
    chroms = [c for c in fa.references if re.match(r"^chr[0-9XY]+$", c)]
    chrom_lens = {c: fa.get_reference_length(c) for c in chroms}

    # Pre-compute weighted-by-length sampling probs once (was inside the loop)
    total_len = sum(chrom_lens.values())
    probs = np.array([chrom_lens[c] / total_len for c in chroms])

    # Build positive set for distance checking (chrom -> sorted positions)
    pos_by_chrom: Dict[str, np.ndarray] = {}
    for chrom, grp in positive_df.groupby("chrom"):
        pos_by_chrom[chrom] = np.sort(grp["pos"].values)

    negatives = []
    attempts = 0
    max_attempts = n_negatives * 20

    while len(negatives) < n_negatives and attempts < max_attempts:
        attempts += 1
        chrom = rng.choice(chroms, p=probs)
        pos = int(rng.integers(1000, chrom_lens[chrom] - 1000))

        if chrom in pos_by_chrom:
            dists = np.abs(pos_by_chrom[chrom] - pos)
            if dists.min() < min_distance:
                continue

        strand = rng.choice(["+", "-"])
        negatives.append({"chrom": chrom, "pos": pos, "strand": strand, "score": 0, "label": 0})

    if len(negatives) < n_negatives:
        print(
            f"    WARNING: only {len(negatives)}/{n_negatives} negatives generated "
            f"after {attempts} attempts; consider relaxing min_distance "
            f"(currently {min_distance}bp)"
        )

    fa.close()
    return pd.DataFrame(negatives)
