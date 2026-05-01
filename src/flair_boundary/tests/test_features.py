#!/usr/bin/env python3
"""
Unit tests for sequence feature extraction.

Run with:
    conda activate flair-ml
    python -m pytest src/flair_boundary/tests/test_features.py -v
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

try:
    import pytest
    HAS_PYTEST = True
except ImportError:
    HAS_PYTEST = False
    # Minimal stub so tests can run standalone
    class _Mark:
        @staticmethod
        def skipif(cond, reason=""):
            def decorator(fn):
                if cond:
                    def wrapper(*a, **kw):
                        print(f"  SKIP: {reason}")
                    wrapper.__name__ = fn.__name__
                    return wrapper
                return fn
            return decorator
    class pytest:  # type: ignore
        mark = _Mark()

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from flair_boundary.extract.features import (
    SequenceFeatureExtractor,
    _all_kmers,
    _reverse_complement,
)

# Path to genome (adjust if needed)
GENOME = Path(__file__).resolve().parents[3] / "data/ref/GRCh38.primary_assembly.genome.fa"
SKIP_GENOME = not GENOME.exists()


# ---------------------------------------------------------------------------
# K-mer utilities
# ---------------------------------------------------------------------------

def test_kmer_counts():
    assert len(_all_kmers(1)) == 4
    assert len(_all_kmers(2)) == 16
    assert len(_all_kmers(3)) == 64
    assert len(_all_kmers(4)) == 256


def test_reverse_complement():
    assert _reverse_complement("AACGT") == "ACGTT"
    assert _reverse_complement("ATCG") == "CGAT"
    assert _reverse_complement("") == ""


# ---------------------------------------------------------------------------
# Feature extraction (requires genome)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(SKIP_GENOME, reason="Genome FASTA not found")
def test_tss_feature_shape():
    ext = SequenceFeatureExtractor(str(GENOME), "tss")
    names = ext.get_feature_names()
    df = pd.DataFrame([{"chrom": "chr1", "pos": 1_000_000, "strand": "+"}])
    X = ext.extract_batch(df)
    assert X.shape == (1, len(names)), f"Shape mismatch: {X.shape} vs (1, {len(names)})"
    assert not np.isnan(X).any(), "NaN in feature matrix"


@pytest.mark.skipif(SKIP_GENOME, reason="Genome FASTA not found")
def test_tts_feature_shape():
    ext = SequenceFeatureExtractor(str(GENOME), "tts")
    names = ext.get_feature_names()
    df = pd.DataFrame([{"chrom": "chr1", "pos": 5_000_000, "strand": "-"}])
    X = ext.extract_batch(df)
    assert X.shape == (1, len(names))
    assert not np.isnan(X).any()


@pytest.mark.skipif(SKIP_GENOME, reason="Genome FASTA not found")
def test_out_of_bounds_returns_nan():
    """Positions too close to chromosome edge should return NaN row."""
    ext = SequenceFeatureExtractor(str(GENOME), "tss")
    df = pd.DataFrame([{"chrom": "chr1", "pos": 10, "strand": "+"}])
    X = ext.extract_batch(df)
    assert np.isnan(X).all(), "Expected all NaN for out-of-bounds position"


@pytest.mark.skipif(SKIP_GENOME, reason="Genome FASTA not found")
def test_strand_symmetry():
    """Features for + and - strand at same position should differ (strand-aware)."""
    ext = SequenceFeatureExtractor(str(GENOME), "tss")
    pos = 2_000_000
    df_plus = pd.DataFrame([{"chrom": "chr1", "pos": pos, "strand": "+"}])
    df_minus = pd.DataFrame([{"chrom": "chr1", "pos": pos, "strand": "-"}])
    X_plus = ext.extract_batch(df_plus)
    X_minus = ext.extract_batch(df_minus)
    # Features should differ due to reverse complement
    assert not np.allclose(X_plus, X_minus), "Plus and minus strand features should differ"


@pytest.mark.skipif(SKIP_GENOME, reason="Genome FASTA not found")
def test_batch_consistent():
    """Batch extraction should give same result as individual calls."""
    ext = SequenceFeatureExtractor(str(GENOME), "tss")
    positions = [1_000_000, 2_000_000, 3_000_000]
    df = pd.DataFrame([{"chrom": "chr1", "pos": p, "strand": "+"} for p in positions])
    X_batch = ext.extract_batch(df)

    for i, p in enumerate(positions):
        df_single = pd.DataFrame([{"chrom": "chr1", "pos": p, "strand": "+"}])
        X_single = ext.extract_batch(df_single)
        np.testing.assert_array_almost_equal(X_batch[i], X_single[0])


@pytest.mark.skipif(SKIP_GENOME, reason="Genome FASTA not found")
def test_kmer_composition_sums_to_one():
    """K-mer frequency vectors should sum to ~1."""
    ext = SequenceFeatureExtractor(str(GENOME), "tss")
    df = pd.DataFrame([{"chrom": "chr1", "pos": 1_500_000, "strand": "+"}])
    X = ext.extract_batch(df)
    names = ext.get_feature_names()

    # Find indices of 2-mer features
    kmer2_idx = [i for i, n in enumerate(names) if n.startswith("kmer_") and len(n) == 7]
    kmer2_sum = X[0, kmer2_idx].sum()
    assert abs(kmer2_sum - 1.0) < 0.01, f"2-mer frequencies don't sum to 1: {kmer2_sum}"


if __name__ == "__main__":
    if HAS_PYTEST:
        pytest.main([__file__, "-v"])
    else:
        # Standalone runner
        import traceback
        tests = [
            test_kmer_counts,
            test_reverse_complement,
            test_tss_feature_shape,
            test_tts_feature_shape,
            test_out_of_bounds_returns_nan,
            test_strand_symmetry,
            test_batch_consistent,
            test_kmer_composition_sums_to_one,
        ]
        passed = failed = skipped = 0
        for fn in tests:
            try:
                fn()
                print(f"  PASS: {fn.__name__}")
                passed += 1
            except Exception as e:
                if "SKIP" in str(e):
                    skipped += 1
                else:
                    print(f"  FAIL: {fn.__name__}: {e}")
                    traceback.print_exc()
                    failed += 1
        print(f"\n{passed} passed, {failed} failed, {skipped} skipped")
