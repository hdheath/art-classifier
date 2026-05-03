"""Predict TSS initiation-zone width from genomic position.

Public API consumed by FLAIR's ted.py to set per-locus HDBSCAN clustering
parameters. ted.py currently uses one global ``tss_slack`` for all junction
chains; with width prediction it can scale ``tss_slack`` per-locus based on
the predicted breadth of the initiation zone.

Biological motivation
---------------------
Promoter architecture is bimodal:
  - **Sharp / TATA-driven**: narrow initiation zone (often <50bp).
  - **Broad / CpG-island**: wide initiation zone (often 100-1000bp).

A single global HDBSCAN epsilon can't capture both — this predictor lets
ted.py adapt per-locus.

PolyA cleavage (TTS) is biologically narrow and PolyA_DB is single-bp
resolution; there is no usable TTS width target, so this module exposes
only TSS width.
"""
from __future__ import annotations

import pickle
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from flair_boundary.extract.features import SequenceFeatureExtractor


# Population median of the aggregated FANTOM promoter-region width (bp) —
# used to normalize predictions into a slack multiplier centered at 1.0.
# Refine from training metadata after the first end-to-end run.
DEFAULT_MEDIAN_WIDTH_BP = 200


class WidthPredictor:
    """Loads the trained TSS width regressor and exposes per-position predictions."""

    def __init__(
        self,
        genome_fasta: str,
        tss_width_model: Optional[str] = None,
        median_width: float = DEFAULT_MEDIAN_WIDTH_BP,
    ):
        self.median_width = float(median_width)
        self._model = None
        self._extractor = None

        if tss_width_model and Path(tss_width_model).exists():
            with open(tss_width_model, "rb") as f:
                self._model = pickle.load(f)
            self._extractor = SequenceFeatureExtractor(genome_fasta, "tss")

    def predict(self, chrom: str, pos: int, strand: str) -> float:
        """Predict TSS initiation-zone width in bp for a single position.

        Returns the population median if no model is loaded or feature
        extraction fails (out-of-bounds, missing chromosome).
        """
        if self._model is None:
            return self.median_width
        import pandas as pd
        df = pd.DataFrame([{"chrom": chrom, "pos": int(pos), "strand": strand}])
        X = self._extractor.extract_batch(df)
        if np.any(np.isnan(X)):
            return self.median_width
        log_pred = float(self._model.predict(X)[0])
        return float(np.expm1(log_pred))

    def predict_batch(self, records: List[Tuple[str, int, str]]) -> np.ndarray:
        """Vectorised version of :meth:`predict`. Returns widths in bp."""
        if self._model is None:
            return np.full(len(records), self.median_width, dtype=float)
        import pandas as pd
        df = pd.DataFrame(records, columns=["chrom", "pos", "strand"])
        X = self._extractor.extract_batch(df)
        nan_rows = np.isnan(X).any(axis=1)
        X = np.nan_to_num(X, nan=0.0)
        log_pred = self._model.predict(X)
        widths_bp = np.expm1(log_pred)
        widths_bp[nan_rows] = self.median_width
        return widths_bp

    def tss_slack(self, chrom: str, pos: int, strand: str) -> float:
        """Return a multiplicative HDBSCAN slack factor for TSS clustering.

        Slack > 1 -> broader-than-median promoter, larger HDBSCAN radius.
        Slack < 1 -> sharper-than-median promoter, tighter clustering.
        Clipped to [0.5, 5.0].
        """
        w = self.predict(chrom, pos, strand)
        slack = w / max(self.median_width, 1.0)
        return float(np.clip(slack, 0.5, 5.0))

    # No tts_slack: PolyA_DB is single-bp resolution; ted.py keeps its
    # existing narrow default (1.0) on the 3' axis.
