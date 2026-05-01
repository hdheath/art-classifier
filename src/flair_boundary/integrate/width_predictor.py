"""Predict TSS/TTS initiation/termination zone width from genomic position.

This is the public API consumed by FLAIR's ted.py to set per-locus HDBSCAN
clustering parameters. ted.py currently uses one global ``tss_slack`` for all
junction chains; with width prediction it can scale ``tss_slack`` per-locus
based on the predicted breadth of the initiation zone.

Biological motivation
---------------------
Promoter architecture is bimodal:
  - **Sharp / TATA-driven**: narrow initiation zone (often <20bp), e.g.
    housekeeping genes, tissue-specific genes with TATA boxes.
  - **Broad / CpG-island**: wide initiation zone (often 100–500bp), e.g.
    most CpG-island promoters, dispersed initiation.

Termination is generally narrower (polyA cleavage is more discrete), but
also shows variation. A single global HDBSCAN epsilon can't capture both —
this predictor lets ted.py adapt per-locus.

Usage
-----
::

    from flair_boundary.integrate import WidthPredictor

    wp = WidthPredictor(
        genome_fasta="GRCh38.fa",
        tss_width_model="model_weights/tss/xgboost_width.pkl",
        tts_width_model="model_weights/tts/xgboost_width.pkl",
    )

    # Single position
    width_bp = wp.predict("chr7", 100_000, "+", "tss")

    # Convert width to HDBSCAN slack: ratio of predicted to median width
    slack = wp.tss_slack("chr7", 100_000, "+")
"""
from __future__ import annotations

import pickle
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from flair_boundary.extract.features import SequenceFeatureExtractor


# Population medians of peak width (bp) — used to normalize predictions into
# a slack multiplier centered at 1.0. Estimated from FANTOM CAGE peaks
# (TSS) and PolyA_DB (TTS) on hg38; refine after first training run.
DEFAULT_MEDIAN_WIDTH_BP = {"tss": 50, "tts": 1}


class WidthPredictor:
    """Loads trained width regressors and exposes per-position predictions.

    Parameters
    ----------
    genome_fasta : str
        Path to indexed GRCh38 FASTA.
    tss_width_model : str or None
        Path to trained TSS width regressor (.pkl). If None, TSS width
        falls back to ``DEFAULT_MEDIAN_WIDTH_BP['tss']``.
    tts_width_model : str or None
        Path to trained TTS width regressor (.pkl).
    median_widths : dict or None
        Per-type median widths used to normalize predictions into a slack
        multiplier. Defaults to :data:`DEFAULT_MEDIAN_WIDTH_BP`.
    """

    def __init__(
        self,
        genome_fasta: str,
        tss_width_model: Optional[str] = None,
        tts_width_model: Optional[str] = None,
        median_widths: Optional[dict] = None,
    ):
        self.median_widths = {**DEFAULT_MEDIAN_WIDTH_BP, **(median_widths or {})}
        self._models: dict = {}
        self._extractors: dict = {}

        for btype, path in [("tss", tss_width_model), ("tts", tts_width_model)]:
            if path and Path(path).exists():
                with open(path, "rb") as f:
                    self._models[btype] = pickle.load(f)
                self._extractors[btype] = SequenceFeatureExtractor(genome_fasta, btype)

    # -----------------------------------------------------------------------
    # Single-position queries
    # -----------------------------------------------------------------------

    def predict(self, chrom: str, pos: int, strand: str, btype: str) -> float:
        """Predict peak width in bp for a single position.

        Returns the population median if no model is loaded for this type
        or feature extraction fails (out-of-bounds, missing chromosome).
        """
        if btype not in self._models:
            return float(self.median_widths[btype])
        import pandas as pd
        df = pd.DataFrame([{"chrom": chrom, "pos": int(pos), "strand": strand}])
        X = self._extractors[btype].extract_batch(df)
        if np.any(np.isnan(X)):
            return float(self.median_widths[btype])
        log_pred = float(self._models[btype].predict(X)[0])
        return float(np.expm1(log_pred))

    def predict_batch(
        self,
        records: List[Tuple[str, int, str]],
        btype: str,
    ) -> np.ndarray:
        """Vectorised version of :meth:`predict`. Returns widths in bp."""
        if btype not in self._models:
            return np.full(len(records), self.median_widths[btype], dtype=float)
        import pandas as pd
        df = pd.DataFrame(records, columns=["chrom", "pos", "strand"])
        X = self._extractors[btype].extract_batch(df)
        # Fill NaN feature rows with population median to avoid prediction errors
        nan_rows = np.isnan(X).any(axis=1)
        X = np.nan_to_num(X, nan=0.0)
        log_pred = self._models[btype].predict(X)
        widths_bp = np.expm1(log_pred)
        widths_bp[nan_rows] = self.median_widths[btype]
        return widths_bp

    # -----------------------------------------------------------------------
    # HDBSCAN slack helpers (the API ted.py actually calls)
    # -----------------------------------------------------------------------

    def tss_slack(self, chrom: str, pos: int, strand: str) -> float:
        """Return a multiplicative slack factor for HDBSCAN clustering at TSS.

        Slack > 1 means "this promoter is broader than the median, give
        HDBSCAN a larger reachability radius along the TSS axis." Slack < 1
        means "this is a sharp promoter, cluster tighter." Used by ted.py
        as the per-locus replacement for the global ``--ted_tss_slack``.

        Slack is clipped to [0.5, 5.0] to keep clustering numerically sane —
        a width prediction 10× the median is more likely a noisy outlier
        than a truly enormous initiation zone.
        """
        w = self.predict(chrom, pos, strand, "tss")
        slack = w / max(self.median_widths["tss"], 1.0)
        return float(np.clip(slack, 0.5, 5.0))

    def tts_slack(self, chrom: str, pos: int, strand: str) -> float:
        """Same as :meth:`tss_slack` but for the 3' end. Tighter clip
        because polyA sites are mostly narrow."""
        w = self.predict(chrom, pos, strand, "tts")
        slack = w / max(self.median_widths["tts"], 1.0)
        return float(np.clip(slack, 0.5, 3.0))
