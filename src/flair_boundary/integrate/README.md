# Integration with FLAIR ted.py

This module exposes the trained boundary models to FLAIR's TED (Transcript
End Determination) module so that:

1. Per-locus width predictions tune HDBSCAN clustering parameters
   (`--ted_tss_slack`, `--ted_end_window`).
2. Boundary classifier probabilities feed into TED's `model_score`
   reality term (this already works in ted.py — `TedScorer._tss_model`).

## Files in `model_weights/{tss,tts}/`

| File | Purpose |
|---|---|
| `xgboost.pkl` | Classifier — P(true boundary \| sequence). Consumed by ted.py via `TedScorer._tss_model`. |
| `xgboost_width.pkl` | **NEW.** Regressor — predicts initiation/termination zone width in bp. Consumed by `WidthPredictor` (this module). |

## Patch needed in `flair4/src/flair/ted.py`

Currently, `_hdbscan_cluster_2d` takes a scalar `tss_slack`. To use
per-locus slack, change the signature so `tss_slack` accepts a callable:

```python
def _hdbscan_cluster_2d(starts, ends, min_cluster_size, min_samples,
                        tss_slack=1.0, end_window=100):
    # If tss_slack is callable, evaluate per-position and use the median
    # for the cluster (we still need a scalar slack inside HDBSCAN).
    if callable(tss_slack):
        slacks = np.array([tss_slack(s) for s in starts])
        slack_value = float(np.median(slacks))
    else:
        slack_value = float(tss_slack)
    scaled_starts = starts.astype(np.float64) / slack_value
    ...
```

Then in `flair_transcriptome.py` (or wherever ted is called):

```python
from flair_boundary.integrate import WidthPredictor

wp = WidthPredictor(genome_fasta=args.genome,
                    tss_width_model=args.tss_width_model,
                    tts_width_model=args.tts_width_model)

ted_collapse_end_groups(
    ...,
    tss_slack=lambda pos: wp.tss_slack(chrom, pos, strand),
    ...,
)
```

## Why this isn't done as a single repo

ted.py lives in `brookslab/tools/flair4` and is shared infrastructure.
The width model lives here in `art-classifier` because it's a research
artifact under active iteration. The integration is one small patch to
ted.py + an import — keeping them separate avoids tangling release
cycles.
