# flAIr Boundary Classifier

Hybrid TSS/TTS boundary classifier for FLAIR transcriptome that combines:
1. **Sequence features** — k-mer composition, positional profiles, and biological motifs from GRCh38
2. **Read-level evidence** — soft-clipping, A-richness, read count, and boundary dispersion from BAM
3. **Composite scoring** — weighted product with artifact penalty

## Environment

```bash
conda activate flair-ml   # has xgboost 3.1, pysam 0.23, sklearn 1.8
```

## Quick Start (full training run)

```bash
conda activate flair-ml
cd /private/groups/brookslab/hdheath/projects/flAIr
bash src/run_training.sh                   # trains both TSS and TTS
bash src/run_training.sh --type tss        # TSS only
bash src/run_training.sh --type tts        # TTS only
bash src/run_training.sh --models "xgboost logistic rf xgboost_calibrated"
```

Outputs:
- `data/derived/tss_features.npz` / `tts_features.npz` — feature matrices
- `model_weights/tss/xgboost.pkl` etc. — trained models
- `data/derived/eval/tss/` — PR curves, feature importance plots, per-chrom breakdown

## Module Structure

```
src/flair_boundary/
├── extract/
│   ├── features.py         # SequenceFeatureExtractor, negative site generation
│   └── prepare_data.py     # CLI: load databases → extract features → save .npz
├── train/
│   └── train.py            # CLI: XGBoost / LogReg / RF training with CV
├── infer/
│   └── infer.py            # CLI+API: HybridBoundaryScorer, composite scoring
├── eval/
│   └── evaluate.py         # CLI: PR/ROC curves, feature importance, per-chrom eval
└── tests/
    └── test_features.py    # Unit tests for feature extraction
```

## Databases

| Type | File | Entries | Build |
|------|------|---------|-------|
| TSS | `TSS_db/FANTOM_TSS_human.bed` | 1,048,124 | hg38 |
| TSS | `TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt` | 241,049 | hg38 |
| TTS | `TTS_db/polyadb.hg38.weighted.bed` | 311,500 | hg38 |
| TTS | `TTS_db/polyAdb_human.PAS.txt` | 311,595 | hg38 (PAS annotations) |

## Features (TSS: 1000, TTS: 1015)

- **K-mer composition** (2,3,4-mers): 16+64+256 = 336 features
- **Positional k-mer profile** (3-mers in 10 bins): 640 features
- **Nucleotide composition** (per sub-window): 12 features
- **TSS motifs**: TATA-box, INR, DPE, BRE, CpG density + positions
- **TTS motifs**: 12 PAS hexamers (AATAAA, ATTAAA, ...) + A-richness

## Inference API

```python
from flair_boundary.infer.infer import HybridBoundaryScorer
import pandas as pd

scorer = HybridBoundaryScorer(
    genome_fasta="data/ref/GRCh38.primary_assembly.genome.fa",
    tss_model_path="model_weights/tss/xgboost.pkl",
    tts_model_path="model_weights/tts/xgboost.pkl",
    seq_weight=0.5,
    read_weight=0.5,
)

# candidates: DataFrame with chrom, pos, strand, [splice_chain]
scored = scorer.score_candidates(candidates, boundary_type="tss", bam_path="reads.bam")
selected = scorer.select_boundaries(scored, min_score=0.3, multi_boundary_threshold=0.5)
```

## Composite Score

```
score(b) = (0.5 * seq_conf + 0.5 * read_support) * (1 - artifact_penalty)
```

- `seq_conf`: XGBoost probability (sequence classifier)
- `read_support`: log1p(n_reads) normalized within cluster
- `artifact_penalty`: soft-clip fraction (TSS) or A-richness (TTS)

Boundaries with score ≥ 0.3 are reported. If two candidates exceed 0.5 and are >50nt apart, both are reported as distinct isoforms.

## Running Tests

```bash
cd /private/groups/brookslab/hdheath/projects/flAIr
PYTHONPATH=src python src/flair_boundary/tests/test_features.py
```
