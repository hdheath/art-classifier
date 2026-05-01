# flAIr
## **F**unctional **L**ong-read **A**nalyzer with **I**ntelligent **R**easoning

Generative transcriptome inference using CNN-based priors for TSS, TES, and splice junctions.

---

## Overview

flAIr replaces heuristic "collapse + filter" transcript inference with a probabilistic framework that:
- Uses **CNN sequence models** to score TSS, TES, and splice sites
- Incorporates **public priors** (CAGE, polyA-DB, GTEx junctions)
- Performs **probabilistic read assignment** via EM algorithm
- Provides **posterior probabilities** and **end uncertainty** for each isoform
- Explicitly models **noise** to reduce artifacts (e.g., internal priming, low-quality reads)

---

## Key Features

- **Probabilistic inference**: Posterior support for each isoform, not just read counts
- **CNN-based scoring**: Sequence-aware evaluation of TSS, TES, and splice sites
- **Artifact suppression**: Internal priming detection, noise component modeling
- **End uncertainty**: Not just point estimates—get confidence intervals for TSS/TES
- **Novel isoform discovery**: Not reference-dependent—finds novel isoforms if read support + sequence plausibility is high

---

## Installation

### Requirements
- Python 3.9+
- PyTorch 2.0+
- pysam, numpy, pandas, pyarrow, scipy
- samtools, bedtools, tabix/bgzip

### Quick Install
```bash
# Clone repository
git clone https://github.com/yourusername/flAIr.git
cd flAIr

# Install dependencies
pip install -r requirements.txt

# Download pre-trained models and priors
./download_data.sh
```

---

## Quick Start

```bash
# 1. Standardize prior files (one-time setup)
python src/extract/standardize_priors.py \
  --tss TSS_db/ \
  --tes TTS_db/polyadb.hg38.weighted.bed \
  --junction splice_db/detail_Linear_splice_annotation.txt \
  --genome GRCh38.fa \
  --output data/ref/

# 2. Extract read evidence from BAM
python src/extract/bam_to_reads.py \
  --bam sample.bam \
  --genes data/ref/genes.bed.gz \
  --output data/derived/reads/sample.parquet

# 3. Run genome-wide inference
python src/infer/infer_genome.py \
  --reads data/derived/reads/sample.parquet \
  --priors data/ref/ \
  --cnn-weights model_weights/cnn_multitask_v1.pt \
  --output data/derived/infer_out/ \
  --threads 16

# 4. Outputs
# - isoforms.gtf: GTF format with posterior in score field
# - isoforms.tsv: Full table with posteriors, abundances, end uncertainty
# - read_assignments.tsv.gz: Per-read assignment probabilities
```

---

## Input Requirements

### Required
- **Aligned BAM**: Sorted and indexed long-read BAM (PacBio or ONT)
- **Reference genome**: GRCh38 FASTA + index

### Priors (provided or build your own)
- **TSS prior**: CAGE peaks (FANTOM5, refTSS)
- **TES prior**: PolyA sites (PolyA_DB, PolyASite)
- **Splice junction prior**: GTEx/Snaptron/RJunBase

### Optional
- **GENCODE GTF**: For gene boundaries and evaluation

---

## Project Structure

```
flAIr/
├── data/
│   ├── ref/                 # Reference files and priors
│   │   ├── genome.fa
│   │   ├── genes.bed.gz
│   │   ├── tss_prior.bed.gz
│   │   ├── tes_prior.bed.gz
│   │   └── junction_prior.tsv.gz
│   ├── bam/                 # Input BAM files
│   └── derived/             # Intermediate and output files
│       ├── reads/           # Extracted read evidence
│       ├── cnn_train/       # CNN training data
│       ├── infer_out/       # Final outputs
│       └── qc/              # QC reports
├── src/
│   ├── extract/             # Data extraction and standardization
│   ├── cnn/                 # CNN model training and inference
│   ├── infer/               # EM inference algorithms
│   └── eval/                # Evaluation and benchmarking
├── model_weights/           # Pre-trained CNN weights
├── tests/                   # Unit and integration tests
├── docs/                    # Documentation
└── examples/                # Tutorial workflows
```

---

## Workflow

### Phase 1: Data Preparation
1. Standardize TSS/TES/junction priors to hg38
2. Extract gene boundaries from GENCODE
3. Validate and index all files

### Phase 2: CNN Training (optional—use pre-trained weights)
1. Generate training windows from priors
2. Train multi-task CNN (TSS, TES, donor, acceptor)
3. Validate model performance

### Phase 3: Read Evidence Extraction
1. Parse BAM CIGAR strings for junctions
2. Extract strand-aware 5' and 3' ends
3. Assign reads to genes
4. Collect QC metrics

### Phase 4: Probabilistic Inference
1. Generate candidate isoforms per gene
2. Score reads against isoforms using CNN + priors
3. Run EM to assign reads probabilistically
4. Compute posterior support for each isoform

### Phase 5: Evaluation
1. Compare to baselines (FLAIR, SQANTI3)
2. Validate TSS/TES accuracy (CAGE, polyA)
3. Assess artifact suppression

---

## Output Files

### isoforms.gtf
Standard GTF with posterior probability in the score field (column 6).

### isoforms.tsv
| Column | Description |
|--------|-------------|
| gene_id | Gene identifier |
| isoform_id | Isoform identifier |
| splice_chain | Ordered junction coordinates |
| posterior | Posterior probability (0-1) |
| abundance | Estimated abundance (TPM-like) |
| tss_mode | Most likely TSS position |
| tss_sd | TSS uncertainty (standard deviation) |
| tes_mode | Most likely TES position |
| tes_sd | TES uncertainty |

### read_assignments.tsv.gz
Per-read assignment probabilities (useful for QC and debugging).

---

## Key Parameters

### CNN Model
- `--tss-window`: Window size for TSS (default: 512bp)
- `--tes-window`: Window size for TES (default: 512bp)
- `--splice-window`: Window size for splice sites (default: 128bp)

### EM Inference
- `--em-iterations`: Maximum EM iterations (default: 50)
- `--convergence-tol`: Convergence tolerance (default: 0.001)
- `--noise-prior`: Prior weight for noise component (default: 0.1)
- `--min-reads`: Minimum reads to infer gene (default: 10)

### Filtering
- `--min-posterior`: Minimum posterior to report isoform (default: 0.1)
- `--min-cnn-score`: Minimum CNN score for novel junctions (default: 0.5)

---

## Performance

### Benchmarks (50M reads, 16 cores)
- **Runtime**: ~3.5 hours
- **Memory**: ~28GB RAM
- **Genes processed**: ~19,000
- **Isoforms inferred**: ~45,000

### Accuracy vs FLAIR
- **Internal priming rate**: 4.2% (flAIr) vs 18.7% (FLAIR)
- **Singleton suppression**: 16% (flAIr) vs 52% (FLAIR)
- **TSS accuracy**: 87bp median to CAGE vs 145bp (FLAIR)

---

## Citation

If you use flAIr in your research, please cite:

```
[Publication pending]
```

---

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

---

## License

MIT License - see [LICENSE](LICENSE) for details.

---

## Contact

- **Issues**: [GitHub Issues](https://github.com/yourusername/flAIr/issues)
- **Email**: hdheath@stanford.edu

---

## Acknowledgments

- FANTOM5 consortium for CAGE data
- PolyA_DB and PolyASite for TES data
- GTEx and RJunBase for splice junction data
- Brooks Lab for computational resources

---

## References

1. FANTOM5: Forrest et al. Nature (2014)
2. PolyA_DB: Wang et al. NAR (2017)
3. refTSS: Zhao et al. NAR (2021)
4. FLAIR: Tang et al. Nat Commun (2020)
5. SQANTI3: Tardaguila et al. Genome Res (2018)
