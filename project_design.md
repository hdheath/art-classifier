```markdown
# Generative Transcriptome Inference with CNN Priors (TSS / TES / Splice)
A practical blueprint for building a probabilistic long-read transcript inference system that uses:
- an aligned long-read BAM (PacBio or ONT),
- public priors for **TSS**, **TES/polyA**, and **splice junctions**,
- and **CNN sequence models** to score candidate 5′ ends, 3′ ends, and splice sites.

The goal is to replace heuristic “collapse + filter” with:
- **probabilistic read assignment** to isoforms,
- **explicit noise modeling**,
- and **posterior probabilities** on inferred isoforms.

This README describes an end-to-end design and an implementable step plan.

---

## 0) Inputs

### Required
- **Aligned BAM** (sorted + indexed)
  - `sample.bam`
  - `sample.bam.bai`
- **Reference genome (hg38)** FASTA + index
  - `GRCh38.primary_assembly.genome.fa`
  - `GRCh38.primary_assembly.genome.fa.fai`

### Priors (you already assembled these)
- **TSS prior** (refTSS/CAGE-derived peaks), hg38 coordinates
  - `tss_prior.bed.gz` (+ `.tbi`)
- **TES / PAS prior** (PolyA_DB / PolyASite), hg38 coordinates
  - `tes_prior.bed.gz` (+ `.tbi`)
- **Splice junction prior** (GTEx/Snaptron/recount-style), hg38 coordinates
  - `junction_prior.tsv.gz` (+ `.tbi`)  
  Format recommendation:
```

chr  donor  acceptor  strand  weight

```

### Optional (highly recommended)
- GENCODE/MANE GTF (for gene boundaries / evaluation)
- `gencode.gtf.gz`

---

## 1) Repository layout (suggested)
```

.
├── data/
│   ├── ref/
│   │   ├── genome.fa
│   │   ├── genome.fa.fai
│   │   ├── genes.bed.gz
│   │   ├── tss_prior.bed.gz
│   │   ├── tes_prior.bed.gz
│   │   └── junction_prior.tsv.gz
│   ├── bam/
│   │   ├── sample.bam
│   │   └── sample.bam.bai
│   └── derived/
│       ├── reads.parquet
│       ├── candidates/
│       ├── cnn_train/
│       └── cnn_infer/
├── src/
│   ├── extract/
│   ├── candidates/
│   ├── cnn/
│   ├── infer/
│   └── eval/
└── README.md

```

---

## 2) Overview of the approach

### A) Build CNN “site plausibility” models from public priors + sequence
We train CNNs to estimate:
- `P_TSS(seq_window)` = probability a position is a real TSS
- `P_TES(seq_window)` = probability a position is a real TES / cleavage+polyA site
- `P_donor(seq_window)` and `P_acceptor(seq_window)` = splice site plausibility
- optionally: `P_internal_priming(seq_window)` = probability a 3′ end is internal priming artifact

These CNNs are trained **once** from public priors (and can be shipped as weights).

### B) Extract long-read evidence from BAM
From each aligned read:
- gene/locus assignment
- splice junction chain (from CIGAR `N`)
- strand-aware 5′ end and 3′ end
- QC covariates (MAPQ, softclip, NM, aligned length, etc.)

### C) Probabilistic isoform inference (per-gene mixture model)
Candidate isoforms come from:
- observed splice chains (unique junction sequences per gene)
- plus optional annotated splice chains

Reads are assigned probabilistically to isoforms (and to a **noise component**).
The likelihood uses:
- splice-chain match (with splice CNN scores + junction prior weights)
- 5′ end likelihood (TSS prior + TSS CNN score)
- 3′ end likelihood (TES prior + TES CNN score + internal priming CNN penalty)
- read QC likelihood (platform/run dependent)

Outputs:
- isoforms with posterior support and abundance
- inferred TSS/TES distributions (uncertainty)
- per-read assignment probabilities
- a GTF/BED export (a *projection* of the posterior)

---

## 3) Step-by-step implementation plan

### Step 1 — Standardize and index priors
All priors should be:
- **hg38**
- **bgzip + tabix** indexed
- consistent `chr` naming

Recommended:
- `tss_prior.bed.gz` is BED6 with weight in score:
```

chr start end name score strand

```
- `tes_prior.bed.gz` is BED6 similarly
- `junction_prior.tsv.gz` is tab-separated, tabix-indexed (bgzip) with `chr` as first column.

---

### Step 2 — Create gene/locus boundaries (`genes.bed.gz`)
Use a GTF to create per-gene intervals (strand-aware). Minimal required fields:
```

chr  start  end  gene_id  score  strand

```

This is used only to make inference tractable (per gene/locus). It does NOT force isoforms to match annotation.

---

### Step 3 — Extract long-read evidence table from BAM
Create `reads.parquet` with one row per read:

Required columns:
- `read_id`
- `chr`, `strand`
- `gene_id` (overlap with `genes.bed`, strand-aware)
- `tss_pos` (5′ end coordinate, strand-aware)
- `tes_pos` (3′ end coordinate, strand-aware)
- `junctions` (ordered list: `(donor, acceptor)` pairs)
- `mapq`, `aligned_len`, `softclip_5p`, `softclip_3p`, `nm` (if present)
- `platform`, `sample_id`

Notes:
- For long reads, `tss_pos` and `tes_pos` are taken from the aligned read ends:
  - if `strand == '+'`: 5′ = alignment start, 3′ = alignment end
  - if `strand == '-'`: 5′ = alignment end, 3′ = alignment start
- Extract junctions from `N` operations in the CIGAR.

---

## 4) CNN training data construction (from your priors)

### 4.1 Define sequence windows
Typical window sizes:
- TSS: ±256 bp (512 bp window)
- TES: ±256 bp
- Donor/acceptor: ±64 bp centered on splice site (128 bp window)
  - donor centered at exon-intron boundary
  - acceptor centered at intron-exon boundary

Represent each window as one-hot encoding (A,C,G,T,N).

### 4.2 Positive labels (from priors)
- TSS positives: centers at TSS prior peak summits
- TES positives: centers at TES prior / PAS sites
- Donor/acceptor positives:
  - from junction priors (GTEx) and/or GENCODE junctions
  - donor site = junction donor coordinate
  - acceptor site = junction acceptor coordinate

### 4.3 Negative labels (critical)
Negatives should be realistic and matched for GC content / mappability when possible.

Recommended negative sampling:
- TSS negatives:
  - same chromosome and strand, >1 kb away from any known TSS peak
  - optionally matched GC distribution
- TES negatives:
  - >1 kb away from known PAS sites
  - include A-rich regions to teach internal priming separation
- Splice negatives:
  - random positions within introns/exons
  - “near-miss” negatives (shifted ±3–20 bp) to help the model learn sharp boundaries

### 4.4 Multi-task vs separate CNNs
Recommended: **multi-task CNN** with four heads:
- `TSS` (binary)
- `TES` (binary)
- `DONOR` (binary)
- `ACCEPTOR` (binary)

Optional fifth head:
- `INTERNAL_PRIMING` (binary)  
  Positives can be constructed from:
  - 3′ ends in A-rich motifs without PAS signals (heuristic + priors)
  - or known internal priming annotations if available

---

## 5) CNN architecture (baseline)
A simple 1D CNN works well for motif-scale biology.

Example design:
- Input: one-hot `L x 4` (L=512 for TSS/TES; L=128 for splice sites)
- 3–6 convolution blocks (Conv1D → BatchNorm → ReLU → Dropout)
- Global max pooling
- Dense layer
- Output heads (sigmoid) for tasks

Loss:
- sum of binary cross-entropies with class weights (handle imbalance)

Metrics:
- AUROC, AUPRC per head
- calibration curve (important for probabilistic inference)

---

## 6) How CNN scores enter the probabilistic isoform model

### 6.1 Candidate isoforms
Per gene:
- Build candidate isoforms from **unique observed splice chains** in reads
- Optionally add annotated splice chains from GTF
- Do NOT explode the space; prune extremely rare chains early (e.g., only keep chains seen ≥2 reads, but allow 1-read chains to remain as “candidates” if they have high splice plausibility)

### 6.2 Likelihood components per read r and isoform k
We compute:
- `L_splice(r,k)`:
  - For each junction in read, score donor+acceptor with CNN
  - Combine with junction prior weight if junction is known
  - Penalize junctions that are both:
    - not in prior, and
    - have weak donor/acceptor CNN scores
- `L_tss(r,k)`:
  - score read’s 5′ end with:
    - TSS prior density near position
    - CNN TSS score at that position
- `L_tes(r,k)`:
  - score read’s 3′ end with:
    - TES prior density near position
    - CNN TES score
    - internal priming penalty (if used)
- `L_qc(r)`:
  - MAPQ, softclip, NM, platform/run parameters

Then:
```

log P(r | k) =
log L_splice(r,k) +
log L_tss(r,k) +
log L_tes(r,k) +
log L_qc(r)

```

### 6.3 Noise component (mandatory)
Add k=0:
- broad end distributions
- weak splice constraints
- QC favors low-quality reads

This prevents inventing isoforms to explain junk.

---

## 7) Inference algorithm (first working version)
Use per-gene EM:

1) Initialize isoforms and mixture weights
2) E-step: responsibilities
```

gamma[r,k] ∝ pi[k] * P(r | k)

```
include k=0 noise
3) M-step:
- update pi[k]
- optionally update isoform-specific end distributions
4) prune low-support isoforms
5) repeat until convergence

Outputs per gene:
- isoform splice chains
- posterior support / abundance
- end distributions (mode ± uncertainty)
- per-read assignment probabilities (including noise probability)

---

## 8) Evaluation plan (must-have)
### Internal consistency (no external truth required)
- Replicate reproducibility:
- high-posterior isoforms should recur across libraries
- Calibration:
- isoforms with posterior 0.9 should reproduce much more than 0.2
- Artifact suppression:
- fraction of TES assigned to A-rich internal priming regions should drop

### External agreement (optional)
- Compare inferred TSS to refTSS/CAGE peaks: distance distribution
- Compare inferred TES to PAS atlas: distance distribution
- Compare junctions to GTEx known junction prevalence

Baseline comparisons:
- FLAIR collapse + SQANTI filtering
- IsoSeq3/TALON (as appropriate)

---

## 9) Practical CLI “milestones” (what to build first)

### Milestone 1: CNN priors
- Build training windows from priors + hg38 sequence
- Train CNNs and export weights
- Provide a script to score any BED position list

### Milestone 2: Evidence table extraction
- BAM → reads.parquet (ends + junctions + QC + gene_id)

### Milestone 3: Probabilistic read assignment
- For a single gene:
- candidate isoforms from observed splice chains
- compute P(r|k) with CNN + priors
- EM to assign reads and infer isoforms

### Milestone 4: Scale genome-wide
- Run per-gene inference across all genes/loci
- Export GTF + posterior tables

---

## 10) Notes on “reference reliance”
This system is NOT “in reference → keep, else drop”.
It is “prior + sequence + read consistency → posterior probability”.

Novel biology is allowed if:
- multiple reads support it consistently
- splice sites are plausible by CNN
- ends are plausible by CNN even if not in atlas
- and the overall likelihood beats the noise model

---

## 11) Minimum dependencies (recommended)
- Python: `pysam`, `numpy`, `pandas`, `pyarrow`, `scipy`
- Genomics: `bedtools`, `samtools`, `tabix/bgzip`
- Deep learning: `pytorch` (or `tensorflow`)
- Optional: `pyranges` for interval joins

---

## 12) Data artifacts you should produce (for users)
- `model_weights/` for CNN heads
- `tss_prior`, `tes_prior`, `junction_prior` packaged as indexed files
- `infer_out/`
- `isoforms.gtf`
- `isoforms.tsv` (posterior, abundance, end uncertainty)
- `read_assignments.tsv.gz` (optional)
- `qc_summary.json`

---

## 13) Suggested default hyperparameters
- TSS/TES window: 512 bp
- splice window: 128 bp
- EM iterations: 20–50
- noise mixture prior: 0.05–0.30 (depends on platform/library quality)
- junction jitter tolerance (scoring): ±3–10 bp (platform-specific)
- end clustering smoothing for priors: ±10–50 bp

---

## 14) What you should store as “priors” after CNN training
For fast inference, precompute:
- genome-wide CNN score tracks at candidate sites (optional)
- or score on-the-fly only at:
- read ends
- junction donor/acceptor coordinates

On-the-fly scoring is usually sufficient.

---

## 15) Implementation caution: avoid exploding candidate space
Candidate isoforms should be restricted initially to:
- observed splice chains per gene
- plus a small set of annotated chains
Then use the posterior to decide which survive.

Do NOT enumerate all paths of a splice graph at the start.

---

## 16) Deliverable definition (what "done" looks like)
For any input BAM:
- produce isoforms with posterior support and end uncertainty
- reduce singleton/garbage isoforms vs heuristic collapse
- improve end plausibility (TSS/TES distance to priors)
- provide per-read noise probabilities

---

## 17) Implementation Phases with Deliverables and Testing

### Phase 1: Data Preparation and Validation

#### Deliverables:
1. **Prior file standardization pipeline** (`src/extract/standardize_priors.py`)
   - Input: Raw TSS/TES/junction files
   - Output: Standardized, indexed BED/TSV files with consistent chromosome naming
   - File: `data/ref/{tss,tes,junction}_prior.bed.gz` + `.tbi`

2. **Gene boundary extraction** (`src/extract/create_gene_loci.py`)
   - Input: GENCODE GTF
   - Output: `data/ref/genes.bed.gz` (strand-aware gene boundaries)

3. **Data quality report** (`qc/prior_summary.html`)
   - TSS/TES peak counts per chromosome
   - Junction count statistics
   - Genome coverage statistics
   - Strand distribution

#### Testing Methods:
- **Unit tests:**
  - Verify tabix indexing succeeds for all priors
  - Check chromosome naming consistency (all `chr*` or all without `chr`)
  - Validate BED coordinate sorting
  - Ensure no overlapping gene boundaries on same strand
- **Integration tests:**
  - Query random genomic regions with tabix
  - Verify all chromosomes present in both genome and priors
  - Check for edge cases (chrM, chrX, chrY)
- **Validation metrics:**
  - TSS prior coverage: >80% of GENCODE annotated genes within 500bp
  - Junction prior completeness: >90% of GENCODE canonical junctions present
  - No duplicate entries in any prior file

---

### Phase 2: BAM Evidence Extraction

#### Deliverables:
1. **Read extraction pipeline** (`src/extract/bam_to_reads.py`)
   - Input: Aligned BAM, gene boundaries
   - Output: `data/derived/reads.parquet`
   - Required columns: `read_id`, `chr`, `strand`, `gene_id`, `tss_pos`, `tes_pos`, `junctions`, `mapq`, `aligned_len`, `softclip_5p`, `softclip_3p`, `nm`, `platform`, `sample_id`

2. **Read QC report** (`qc/read_extraction_summary.html`)
   - Reads per gene distribution
   - Strand assignment statistics
   - Junction complexity distribution (0, 1, 2-5, 6+ junctions per read)
   - MAPQ distribution
   - Softclipping statistics

#### Testing Methods:
- **Unit tests:**
  - Parse CIGAR strings correctly (test cases with various `M`, `N`, `D`, `I`, `S` operations)
  - Strand-aware coordinate extraction (compare + vs - strand reads)
  - Gene assignment overlap logic
- **Integration tests:**
  - Compare extracted junction coordinates to manual IGV inspection (sample 20 reads)
  - Verify read count matches `samtools view -c` for specific genes
  - Check that all reads have valid gene assignments or are marked as intergenic
- **Validation metrics:**
  - Junction extraction accuracy: >99% match to ground truth on test set
  - Gene assignment rate: >85% of mapped reads assigned to genes
  - Coordinate consistency: 0 cases where `tss_pos > tes_pos` on + strand

---

### Phase 3: CNN Training Data Construction

#### Deliverables:
1. **Sequence window extractor** (`src/cnn/build_training_windows.py`)
   - Input: Priors, reference genome
   - Output:
     - `data/derived/cnn_train/tss_windows.npz` (one-hot encoded, labels)
     - `data/derived/cnn_train/tes_windows.npz`
     - `data/derived/cnn_train/donor_windows.npz`
     - `data/derived/cnn_train/acceptor_windows.npz`

2. **Training set statistics** (`qc/cnn_training_data.json`)
   - Positive/negative counts per task
   - GC content distributions
   - Chromosome splits (train/val/test)

3. **Data augmentation pipeline** (optional)
   - Reverse complement augmentation
   - Hard negative mining strategy

#### Testing Methods:
- **Unit tests:**
  - One-hot encoding correctness (check A=[1,0,0,0], etc.)
  - Window centering on exact prior coordinates
  - Negative sampling constraints (>1kb from positives)
  - Sequence extraction matches reference genome (spot checks)
- **Integration tests:**
  - GC content matching between positives and negatives
  - No label leakage between train/val/test chromosomes
  - Verify canonical splice motifs enriched in positive donor/acceptor examples
- **Validation metrics:**
  - TSS positives: >70% have TATA box or Inr motifs within window
  - TES positives: >60% have AATAAA or ATTAAA within 50bp
  - Donor positives: >95% have GT dinucleotide at junction
  - Acceptor positives: >95% have AG dinucleotide at junction
  - Class balance: negatives outnumber positives by <100:1

---

### Phase 4: CNN Model Training and Validation

#### Deliverables:
1. **Multi-task CNN implementation** (`src/cnn/model.py`)
   - Architecture definition
   - Training loop with early stopping
   - Calibration utilities

2. **Trained model weights** (`model_weights/cnn_multitask_v1.pt`)

3. **Model performance report** (`qc/cnn_performance.html`)
   - AUROC/AUPRC per task on held-out test set
   - Calibration curves (reliability diagrams)
   - Confusion matrices at multiple thresholds
   - Attention/gradient visualizations for motif discovery

4. **Inference script** (`src/cnn/score_positions.py`)
   - Score arbitrary BED positions with trained model

#### Testing Methods:
- **Unit tests:**
  - Model forward pass shape correctness
  - Loss computation (multi-task BCE)
  - Gradient flow (no vanishing/exploding)
- **Cross-validation:**
  - 5-fold CV on training chromosomes
  - Chromosome-holdout validation (test generalization to unseen chromosomes)
- **Validation metrics:**
  - TSS task: AUROC >0.90, AUPRC >0.75
  - TES task: AUROC >0.88, AUPRC >0.70
  - Donor task: AUROC >0.95, AUPRC >0.85
  - Acceptor task: AUROC >0.95, AUPRC >0.85
  - Calibration: expected calibration error (ECE) <0.05 for all tasks
- **Biological validation:**
  - Top-scoring TSS regions should overlap CAGE peaks
  - Top-scoring splice sites should match canonical GT-AG motifs
  - Model should assign low scores to frameshifting junctions

---

### Phase 5: Likelihood Model Implementation

#### Deliverables:
1. **Likelihood components** (`src/infer/likelihood.py`)
   - `compute_splice_likelihood(read, isoform, cnn_scores, junction_prior)`
   - `compute_tss_likelihood(position, cnn_score, tss_prior)`
   - `compute_tes_likelihood(position, cnn_score, tes_prior, internal_priming_score)`
   - `compute_qc_likelihood(read_metrics, platform_params)`

2. **Noise model definition** (`src/infer/noise_model.py`)
   - Uniform splice model
   - Broad end distributions
   - Low-quality read scoring

3. **Unit test suite** (`tests/test_likelihood.py`)

#### Testing Methods:
- **Unit tests:**
  - Likelihood values are probabilities (0-1 range after exp)
  - Log-likelihood is negative or zero
  - Known high-quality reads score higher than low-quality reads
  - Noise model assigns broad probability
- **Synthetic data tests:**
  - Generate synthetic reads from known isoforms
  - Verify correct isoform receives highest likelihood
  - Verify junk reads assigned to noise component
- **Validation metrics:**
  - Likelihood ordering: annotated isoform match >0.9 for clean reads
  - Noise detection: artifactual reads (low MAPQ + A-rich ends) >0.7 noise probability

---

### Phase 6: Per-Gene EM Inference

#### Deliverables:
1. **Candidate isoform generator** (`src/infer/candidates.py`)
   - Extract unique splice chains per gene from reads
   - Add annotated isoforms from GTF
   - Pruning logic for rare chains

2. **EM algorithm implementation** (`src/infer/em_inference.py`)
   - E-step: compute responsibilities
   - M-step: update mixture weights and end distributions
   - Convergence detection
   - Posterior computation

3. **Single-gene inference CLI** (`src/infer/infer_gene.py`)
   - Run inference on one gene for debugging

#### Testing Methods:
- **Unit tests:**
  - EM convergence on toy data (3 isoforms, 100 reads)
  - Responsibilities sum to 1 across isoforms + noise
  - Mixture weights sum to 1
  - Log-likelihood increases monotonically
- **Simulation tests:**
  - Generate reads from known mixture (e.g., 3 isoforms at 50%/30%/20%)
  - Verify EM recovers true abundances within 5%
  - Test noise component captures artifactual reads
- **Biological validation (single gene):**
  - Run on well-characterized gene (e.g., ACTB, GAPDH)
  - Compare inferred isoforms to MANE transcripts
  - Verify major isoforms have posterior >0.8
- **Validation metrics:**
  - EM convergence: <0.001 log-likelihood change in final 5 iterations
  - Abundance recovery: Pearson R >0.95 on simulated data
  - Isoform precision: annotated isoforms recovered with >0.7 posterior

---

### Phase 7: Genome-Wide Inference and Export

#### Deliverables:
1. **Parallel genome-wide inference** (`src/infer/infer_genome.py`)
   - Per-gene parallelization (multiprocessing or Snakemake)
   - Progress logging
   - Error handling for edge cases

2. **Output files** (`data/derived/infer_out/`)
   - `isoforms.gtf` (GTF format with posterior in score field)
   - `isoforms.tsv` (gene_id, isoform_id, splice_chain, posterior, abundance, tss_mode, tss_sd, tes_mode, tes_sd)
   - `read_assignments.tsv.gz` (read_id, gene_id, isoform_id, responsibility, noise_prob)
   - `qc_summary.json` (genes processed, isoforms inferred, total abundance, etc.)

3. **Summary report** (`qc/genome_inference_summary.html`)
   - Isoforms per gene distribution
   - Posterior distribution across all isoforms
   - Noise assignment rate
   - TSS/TES uncertainty distributions

#### Testing Methods:
- **Integration tests:**
  - Run on subset of genes (e.g., chr21 or 100 random genes)
  - Verify GTF format validity (check with `gtftools` or `gffutils`)
  - Check all reads have assignment probabilities
- **Reproducibility tests:**
  - Run inference twice with same parameters, verify identical results
  - Split-half reliability: correlate isoform posteriors from two random read subsamples
- **Validation metrics:**
  - Genome coverage: >90% of genes with ≥10 reads produce ≥1 isoform with posterior >0.5
  - Noise filtering: <10% of high-MAPQ reads assigned to noise
  - End uncertainty: TSS/TES standard deviation <100bp for high-posterior isoforms

---

### Phase 8: Evaluation and Benchmarking

#### Deliverables:
1. **Internal consistency metrics** (`src/eval/internal_validation.py`)
   - Replicate reproducibility (if multiple samples available)
   - Calibration analysis (posterior vs empirical support)

2. **External validation metrics** (`src/eval/external_validation.py`)
   - TSS distance to CAGE/refTSS peaks
   - TES distance to polyA sites
   - Junction agreement with GTEx
   - Comparison to MANE/GENCODE principal isoforms

3. **Baseline comparison** (`src/eval/compare_baselines.py`)
   - Run FLAIR collapse on same BAM
   - Run SQANTI3 filtering
   - Compare isoform counts, end accuracy, artifact rates

4. **Benchmark report** (`qc/benchmark_results.html`)
   - Side-by-side metrics table
   - Precision-recall curves
   - Artifact suppression rates

#### Testing Methods:
- **Internal validation:**
  - Posterior calibration: bin isoforms by posterior, measure recurrence rate in replicate
  - Expected: isoforms with posterior 0.9 recur in >80% of replicates
- **External validation:**
  - TSS accuracy: median distance to nearest CAGE peak
  - Expected: <100bp for high-posterior isoforms
  - TES accuracy: median distance to polyA site
  - Expected: <50bp for high-posterior isoforms
  - Junction concordance: % of inferred junctions in GTEx
  - Expected: >85% for junctions with CNN score >0.8
- **Baseline comparison:**
  - Internal priming rate: fraction of TES in A-rich regions without PAS
  - Expected: flAIr <5%, FLAIR ~15-25%
  - Singleton suppression: fraction of isoforms supported by 1 read
  - Expected: flAIr <20%, FLAIR ~40-60%

---

### Phase 9: Documentation and Packaging

#### Deliverables:
1. **User documentation** (`docs/`)
   - Installation instructions
   - Quick start tutorial
   - Full CLI reference
   - Input file format specifications
   - Output file descriptions
   - Hyperparameter tuning guide

2. **Example workflow** (`examples/tutorial.sh`)
   - End-to-end example on public data (e.g., NA12878 PacBio)

3. **Docker container** (`Dockerfile`)
   - Pre-packaged environment with all dependencies
   - Pre-trained CNN weights included

4. **Conda package** (optional)
   - Publish to bioconda

#### Testing Methods:
- **Documentation testing:**
  - Run all example commands in fresh environment
  - Verify tutorial completes without errors
- **Installation testing:**
  - Test pip install on clean Python 3.9, 3.10, 3.11 environments
  - Test Docker build and run
- **Continuous integration:**
  - GitHub Actions for unit tests on every commit
  - Integration tests on smaller test dataset

---

## 18) Testing Strategy Summary

### Unit Tests (per phase)
- Test individual functions in isolation
- Mock external dependencies (file I/O, genome access)
- Fast execution (<1 min total)
- Target: >80% code coverage

### Integration Tests
- Test workflows across multiple components
- Use small real datasets (e.g., single chromosome, subset of genes)
- Verify file format compatibility
- Target: All major workflows pass

### Validation Tests
- Biological correctness checks
- Compare against known biology (CAGE, polyA sites, annotated junctions)
- Use both simulated and real data
- Target: Meet all specified validation metrics

### Regression Tests
- Lock in results on benchmark datasets
- Detect unintended changes in model behavior
- Run on every release

### Performance Tests
- Runtime benchmarking (e.g., 100M reads should process in <X hours)
- Memory profiling (should handle 10GB BAM on 32GB RAM machine)
- Scalability testing (linear scaling with read count)

---

## 19) Success Criteria (Project Completion)

The project is considered successful when:

1. **Functionality:**
   - All 9 phases delivered with passing tests
   - Genome-wide inference completes on real PacBio/ONT dataset
   - Outputs valid GTF and TSV files

2. **Performance:**
   - CNN models achieve target AUROC/AUPRC metrics
   - EM inference converges for >95% of genes
   - Runtime: <4 hours for 50M reads on 16-core machine

3. **Accuracy:**
   - TSS median distance to CAGE peaks <100bp
   - TES median distance to polyA sites <50bp
   - Junction concordance with GTEx >85%

4. **Improvement over baselines:**
   - Internal priming rate reduced by >50% vs FLAIR
   - Singleton isoform rate reduced by >50% vs FLAIR
   - Posterior calibration: ECE <0.1

5. **Usability:**
   - Complete documentation published
   - One-command installation via pip or conda
   - Tutorial completable in <2 hours

6. **Reproducibility:**
   - Pre-trained models publicly available
   - All priors packaged and versioned
   - Example datasets with expected outputs

---

## Next action items (concrete)
1) Confirm your finalized file formats for:
- `tss_prior.bed.gz`
- `tes_prior.bed.gz`
- `junction_prior.tsv.gz`
2) Build training windows for CNN:
- positives from priors
- negatives sampled as described
3) Train multi-task CNN and export weights
4) Extract `reads.parquet` from your BAM
5) Implement per-gene EM inference using CNN+priors

```

