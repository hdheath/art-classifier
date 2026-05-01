# flAIr Project Gantt Chart
## Generative Transcriptome Inference with CNN Priors

**Project Start:** Week 1 (Current)
**Estimated Completion:** Week 24-28

---

## Timeline Overview

```
Phase 1: Data Prep         [████░░░░░░] Weeks  1-3
Phase 2: BAM Extraction    [░░░░████░░] Weeks  3-5
Phase 3: CNN Training Data [░░░░░░██░░] Weeks  4-6
Phase 4: CNN Training      [░░░░░░░███] Weeks  6-9
Phase 5: Likelihood Model  [░░░░░░░░██] Weeks  8-10
Phase 6: EM Inference      [░░░░░░░░░█] Weeks 10-14
Phase 7: Genome-wide       [░░░░░░░░░░] Weeks 14-18
Phase 8: Evaluation        [░░░░░░░░░░] Weeks 16-22
Phase 9: Documentation     [░░░░░░░░░░] Weeks 20-24
```

---

## Phase 1: Data Preparation and Validation (Weeks 1-3)

**Duration:** 3 weeks
**Dependencies:** None
**Effort:** 2-3 days/week

### Week 1: Setup and Inventory
- [x] Project design documentation
- [ ] Assess existing database files
- [ ] Verify genome builds (hg19 vs hg38)
- [ ] Download GRCh38 reference genome
- [ ] Download GENCODE v46 GTF

**Deliverable:** DATA_VALIDATION_REPORT.md ✓

### Week 2: Prior File Standardization
- [ ] Write TSS standardization script
- [ ] Write TES standardization script
- [ ] Write splice junction parser
- [ ] Generate standardized BED/TSV files
- [ ] Bgzip and tabix index all files

**Deliverable:** `src/extract/standardize_priors.py`

### Week 3: Gene Boundaries and QC
- [ ] Extract gene boundaries from GENCODE
- [ ] Generate QC report (coverage, overlaps)
- [ ] Validate chromosome consistency
- [ ] Check motif enrichment (GT-AG for splice sites)

**Deliverables:**
- `src/extract/create_gene_loci.py`
- `qc/prior_summary.html`
- All standardized prior files in `data/ref/`

**Milestone 1:** ✓ All priors standardized, indexed, and validated

---

## Phase 2: BAM Evidence Extraction (Weeks 3-5)

**Duration:** 2-3 weeks
**Dependencies:** Phase 1 (gene boundaries)
**Effort:** 3-4 days/week

### Week 3-4: Read Extraction Pipeline
- [ ] Implement CIGAR parsing (junctions, softclips)
- [ ] Strand-aware coordinate extraction
- [ ] Gene assignment logic (overlap with gene boundaries)
- [ ] QC metrics collection (MAPQ, NM, alignment length)

**Deliverable:** `src/extract/bam_to_reads.py`

### Week 4-5: BAM Processing and Validation
- [ ] Process sample BAM file
- [ ] Generate `reads.parquet`
- [ ] Create read QC report
- [ ] Validate junction extraction (IGV spot checks)
- [ ] Profile memory usage and runtime

**Deliverables:**
- `data/derived/reads.parquet`
- `qc/read_extraction_summary.html`

**Milestone 2:** ✓ Read evidence table extracted with QC passing

---

## Phase 3: CNN Training Data Construction (Weeks 4-6)

**Duration:** 3 weeks
**Dependencies:** Phase 1 (priors, genome)
**Effort:** 4-5 days/week
**Can start:** Parallel with Phase 2

### Week 4-5: Sequence Window Extraction
- [ ] Implement one-hot encoding
- [ ] Extract positive windows from priors
- [ ] Implement negative sampling strategy
- [ ] GC-content matched negatives
- [ ] Train/val/test chromosome splits

**Deliverable:** `src/cnn/build_training_windows.py`

### Week 5-6: Training Set Generation and QC
- [ ] Generate all window NPZ files
- [ ] Validate motif enrichment (TATA, polyA, GT-AG)
- [ ] Check class balance
- [ ] Generate training data statistics

**Deliverables:**
- `data/derived/cnn_train/tss_windows.npz`
- `data/derived/cnn_train/tes_windows.npz`
- `data/derived/cnn_train/donor_windows.npz`
- `data/derived/cnn_train/acceptor_windows.npz`
- `qc/cnn_training_data.json`

**Milestone 3:** ✓ Training data ready with biological validation passing

---

## Phase 4: CNN Model Training and Validation (Weeks 6-9)

**Duration:** 4 weeks
**Dependencies:** Phase 3 (training data)
**Effort:** Full time (5 days/week)

### Week 6-7: Model Implementation
- [ ] Design multi-task CNN architecture
- [ ] Implement training loop with early stopping
- [ ] Set up cross-validation framework
- [ ] Configure GPU environment

**Deliverable:** `src/cnn/model.py`

### Week 7-8: Model Training
- [ ] Train on full dataset
- [ ] Monitor AUROC/AUPRC metrics
- [ ] Tune hyperparameters
- [ ] Run 5-fold cross-validation
- [ ] Chromosome-holdout validation

**Target Metrics:**
- TSS: AUROC >0.90, AUPRC >0.75
- TES: AUROC >0.88, AUPRC >0.70
- Donor: AUROC >0.95, AUPRC >0.85
- Acceptor: AUROC >0.95, AUPRC >0.85

### Week 8-9: Model Validation and Export
- [ ] Calibration analysis (reliability diagrams)
- [ ] Gradient visualization for motif discovery
- [ ] Biological validation (CAGE overlap, GT-AG motifs)
- [ ] Export trained weights

**Deliverables:**
- `model_weights/cnn_multitask_v1.pt`
- `qc/cnn_performance.html`
- `src/cnn/score_positions.py`

**Milestone 4:** ✓ CNN models trained and validated, AUROC targets met

---

## Phase 5: Likelihood Model Implementation (Weeks 8-10)

**Duration:** 2-3 weeks
**Dependencies:** Phase 4 (CNN weights)
**Effort:** 3-4 days/week
**Can start:** Parallel with Phase 4 validation

### Week 8-9: Likelihood Components
- [ ] Implement splice likelihood function
- [ ] Implement TSS likelihood function
- [ ] Implement TES likelihood function
- [ ] Implement QC likelihood function
- [ ] Implement noise model

**Deliverables:**
- `src/infer/likelihood.py`
- `src/infer/noise_model.py`

### Week 9-10: Testing and Validation
- [ ] Unit tests for all likelihood functions
- [ ] Synthetic data generation
- [ ] Test on synthetic reads
- [ ] Validate likelihood ordering
- [ ] Validate noise detection

**Deliverable:** `tests/test_likelihood.py`

**Milestone 5:** ✓ Likelihood model implemented with synthetic tests passing

---

## Phase 6: Per-Gene EM Inference (Weeks 10-14)

**Duration:** 4 weeks
**Dependencies:** Phase 5 (likelihood), Phase 2 (reads)
**Effort:** Full time

### Week 10-11: EM Algorithm Implementation
- [ ] Candidate isoform generator
- [ ] E-step implementation
- [ ] M-step implementation
- [ ] Convergence detection
- [ ] Posterior computation

**Deliverables:**
- `src/infer/candidates.py`
- `src/infer/em_inference.py`

### Week 11-12: Single-Gene Testing
- [ ] Implement single-gene CLI
- [ ] Test on well-characterized genes (ACTB, GAPDH)
- [ ] Verify EM convergence
- [ ] Compare to GENCODE/MANE isoforms
- [ ] Profile runtime

**Deliverable:** `src/infer/infer_gene.py`

### Week 12-14: Simulation and Validation
- [ ] Generate simulated reads from known mixtures
- [ ] Test abundance recovery (target: R >0.95)
- [ ] Test noise component performance
- [ ] Optimize hyperparameters
- [ ] Debug edge cases

**Validation Targets:**
- EM convergence: <0.001 log-likelihood change
- Abundance recovery: Pearson R >0.95
- Isoform precision: >0.7 posterior for annotated isoforms

**Milestone 6:** ✓ EM inference working on single genes, validation targets met

---

## Phase 7: Genome-Wide Inference and Export (Weeks 14-18)

**Duration:** 4 weeks
**Dependencies:** Phase 6 (EM)
**Effort:** Full time

### Week 14-15: Parallelization
- [ ] Implement per-gene parallelization
- [ ] Set up progress logging
- [ ] Error handling for edge cases
- [ ] Memory optimization

**Deliverable:** `src/infer/infer_genome.py`

### Week 15-16: Pilot Run (Chromosome 21)
- [ ] Run on chr21 genes (~200 genes)
- [ ] Validate output formats
- [ ] Profile runtime and memory
- [ ] Debug issues

### Week 16-18: Full Genome Run
- [ ] Run on all chromosomes
- [ ] Generate all output files
- [ ] Create summary report
- [ ] Reproducibility testing

**Deliverables:**
- `data/derived/infer_out/isoforms.gtf`
- `data/derived/infer_out/isoforms.tsv`
- `data/derived/infer_out/read_assignments.tsv.gz`
- `data/derived/infer_out/qc_summary.json`
- `qc/genome_inference_summary.html`

**Target Performance:**
- Runtime: <4 hours for 50M reads (16 cores)
- Memory: <32GB RAM
- Coverage: >90% of genes with ≥10 reads produce isoforms

**Milestone 7:** ✓ Genome-wide inference complete, outputs generated

---

## Phase 8: Evaluation and Benchmarking (Weeks 16-22)

**Duration:** 6 weeks
**Dependencies:** Phase 7 (outputs)
**Effort:** 3-4 days/week
**Can start:** Week 16 (parallel with Phase 7 completion)

### Week 16-18: Internal Validation
- [ ] Replicate reproducibility analysis
- [ ] Posterior calibration analysis
- [ ] Split-half reliability testing
- [ ] TSS/TES distance to priors

**Deliverable:** `src/eval/internal_validation.py`

### Week 18-20: External Validation
- [ ] Compare TSS to CAGE peaks
- [ ] Compare TES to polyA atlas
- [ ] Junction concordance with GTEx
- [ ] Compare to MANE/GENCODE principal isoforms

**Deliverable:** `src/eval/external_validation.py`

### Week 20-22: Baseline Comparison
- [ ] Run FLAIR on same BAM
- [ ] Run SQANTI3 filtering
- [ ] Compare isoform counts
- [ ] Compare artifact rates
- [ ] Generate benchmark report

**Deliverables:**
- `src/eval/compare_baselines.py`
- `qc/benchmark_results.html`

**Validation Targets:**
- TSS accuracy: <100bp median distance to CAGE
- TES accuracy: <50bp median distance to polyA
- Junction concordance: >85% match to GTEx
- Internal priming: <5% (vs FLAIR ~15-25%)
- Singleton rate: <20% (vs FLAIR ~40-60%)

**Milestone 8:** ✓ Evaluation complete, improvement over baselines demonstrated

---

## Phase 9: Documentation and Packaging (Weeks 20-24)

**Duration:** 4 weeks
**Dependencies:** All phases
**Effort:** 3 days/week
**Can start:** Week 20 (parallel with Phase 8)

### Week 20-21: User Documentation
- [ ] Installation instructions
- [ ] Quick start tutorial
- [ ] Full CLI reference
- [ ] Input/output format specs
- [ ] Hyperparameter tuning guide

**Deliverable:** `docs/` directory

### Week 21-22: Example Workflow
- [ ] Download public dataset (NA12878)
- [ ] Create end-to-end tutorial
- [ ] Test tutorial in clean environment
- [ ] Record expected outputs

**Deliverable:** `examples/tutorial.sh`

### Week 22-23: Containerization
- [ ] Create Dockerfile
- [ ] Include pre-trained weights
- [ ] Test Docker build
- [ ] Create conda environment file
- [ ] Set up GitHub Actions CI

**Deliverables:**
- `Dockerfile`
- `environment.yml`
- `.github/workflows/tests.yml`

### Week 23-24: Release Preparation
- [ ] Final testing on multiple datasets
- [ ] Performance benchmarking
- [ ] Code cleanup and documentation
- [ ] Create release notes
- [ ] Tag v1.0.0 release

**Milestone 9:** ✓ Project complete, ready for publication and distribution

---

## Critical Path

The longest dependency chain (critical path):

```
Phase 1 → Phase 2 → Phase 6 → Phase 7 → Phase 8
(Week 1-3) (3-5)    (10-14)   (14-18)   (16-22)
```

**Critical Path Duration:** ~22 weeks

**Parallelizable work:**
- Phase 3 (CNN data) can start during Phase 2
- Phase 4 (CNN training) can start during Phase 2
- Phase 5 (Likelihood) can start during Phase 4
- Phase 8 (Evaluation) can start during Phase 7
- Phase 9 (Documentation) can start during Phase 8

---

## Risk Mitigation

### High-Risk Items:
1. **CNN not achieving target metrics** (Phase 4)
   - Mitigation: Start early (Week 6), iterate on architecture
   - Fallback: Use pre-trained sequence models (SpliceAI, Enformer)

2. **EM not converging** (Phase 6)
   - Mitigation: Extensive synthetic testing
   - Fallback: Use simpler assignment (hard assignment + posterior)

3. **Memory issues genome-wide** (Phase 7)
   - Mitigation: Profile early on chr21
   - Fallback: Process genes in batches, disk-based intermediate storage

4. **Poor baseline comparison** (Phase 8)
   - Mitigation: Run baseline early (Week 18)
   - Response: Iterate on noise model, hyperparameters

---

## Resource Requirements

### Computational:
- **CPU:** 16+ cores for parallelization
- **RAM:** 32-64GB for genome-wide inference
- **GPU:** 1x GPU (V100/A100) for CNN training (Weeks 6-9)
- **Storage:** 500GB for data, models, outputs

### Time Allocation:
- **Weeks 1-10:** Part-time (2-4 days/week) = ~60 days
- **Weeks 10-18:** Full-time = ~40 days
- **Weeks 18-24:** Part-time (3 days/week) = ~20 days
- **Total:** ~120 days (~6 months calendar time)

---

## Success Metrics by Phase

| Phase | Success Criteria |
|-------|------------------|
| 1 | All priors hg38, indexed, >80% GENCODE coverage |
| 2 | >85% reads assigned to genes, >99% junction accuracy |
| 3 | >95% splice sites have GT-AG, class balance <100:1 |
| 4 | AUROC >0.90 (TSS), >0.95 (splice), ECE <0.05 |
| 5 | Likelihood correctly orders synthetic reads |
| 6 | Abundance recovery R >0.95 on simulated data |
| 7 | <4 hours runtime for 50M reads, >90% gene coverage |
| 8 | TSS <100bp to CAGE, internal priming <5% |
| 9 | Tutorial completes without errors in clean env |

---

## Current Status

**Week:** 1
**Phase:** 1 (Data Preparation)
**Progress:** 10%

**Completed:**
- [x] Project design document
- [x] Data inventory assessment
- [x] Validation report

**Next Steps (Week 1-2):**
1. Verify RJunBase genome build (hg19 vs hg38)
2. Download GRCh38 reference and GENCODE v46
3. Write prior standardization scripts
4. Generate standardized BED/TSV files
