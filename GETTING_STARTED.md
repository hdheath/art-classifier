# Getting Started with flAIr

Welcome to the flAIr project! This guide will help you set up your development environment and begin Phase 1.

---

## What Has Been Set Up

### Documentation Created ✓
1. **[project_design.md](project_design.md)** - Complete technical specification with 9 implementation phases
2. **[DATA_VALIDATION_REPORT.md](DATA_VALIDATION_REPORT.md)** - Assessment of your existing TSS/TES/junction databases
3. **[GANTT_CHART.md](GANTT_CHART.md)** - 24-week timeline with deliverables and testing methods
4. **[README.md](README.md)** - Project overview and user documentation
5. **This file** - Quick start guide

### Directory Structure ✓
```
flAIr/
├── data/
│   ├── ref/                    # Reference genomes and priors (populate this)
│   ├── bam/                    # Input BAM files (add your BAMs here)
│   └── derived/                # Intermediate and output files (generated)
├── src/
│   ├── extract/                # Data extraction scripts
│   │   └── standardize_priors.py ✓  # Ready to use!
│   ├── cnn/                    # CNN training (to be implemented)
│   ├── infer/                  # EM inference (to be implemented)
│   └── eval/                   # Evaluation (to be implemented)
├── model_weights/              # CNN weights (will be populated)
├── qc/                         # QC reports (generated)
├── tests/                      # Unit tests (to be implemented)
├── docs/                       # Extended documentation
└── examples/                   # Tutorial workflows
```

### Initial Scripts ✓
- `src/extract/standardize_priors.py` - Ready to standardize your database files
- `requirements.txt` - Python dependencies
- `setup.py` - Package installation script

---

## Critical Findings from Data Assessment

### ✓ Good News:
1. **TES data is ready**: `polyadb.hg38.weighted.bed` is already in hg38 and good format
2. **TSS data looks good**: refTSS and FANTOM BED are hg38-compatible
3. **Splice junctions are comprehensive**: RJunBase has ~800K junctions

### ⚠️ Action Required:
1. **FANTOM CAGE file needs coordinate extraction**: The large .gz file contains expression data, needs parsing
2. **RJunBase genome build verification**: Need to confirm it's hg38 (likely is, based on coordinates)
3. **Missing reference files**: Need to download GRCh38 and GENCODE

---

## Phase 1: Week 1-2 Action Items

### Week 1 (This Week):

#### 1. Verify Genome Builds (Priority: HIGH)
```bash
# Check RJunBase coordinates against GENCODE
# Example: Verify RNF213 junction from your data
# chr17:80328477|80332006:+

# Should match GENCODE v46 RNF213 gene coordinates
```

#### 2. Download Required Files (Priority: HIGH)
```bash
cd /private/groups/brookslab/hdheath/projects/flAIr/data/ref/

# GRCh38 primary assembly
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
samtools faidx GRCh38.primary_assembly.genome.fa

# GENCODE v46 annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz
```

#### 3. Run Prior Standardization (Priority: MEDIUM)
```bash
cd /private/groups/brookslab/hdheath/projects/flAIr

# Install dependencies first
pip install pandas numpy pysam

# Run standardization
python src/extract/standardize_priors.py \
    --tss-fantom TSS_db/FANTOM_TSS_human.bed \
    --tss-reftss TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt \
    --tes TTS_db/polyadb.hg38.weighted.bed \
    --junction splice_db/detail_Linear_splice_annotation.txt \
    --output data/ref/ \
    --threads 4
```

**Expected outputs:**
- `data/ref/tss_prior.bed.gz` + `.tbi`
- `data/ref/tes_prior.bed.gz` + `.tbi`
- `data/ref/junction_prior.tsv.gz` + `.tbi`

#### 4. Create Gene Boundaries (Priority: MEDIUM)
```bash
# After downloading GENCODE, extract gene boundaries
# This script will be created in Week 2
```

---

### Week 2:

#### 1. Implement Gene Boundary Extraction
Create `src/extract/create_gene_loci.py`:
- Parse GENCODE GTF
- Extract gene-level intervals (strand-aware)
- Output `data/ref/genes.bed.gz`

#### 2. QC Validation
Create `src/extract/validate_priors.py`:
- Check TSS overlap with GENCODE gene starts (>80% within 500bp)
- Check junction overlap with GENCODE junctions (>90%)
- Verify motif enrichment (GT-AG for splice donors/acceptors)
- Generate QC report: `qc/prior_summary.html`

#### 3. Test on Sample Data
- Find a small test BAM (e.g., chr21 only)
- Verify all priors are queryable with tabix
- Spot-check coordinates in IGV

---

## Quick Commands Reference

### Examine Your Current Data
```bash
cd /private/groups/brookslab/hdheath/projects/flAIr

# Count TSS peaks
wc -l TSS_db/*.bed TSS_db/*.txt

# Count TES sites
wc -l TTS_db/polyadb.hg38.weighted.bed

# Count junctions
wc -l splice_db/detail_Linear_splice_annotation.txt

# Check chromosome naming
head -5 TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt
head -5 TTS_db/polyadb.hg38.weighted.bed
```

### Validate Tabix Indices (After Standardization)
```bash
# Query a region
tabix data/ref/tss_prior.bed.gz chr1:1000000-2000000
tabix data/ref/tes_prior.bed.gz chr1:1000000-2000000
tabix data/ref/junction_prior.tsv.gz chr1:1000000-2000000
```

### Check Genome Build
```bash
# Compare coordinates
# Get gene coordinates from GENCODE
zcat data/ref/gencode.v46.annotation.gtf.gz | \
    awk '$3=="gene" && $10=="\"RNF213\""' | head -1

# Should match junction:
# chr17:80328477|80332006:+ (from your data)
# chr17 29584646 29859590 (GENCODE gene span)
```

---

## Development Workflow

### 1. Install Development Environment
```bash
cd /private/groups/brookslab/hdheath/projects/flAIr

# Create conda environment
conda create -n flair python=3.11
conda activate flair

# Install dependencies
pip install -r requirements.txt

# Install in development mode
pip install -e .
```

### 2. Run Tests (Once Implemented)
```bash
pytest tests/ -v
```

### 3. Format Code
```bash
black src/
flake8 src/
```

---

## Expected Timeline

### Phase 1 (Weeks 1-3): Data Prep
- **Week 1**: Download files, verify builds, run standardization
- **Week 2**: Gene boundaries, QC validation
- **Week 3**: Finalize all priors, comprehensive testing

### Next Milestone
**Phase 1 Complete** when:
- [ ] All priors in hg38, bgzipped, tabix-indexed
- [ ] TSS coverage >80% of GENCODE genes
- [ ] Junction coverage >90% of GENCODE junctions
- [ ] Motif enrichment validated (GT-AG >95%)
- [ ] QC report generated

---

## Troubleshooting

### Common Issues:

#### "tabix not found"
```bash
# Install samtools/tabix
conda install -c bioconda samtools htslib
```

#### "pandas not found"
```bash
pip install pandas numpy
```

#### "Permission denied on script"
```bash
chmod +x src/extract/standardize_priors.py
```

#### "Coordinate mismatch"
Check genome build:
```bash
# All files should use same chromosome naming
# Either: chr1, chr2, ... (UCSC)
# Or: 1, 2, ... (Ensembl)

# Your data appears to use UCSC (chr prefix) which matches GRCh38
```

---

## Getting Help

### Resources:
1. **Project design**: See [project_design.md](project_design.md) for full technical specs
2. **Timeline**: See [GANTT_CHART.md](GANTT_CHART.md) for detailed schedule
3. **Data validation**: See [DATA_VALIDATION_REPORT.md](DATA_VALIDATION_REPORT.md) for current state

### Questions to Resolve:
1. Do you have a test BAM file ready? (Ideally <10GB for development)
2. What compute resources are available? (GPU for CNN training in Phase 4)
3. Target dataset for final evaluation? (Needed for Phase 8)

---

## Success Criteria for Week 1

By end of Week 1, you should have:
- [x] Project structure set up ✓ (Done)
- [ ] GRCh38 reference downloaded
- [ ] GENCODE v46 GTF downloaded
- [ ] RJunBase genome build verified
- [ ] Standardization script run successfully
- [ ] All three prior files (.bed.gz/.tsv.gz) with tabix indices

---

## Next Steps After Week 1

1. **Implement gene boundary extraction** (src/extract/create_gene_loci.py)
2. **Create QC validation script** (src/extract/validate_priors.py)
3. **Begin Phase 2 planning**: BAM extraction pipeline design
4. **Start Phase 3 in parallel**: CNN training data construction

---

## Quick Start Summary

```bash
# 1. Go to project directory
cd /private/groups/brookslab/hdheath/projects/flAIr

# 2. Download reference files (Week 1)
cd data/ref
wget [GRCh38 URL]
wget [GENCODE URL]

# 3. Run standardization (Week 1)
cd ../..
python src/extract/standardize_priors.py [options]

# 4. Validate (Week 2)
python src/extract/validate_priors.py --priors data/ref/

# 5. Extract gene boundaries (Week 2)
python src/extract/create_gene_loci.py \
    --gtf data/ref/gencode.v46.annotation.gtf.gz \
    --output data/ref/genes.bed.gz

# 6. Ready for Phase 2!
```

---

Good luck! This is an ambitious and well-designed project. Take it one phase at a time, and don't hesitate to iterate on the design as you learn more from the data.

**Remember**: The goal is not just to build a tool, but to improve transcriptome inference accuracy with principled probabilistic methods. Focus on validation and testing at each phase!
