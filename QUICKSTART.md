# flAIr Quick Start Guide

## One-Command Setup

```bash
cd /private/groups/brookslab/hdheath/projects/flAIr
./setup_environment.sh
```

This automated script will:
1. Create the conda environment with GPU support
2. Install all dependencies
3. Create package structure
4. Install flAIr in development mode
5. Verify everything works

**Time:** 5-10 minutes

---

## Manual Setup (Alternative)

If you prefer manual control:

```bash
# 1. Check your CUDA version first
nvidia-smi

# 2. Edit environment.yml if needed (line 29)
# - CUDA 11.x → keep pytorch-cuda=11.8
# - CUDA 12.x → change to pytorch-cuda=12.1
# - No GPU → remove pytorch-cuda line

# 3. Create environment
conda env create -f environment.yml

# 4. Activate
conda activate flair

# 5. Create package structure
touch src/__init__.py src/extract/__init__.py \
      src/cnn/__init__.py src/infer/__init__.py \
      src/eval/__init__.py tests/__init__.py

# 6. Install package
pip install -e .

# 7. Verify
python -c "import torch; print(f'CUDA: {torch.cuda.is_available()}')"
```

---

## Your Data Status

You already have excellent data ready to use:

| Data Type | File | Status | Count |
|-----------|------|--------|-------|
| **Splice Junctions** | [detail_Linear_splice_annotation.txt](splice_db/detail_Linear_splice_annotation.txt) | ✓ Ready | 682K |
| **TSS (FANTOM)** | [FANTOM_TSS_human.bed](TSS_db/FANTOM_TSS_human.bed) | ✓ Ready | ~100K |
| **TSS (refTSS)** | [refTSS_v4.1_human_coordinate.hg38.bed.txt](TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt) | ✓ Ready | ~60K |
| **TES (polyA)** | [polyadb.hg38.weighted.bed](TTS_db/polyadb.hg38.weighted.bed) | ✓ Ready | ~70K |

**No additional downloads needed for training data!**

---

## Week 1 Checklist

### Day 1: Environment Setup
- [ ] Run `./setup_environment.sh` OR follow manual setup
- [ ] Verify GPU detection: `python -c "import torch; print(torch.cuda.is_available())"`
- [ ] Test script: `python src/extract/standardize_priors.py --help`

### Day 2-3: Process Prior Databases
```bash
conda activate flair

python src/extract/standardize_priors.py \
    --tss-fantom TSS_db/FANTOM_TSS_human.bed \
    --tss-reftss TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt \
    --tes TTS_db/polyadb.hg38.weighted.bed \
    --junction splice_db/detail_Linear_splice_annotation.txt \
    --output data/ref/
```

**Expected output:**
- `data/ref/tss_prior.bed.gz` + `.tbi` (~150K peaks)
- `data/ref/tes_prior.bed.gz` + `.tbi` (~70K sites)
- `data/ref/junction_prior.tsv.gz` + `.tbi` (~650K junctions after filtering)

**Time:** 10-15 minutes

### Day 4: Download Reference Genome

```bash
# Create directory
mkdir -p data/ref/

# Download GRCh38 (3.1GB)
wget -P data/ref/ \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz

# Decompress
gunzip data/ref/GRCh38.primary_assembly.genome.fa.gz

# Index with samtools
samtools faidx data/ref/GRCh38.primary_assembly.genome.fa
```

### Day 5: Download GENCODE Annotation

```bash
# Download GENCODE GTF (54MB compressed)
wget -P data/ref/ \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz

# Index with tabix
tabix -p gff data/ref/gencode.v46.annotation.gtf.gz
```

### Verify Week 1 Completion

```bash
# Check all files exist
ls -lh data/ref/*.gz data/ref/*.tbi data/ref/*.fai

# Test tabix queries
tabix data/ref/tss_prior.bed.gz chr1:1000000-2000000
tabix data/ref/junction_prior.tsv.gz chr1:1000000-2000000

# Verify genome
samtools faidx data/ref/GRCh38.primary_assembly.genome.fa chr1:1000000-1000100
```

---

## Common Issues

### "conda: command not found"
```bash
~/miniconda3/bin/conda init bash
exec bash
```

### "CUDA available: False" (but you have GPU)
```bash
# Check CUDA version
nvidia-smi

# Reinstall PyTorch with correct CUDA
conda activate flair
conda remove pytorch torchvision pytorch-cuda
conda install pytorch torchvision pytorch-cuda=11.8 -c pytorch -c nvidia
```

### "No module named 'extract'"
```bash
# Make sure __init__.py files exist
touch src/__init__.py src/extract/__init__.py

# Reinstall package
pip install -e .
```

---

## What's Next?

After Week 1, you'll be ready for:

**Week 2-3: Phase 2 - BAM Read Extraction**
- Process your long-read BAM files
- Extract aligned segments by gene
- Save to Parquet format

**Week 4-6: Phase 3 - CNN Training Data Generation**
- Build training windows from prior databases
- Generate positive and negative examples
- Create balanced datasets

**Week 7-10: Phase 4 - CNN Training**
- Train multi-task CNN for TSS/TES/splice sites
- Validate on held-out data
- Export trained model weights

See [GANTT_CHART.md](GANTT_CHART.md) for full timeline.

---

## Key Files Reference

| File | Purpose |
|------|---------|
| [README.md](README.md) | Project overview |
| [SETUP_GUIDE.md](SETUP_GUIDE.md) | Detailed setup instructions |
| [GETTING_STARTED.md](GETTING_STARTED.md) | Week 1 walkthrough |
| [project_design.md](project_design.md) | Complete technical design |
| [GANTT_CHART.md](GANTT_CHART.md) | 24-week timeline |
| [DATA_VALIDATION_REPORT.md](DATA_VALIDATION_REPORT.md) | Your database assessment |
| [TSS_DATA_GUIDE.md](TSS_DATA_GUIDE.md) | TSS data recommendations |
| [RJUNBASE_ANALYSIS.md](RJUNBASE_ANALYSIS.md) | Splice junction analysis |

---

## Getting Help

- **Environment issues**: See [SETUP_GUIDE.md](SETUP_GUIDE.md)
- **Data questions**: See [DATA_VALIDATION_REPORT.md](DATA_VALIDATION_REPORT.md)
- **Implementation details**: See [project_design.md](project_design.md)
- **Timeline questions**: See [GANTT_CHART.md](GANTT_CHART.md)

---

## Summary

**You have everything you need to start!**

1. ✓ Excellent training data (682K junctions, 150K TSS, 70K TES)
2. ✓ GPU access for CNN training
3. ✓ Automated setup script
4. ✓ Clear 24-week plan with deliverables

**First command to run:**
```bash
./setup_environment.sh
```

**Total setup time:** ~20-30 minutes including downloads

Good luck! 🚀
