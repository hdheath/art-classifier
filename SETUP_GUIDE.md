# flAIr Environment Setup Guide

## Quick Start (Recommended)

### Option 1: Create conda environment from environment.yml

```bash
# Navigate to project directory
cd /private/groups/brookslab/hdheath/projects/flAIr

# Create conda environment
conda env create -f environment.yml

# Activate environment
conda activate flair

# Install flAIr package in development mode
pip install -e .
```

**Time:** 5-10 minutes
**Result:** Fully configured environment ready to use

---

### Option 2: Manual conda environment setup

```bash
# Create new environment with Python 3.11
conda create -n flair python=3.11 -y

# Activate
conda activate flair

# Install dependencies
conda install -c conda-forge -c bioconda \
    numpy pandas scipy pyarrow \
    pysam samtools bedtools tabix \
    pytorch torchvision \
    matplotlib seaborn plotly \
    tqdm pyyaml pytest black flake8 \
    ipython jupyter -y

# Install flAIr in development mode
pip install -e .
```

---

## What Each Component Does

### Core Dependencies:

| Package | Purpose | Phase Used |
|---------|---------|------------|
| **numpy, pandas** | Data manipulation | All phases |
| **scipy** | EM algorithm, statistics | Phase 6-7 |
| **pyarrow** | Fast parquet I/O | Phase 2 |
| **pysam** | BAM file reading | Phase 2 |
| **pytorch** | CNN training | Phase 3-4 |
| **matplotlib, seaborn** | Visualization, QC reports | All phases |

### Genomics Tools:

| Tool | Purpose | Phase Used |
|------|---------|------------|
| **samtools** | BAM indexing, stats | Phase 2 |
| **bedtools** | Interval operations | Phase 1 |
| **tabix** | Index BED/TSV files | Phase 1 |

---

## GPU Support (Included by Default)

The environment.yml includes GPU support by default with CUDA 11.8.

### Verify GPU is available after setup:
```bash
conda activate flair
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"
python -c "import torch; print(f'GPU count: {torch.cuda.device_count()}')"
```

### If you need a different CUDA version:

**Before creating the environment**, check your CUDA version:
```bash
nvidia-smi  # Look for "CUDA Version: X.X" in the output
```

Then edit environment.yml:
```bash
# For CUDA 12.1 (if nvidia-smi shows 12.x)
# Change line 29 from:
#   - pytorch-cuda=11.8
# to:
#   - pytorch-cuda=12.1

# For CPU-only (no GPU):
# Remove line 29 entirely:
#   - pytorch-cuda=11.8
```

---

## Verify Installation

After setup, test that everything works:

```bash
# Activate environment
conda activate flair

# Test Python packages
python -c "
import numpy as np
import pandas as pd
import pysam
import torch
print('NumPy:', np.__version__)
print('Pandas:', pd.__version__)
print('PyTorch:', torch.__version__)
print('Pysam:', pysam.__version__)
print('All imports successful!')
"

# Test command-line tools
samtools --version
bedtools --version
tabix --version

# Test flAIr package
python -c "
import sys
sys.path.insert(0, 'src')
# Basic import test - will fail until we create __init__.py files
"
```

---

## Current Environment Issues

You tried to run:
```bash
python setup.py
```

**Problem:** `setup.py` is not meant to be run directly. It's used by pip during installation.

**Solution:**
```bash
# Instead, use pip to install in "editable" mode:
pip install -e .

# This runs setup.py automatically and installs the package
# The -e flag means changes to src/ files are immediately reflected
```

---

## File Structure After Setup

```
flair/  (conda environment - in ~/miniconda3/envs/flair/)
├── bin/
│   ├── python
│   ├── samtools
│   ├── bedtools
│   └── tabix
├── lib/
│   └── python3.11/site-packages/
│       ├── numpy/
│       ├── pandas/
│       ├── torch/
│       └── ...
└── ...

flAIr/  (project directory - current location)
├── src/
│   ├── extract/
│   │   ├── __init__.py  (need to create)
│   │   ├── standardize_priors.py
│   │   └── ...
│   ├── cnn/
│   ├── infer/
│   └── eval/
├── setup.py
├── environment.yml
└── ...
```

---

## Creating Missing __init__.py Files

For the package to work properly, we need `__init__.py` files:

```bash
cd /private/groups/brookslab/hdheath/projects/flAIr

# Create __init__.py files
touch src/__init__.py
touch src/extract/__init__.py
touch src/cnn/__init__.py
touch src/infer/__init__.py
touch src/eval/__init__.py
touch tests/__init__.py

# Now install package
pip install -e .
```

---

## Testing Your Scripts

After setup, test the standardization script:

```bash
# Activate environment
conda activate flair

# Test help message
python src/extract/standardize_priors.py --help

# Run on your data
python src/extract/standardize_priors.py \
    --tss-fantom TSS_db/FANTOM_TSS_human.bed \
    --tss-reftss TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt \
    --tes TTS_db/polyadb.hg38.weighted.bed \
    --junction splice_db/detail_Linear_splice_annotation.txt \
    --output data/ref/
```

---

## Common Issues and Solutions

### Issue 1: "conda: command not found"

**Solution:**
```bash
# Initialize conda for your shell
~/miniconda3/bin/conda init bash
# Or if using zsh:
~/miniconda3/bin/conda init zsh

# Restart shell
exec bash
```

### Issue 2: "Solving environment: failed"

**Solution:**
```bash
# Update conda
conda update -n base conda

# Try mamba (faster solver)
conda install -n base mamba -y
mamba env create -f environment.yml
```

### Issue 3: "No module named 'extract'"

**Solution:**
```bash
# Make sure __init__.py files exist
touch src/__init__.py src/extract/__init__.py

# Reinstall in editable mode
pip install -e .

# Or run scripts directly without package install
python src/extract/standardize_priors.py --help
```

### Issue 4: PyTorch CPU vs GPU mismatch

**Solution:**
```bash
# Uninstall current pytorch
conda remove pytorch torchvision

# Install specific version
# For CPU:
conda install pytorch torchvision cpuonly -c pytorch

# For GPU (CUDA 11.8):
conda install pytorch torchvision pytorch-cuda=11.8 -c pytorch -c nvidia
```

---

## Development Workflow

### Daily usage:

```bash
# Activate environment (do this every time)
conda activate flair

# Your environment should show in prompt:
(flair) hdheath@mustard:~/projects/flAIr$

# Run scripts
python src/extract/standardize_priors.py ...

# Run tests
pytest tests/

# Deactivate when done
conda deactivate
```

### Installing new packages:

```bash
# Activate environment
conda activate flair

# Install from conda
conda install package_name

# Or from pip
pip install package_name

# Update environment.yml to remember it
conda env export > environment.yml
```

---

## Recommended Workflow for Week 1

1. **Setup environment** (5-10 mins):
   ```bash
   cd /private/groups/brookslab/hdheath/projects/flAIr
   conda env create -f environment.yml
   conda activate flair
   ```

2. **Create __init__.py files** (1 min):
   ```bash
   touch src/__init__.py src/extract/__init__.py \
         src/cnn/__init__.py src/infer/__init__.py \
         src/eval/__init__.py tests/__init__.py
   ```

3. **Install package** (1 min):
   ```bash
   pip install -e .
   ```

4. **Test installation** (1 min):
   ```bash
   python src/extract/standardize_priors.py --help
   samtools --version
   tabix --version
   ```

5. **Process your data** (5-10 mins):
   ```bash
   python src/extract/standardize_priors.py \
       --tss-fantom TSS_db/FANTOM_TSS_human.bed \
       --tss-reftss TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt \
       --tes TTS_db/polyadb.hg38.weighted.bed \
       --junction splice_db/detail_Linear_splice_annotation.txt \
       --output data/ref/
   ```

**Total time: ~20-30 minutes**

---

## Summary

**Your current error:**
```
python setup.py
error: no commands supplied
```

**What you should do instead:**
```bash
# Step 1: Create environment
conda env create -f environment.yml

# Step 2: Activate
conda activate flair

# Step 3: Install package
pip install -e .

# Step 4: Use the scripts
python src/extract/standardize_priors.py --help
```

**That's it!** You'll have a fully functional development environment.
