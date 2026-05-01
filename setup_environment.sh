#!/bin/bash
# flAIr Environment Setup Script
# This script automates the conda environment setup process

set -e  # Exit on error

echo "=========================================="
echo "flAIr Environment Setup"
echo "=========================================="
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found!"
    echo "Please initialize conda first:"
    echo "  ~/miniconda3/bin/conda init bash"
    echo "  exec bash"
    exit 1
fi

# Check CUDA version if nvidia-smi is available
if command -v nvidia-smi &> /dev/null; then
    echo "Detected GPU system. Checking CUDA version..."
    nvidia-smi | grep "CUDA Version" || echo "CUDA version check skipped"
    echo ""
    echo "environment.yml is configured for CUDA 11.8"
    echo "If you need CUDA 12.x, edit environment.yml line 29 before continuing"
    echo ""
    read -p "Press Enter to continue or Ctrl+C to abort and edit..."
else
    echo "WARNING: nvidia-smi not found. Installing with GPU support anyway."
    echo "If you don't have a GPU, edit environment.yml to remove pytorch-cuda"
    echo ""
    read -p "Press Enter to continue or Ctrl+C to abort..."
fi

# Create conda environment
echo ""
echo "Creating conda environment 'flair'..."
echo "This may take 5-10 minutes..."
echo ""

conda env create -f environment.yml

echo ""
echo "=========================================="
echo "Environment created successfully!"
echo "=========================================="
echo ""

# Activate environment (note: this only works within the script)
echo "Activating environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate flair

# Create __init__.py files
echo "Creating package structure..."
touch src/__init__.py
touch src/extract/__init__.py
touch src/cnn/__init__.py
touch src/infer/__init__.py
touch src/eval/__init__.py
touch tests/__init__.py

# Install package in editable mode
echo "Installing flAIr package in development mode..."
pip install -e .

echo ""
echo "=========================================="
echo "Verifying installation..."
echo "=========================================="
echo ""

# Test imports
python -c "
import numpy as np
import pandas as pd
import pysam
import torch

print('NumPy:', np.__version__)
print('Pandas:', pd.__version__)
print('PyTorch:', torch.__version__)
print('Pysam:', pysam.__version__)
print('CUDA available:', torch.cuda.is_available())
if torch.cuda.is_available():
    print('GPU count:', torch.cuda.device_count())
    print('GPU name:', torch.cuda.get_device_name(0))
print('')
print('All Python packages imported successfully!')
"

# Test command-line tools
echo ""
echo "Command-line tools:"
samtools --version | head -1
bedtools --version
tabix --version 2>&1 | head -1

echo ""
echo "=========================================="
echo "Setup Complete!"
echo "=========================================="
echo ""
echo "To use the environment:"
echo "  conda activate flair"
echo ""
echo "Next steps:"
echo "  1. Process your prior databases:"
echo "     python src/extract/standardize_priors.py \\"
echo "         --tss-fantom TSS_db/FANTOM_TSS_human.bed \\"
echo "         --tss-reftss TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt \\"
echo "         --tes TTS_db/polyadb.hg38.weighted.bed \\"
echo "         --junction splice_db/detail_Linear_splice_annotation.txt \\"
echo "         --output data/ref/"
echo ""
echo "  2. Download reference genome (if needed):"
echo "     See GETTING_STARTED.md for details"
echo ""
