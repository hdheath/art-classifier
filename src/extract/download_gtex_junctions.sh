#!/bin/bash
#
# Download GTEx v8 junction data
#
# GTEx provides comprehensive splice junction data across all tissues.
# This script downloads the junction files and prepares them for use with flAIr.
#
# Usage:
#   bash src/extract/download_gtex_junctions.sh [output_dir]
#
# Default output: data/downloads/gtex/

set -e  # Exit on error
set -u  # Exit on undefined variable

# Configuration
OUTPUT_DIR="${1:-data/downloads/gtex}"
GTEX_BASE="https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data"

# Create output directory
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

echo "================================"
echo "GTEx v8 Junction Download"
echo "================================"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Function to download with retry
download_file() {
    local url=$1
    local output=$2
    local desc=$3

    if [ -f "${output}" ]; then
        echo "✓ Already downloaded: ${output}"
        return 0
    fi

    echo "Downloading: ${desc}"
    echo "  URL: ${url}"
    echo "  Output: ${output}"

    # Try wget first, fall back to curl
    if command -v wget &> /dev/null; then
        wget --continue --progress=bar:force "${url}" -O "${output}.tmp"
    elif command -v curl &> /dev/null; then
        curl -L -C - "${url}" -o "${output}.tmp"
    else
        echo "ERROR: Neither wget nor curl found!"
        exit 1
    fi

    # Move to final location if successful
    mv "${output}.tmp" "${output}"
    echo "✓ Download complete: ${output}"
    echo ""
}

# Option 1: Full GTEx junction file (LARGE - ~3GB)
echo "[1/2] Full GTEx Junction File"
echo "WARNING: This file is ~3GB compressed"
echo ""

JUNCTION_FILE="GTEx_Analysis_2017-06-05_v8_STARjunctions.gct.gz"
JUNCTION_URL="${GTEX_BASE}/${JUNCTION_FILE}"

read -p "Download full junction file? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    download_file "${JUNCTION_URL}" "${JUNCTION_FILE}" "GTEx Full Junctions"
else
    echo "Skipping full junction file"
    echo ""
fi

# Option 2: Gene-level reads (alternative, smaller)
echo "[2/2] GTEx Gene Read Counts"
echo "This file is smaller and can be used to extract junction weights"
echo ""

GENE_FILE="GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
GENE_URL="${GTEX_BASE}/${GENE_FILE}"

read -p "Download gene read counts? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    download_file "${GENE_URL}" "${GENE_FILE}" "GTEx Gene Reads"
else
    echo "Skipping gene read counts"
    echo ""
fi

# Download sample attributes (useful for filtering by tissue)
echo "[3/3] GTEx Sample Attributes"
SAMPLE_FILE="GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
SAMPLE_URL="https://storage.googleapis.com/gtex_analysis_v8/annotations/${SAMPLE_FILE}"

download_file "${SAMPLE_URL}" "${SAMPLE_FILE}" "GTEx Sample Attributes"

echo "================================"
echo "Download Summary"
echo "================================"
echo "Files downloaded to: ${OUTPUT_DIR}"
ls -lh *.gz *.txt 2>/dev/null || echo "No files found"
echo ""

echo "================================"
echo "Next Steps"
echo "================================"
echo ""
echo "1. Parse GTEx junctions (if downloaded):"
echo "   python src/extract/parse_gtex_junctions.py \\"
echo "       --input ${OUTPUT_DIR}/${JUNCTION_FILE} \\"
echo "       --output data/ref/gtex_junctions.tsv \\"
echo "       --min-samples 5 \\"
echo "       --min-reads 10"
echo ""
echo "2. Filter by tissue (optional):"
echo "   python src/extract/parse_gtex_junctions.py \\"
echo "       --input ${OUTPUT_DIR}/${JUNCTION_FILE} \\"
echo "       --samples ${OUTPUT_DIR}/${SAMPLE_FILE} \\"
echo "       --tissues BRAIN HEART LIVER \\"
echo "       --output data/ref/gtex_junctions_tissues.tsv"
echo ""
echo "3. OR use your existing RJunBase data (recommended):"
echo "   python src/extract/standardize_priors.py \\"
echo "       --junction splice_db/detail_Linear_splice_annotation.txt \\"
echo "       --output data/ref/"
echo ""
echo "RECOMMENDATION:"
echo "Your RJunBase data (detail_Linear_splice_annotation.txt) already"
echo "contains ~800K high-quality junctions. You likely don't need"
echo "additional GTEx data unless you want tissue-specific subsets."
echo ""
