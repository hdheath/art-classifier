#!/bin/bash
#
# Download pre-computed splice junctions from recount3
#
# recount3 provides READY-TO-USE junction files - no inference needed!
#
# Available data:
# - GTEx: All tissues, junctions already computed from STAR alignments
# - TCGA: Cancer samples, junctions pre-computed
# - SRA: Public data, junctions pre-computed
#
# Usage:
#   bash src/extract/download_recount3_junctions.sh gtex
#   bash src/extract/download_recount3_junctions.sh tcga
#

set -e

SOURCE="${1:-gtex}"
OUTPUT_DIR="data/downloads/recount3/${SOURCE}"

mkdir -p "${OUTPUT_DIR}"

echo "========================================"
echo "recount3 Pre-computed Junction Download"
echo "========================================"
echo "Source: ${SOURCE}"
echo "Output: ${OUTPUT_DIR}"
echo ""

# recount3 provides junctions via their R package or direct URLs
# The junctions are ALREADY COMPUTED from STAR alignments

case "${SOURCE}" in
    gtex)
        echo "GTEx junctions are available through recount3"
        echo ""
        echo "METHOD 1: Direct R package (RECOMMENDED)"
        echo "========================================="
        echo ""
        cat << 'EOF'
# Install recount3 in R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("recount3")

library(recount3)

# Get GTEx projects (by tissue)
human_projects <- available_projects()
gtex_projects <- subset(human_projects, file_source == "gtex")

# Download junction data for specific project
# Option 1: Get all GTEx data
rse_gtex <- create_rse_manual(
    project = "gtex",
    project_home = "data_sources/gtex",
    type = "jxn",  # JUNCTION data (pre-computed!)
    organism = "human",
    annotation = "gencode_v26"
)

# Extract junction coordinates
jxn_coords <- rowRanges(rse_gtex)
jxn_counts <- assays(rse_gtex)$counts

# Export to BED format
jxn_df <- data.frame(
    chr = as.character(seqnames(jxn_coords)),
    start = start(jxn_coords),
    end = end(jxn_coords),
    strand = as.character(strand(jxn_coords)),
    counts = rowSums(jxn_counts)
)

write.table(jxn_df, "gtex_junctions.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
EOF
        echo ""
        echo "METHOD 2: Direct file download"
        echo "==============================="
        echo ""
        echo "recount3 stores data on AWS S3:"
        echo "  s3://recount-opendata/recount3/..."
        echo ""
        echo "Junction files are in the 'jxn' subdirectory"
        echo ""
        echo "Example download using AWS CLI:"
        echo "  aws s3 ls s3://recount-opendata/recount3/human/data_sources/gtex/"
        echo "  aws s3 cp s3://recount-opendata/recount3/human/data_sources/gtex/jxn/ . --recursive"
        echo ""
        ;;

    tcga)
        echo "TCGA junctions available through recount3"
        echo ""
        cat << 'EOF'
library(recount3)

# Get TCGA projects
human_projects <- available_projects()
tcga_projects <- subset(human_projects, file_source == "tcga")

# Download junction data
rse_tcga <- create_rse_manual(
    project = "BRCA",  # or any TCGA project
    project_home = "data_sources/tcga",
    type = "jxn",  # Pre-computed junctions
    organism = "human"
)
EOF
        ;;

    *)
        echo "Unknown source: ${SOURCE}"
        echo "Available: gtex, tcga"
        exit 1
        ;;
esac

echo ""
echo "========================================"
echo "IMPORTANT DISCOVERY!"
echo "========================================"
echo ""
echo "recount3 junctions are ALREADY COMPUTED!"
echo ""
echo "They used STAR aligner to compute junctions from:"
echo "  - GTEx: ~17,000 samples across all tissues"
echo "  - TCGA: ~11,000 cancer samples"
echo "  - SRA: Public RNA-seq data"
echo ""
echo "Junction format in recount3:"
echo "  - GRanges object with junction coordinates"
echo "  - Read counts per junction per sample"
echo "  - Annotated vs novel junctions"
echo "  - Already on hg38 (GENCODE v26+)"
echo ""
echo "========================================"
echo "HOWEVER..."
echo "========================================"
echo ""
echo "You ALREADY HAVE excellent junction data:"
echo "  splice_db/detail_Linear_splice_annotation.txt"
echo ""
echo "This RJunBase file contains:"
echo "  - ~800,000 junctions"
echo "  - Multi-tissue support"
echo "  - Expression weights"
echo "  - Already processed and filtered"
echo ""
echo "RECOMMENDATION:"
echo "  Skip recount3 download and use your RJunBase data!"
echo ""
echo "To process your existing data:"
echo "  python src/extract/standardize_priors.py \\"
echo "      --junction splice_db/detail_Linear_splice_annotation.txt \\"
echo "      --output data/ref/"
echo ""
