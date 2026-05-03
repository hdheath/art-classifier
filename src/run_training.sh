#!/usr/bin/env bash
# =============================================================================
# flAIr: Hybrid Boundary Classifier Training Pipeline
# =============================================================================
#
# Trains XGBoost (and optionally other) classifiers for TSS and TTS boundary
# detection using reference databases and the GRCh38 genome.
#
# Usage:
#   conda activate flair-ml
#   bash src/run_training.sh [--type tss|tts|both] [options]
#
# Quick start (both models):
#   bash src/run_training.sh
#
# TSS only:
#   bash src/run_training.sh --type tss
#
# Full run with all model variants:
#   bash src/run_training.sh --models "xgboost logistic rf xgboost_calibrated"
#
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults (override via command-line flags)
# ---------------------------------------------------------------------------
BOUNDARY_TYPE="both"        # tss | tts | both
MODELS="xgboost logistic rf"
MAX_SITES=200000            # Max positive sites per type (0 = all)
NEG_RATIO=3                 # Negatives per positive
CV_FOLDS=5
HOLDOUT_CHROM="chr22"
THREADS=4
MERGE_WINDOW=25

# Project paths
PROJ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GENOME="${PROJ_DIR}/data/ref/GRCh38.primary_assembly.genome.fa"
GTF="${PROJ_DIR}/data/ref/gencode.v49.basic.annotation.gtf"
# FANTOM_TSS: the verified-hg38 (re-lifted) permissive peaks. NEVER use the
# original FANTOM_TSS_human.bed — that file is hg19 (build audit Apr 2026).
FANTOM_TSS="${PROJ_DIR}/TSS_db/FANTOM_TSS_human.hg38.bed"
# FANTOM_FAIR: hg38-native robust merged CAGE clusters (Lizio et al. 2017).
# Used as the width-target source via spatial overlap with positives.
FANTOM_FAIR="${PROJ_DIR}/TSS_db/fantom5_hg38/hg38_fair_CAGE_peaks_phase1and2.bed"
POLYADB="${PROJ_DIR}/TTS_db/polyadb.hg38.weighted.bed"
POLYADB_PAS="${PROJ_DIR}/TTS_db/polyAdb_human.PAS.txt"

DATA_DIR="${PROJ_DIR}/data/derived"
MODEL_DIR="${PROJ_DIR}/model_weights"
EVAL_DIR="${PROJ_DIR}/data/derived/eval"

SRC="${PROJ_DIR}/src/flair_boundary"

# ---------------------------------------------------------------------------
# Parse CLI arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --type)       BOUNDARY_TYPE="$2"; shift 2 ;;
        --models)     MODELS="$2"; shift 2 ;;
        --max-sites)  MAX_SITES="$2"; shift 2 ;;
        --neg-ratio)  NEG_RATIO="$2"; shift 2 ;;
        --cv-folds)   CV_FOLDS="$2"; shift 2 ;;
        --holdout)    HOLDOUT_CHROM="$2"; shift 2 ;;
        --threads)    THREADS="$2"; shift 2 ;;
        --genome)     GENOME="$2"; shift 2 ;;
        *)            echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate environment
# ---------------------------------------------------------------------------
echo "============================================================"
echo "flAIr Hybrid Boundary Classifier Training"
echo "============================================================"
echo ""
echo "Project dir: ${PROJ_DIR}"
echo "Genome:      ${GENOME}"
echo "Type(s):     ${BOUNDARY_TYPE}"
echo "Models:      ${MODELS}"
echo "Max sites:   ${MAX_SITES}"
echo "Neg ratio:   ${NEG_RATIO}x"
echo "CV folds:    ${CV_FOLDS}"
echo "Holdout:     ${HOLDOUT_CHROM}"
echo ""

if [ ! -f "${GENOME}" ]; then
    echo "ERROR: Genome not found: ${GENOME}"
    exit 1
fi

if ! python -c "import xgboost, pysam, sklearn" 2>/dev/null; then
    echo "ERROR: Required packages missing. Run: conda activate flair-ml"
    exit 1
fi

mkdir -p "${DATA_DIR}" "${MODEL_DIR}" "${EVAL_DIR}"

# ---------------------------------------------------------------------------
# Helper: run one boundary type
# ---------------------------------------------------------------------------

run_type() {
    local btype="$1"
    echo ""
    echo "============================================================"
    echo "  Processing: ${btype^^}"
    echo "============================================================"

    local features_file="${DATA_DIR}/${btype}_features.npz"
    local model_out="${MODEL_DIR}/${btype}"
    local eval_out="${EVAL_DIR}/${btype}"

    # ---- Step 1: Prepare training data ----
    if [ -f "${features_file}" ]; then
        echo ""
        echo "[Step 1] Features already exist: ${features_file}"
        echo "  Delete to regenerate: rm ${features_file}"
    else
        echo ""
        echo "[Step 1] Preparing ${btype^^} training data..."

        local prepare_args=(
            --type "${btype}"
            --genome "${GENOME}"
            --output "${features_file}"
            --neg-ratio "${NEG_RATIO}"
            --merge-window "${MERGE_WINDOW}"
            --threads "${THREADS}"
        )

        if [ "${MAX_SITES}" -gt 0 ]; then
            prepare_args+=(--max-sites "${MAX_SITES}")
        fi

        if [ "${btype}" == "tss" ]; then
            prepare_args+=(--fantom "${FANTOM_TSS}")
            # FANTOM5-hg38 fair peaks supply the width regression target via
            # spatial overlap with positives. Optional — TSS classifier still
            # trains without it, just no width head.
            if [ -f "${FANTOM_FAIR}" ]; then
                prepare_args+=(--fantom-fair "${FANTOM_FAIR}")
            fi
        else
            prepare_args+=(--polyadb "${POLYADB}")
            if [ -f "${POLYADB_PAS}" ]; then
                prepare_args+=(--polyadb-pas "${POLYADB_PAS}")
            fi
        fi
        # GTF is needed by both types: TSS uses it for proximity context,
        # TTS uses it for hard negatives sampled from inside gene bodies.
        if [ -f "${GTF}" ]; then
            prepare_args+=(--gtf "${GTF}")
        fi

        python "${SRC}/extract/prepare_data.py" "${prepare_args[@]}"
    fi

    # ---- Step 2: Train models ----
    echo ""
    echo "[Step 2] Training ${btype^^} classifiers..."

    python "${SRC}/train/train.py" \
        --features "${features_file}" \
        --output "${model_out}" \
        --type "${btype}" \
        --model ${MODELS} \
        --cv-folds "${CV_FOLDS}" \
        --holdout-chrom "${HOLDOUT_CHROM}" \
        --feature-importance \
        --early-stopping 50 \
        --width-model

    # ---- Step 3: Evaluate ----
    echo ""
    echo "[Step 3] Evaluating ${btype^^} classifiers..."

    python "${SRC}/eval/evaluate.py" \
        --features "${features_file}" \
        --model-dir "${model_out}" \
        --output "${eval_out}" \
        --type "${btype}" \
        --holdout-chrom "${HOLDOUT_CHROM}"

    echo ""
    echo "  ${btype^^} complete!"
    echo "  Features:  ${features_file}"
    echo "  Models:    ${model_out}/"
    echo "  Eval:      ${eval_out}/"
}

# ---------------------------------------------------------------------------
# Run requested type(s)
# ---------------------------------------------------------------------------

if [ "${BOUNDARY_TYPE}" == "tss" ] || [ "${BOUNDARY_TYPE}" == "both" ]; then
    run_type "tss"
fi

if [ "${BOUNDARY_TYPE}" == "tts" ] || [ "${BOUNDARY_TYPE}" == "both" ]; then
    run_type "tts"
fi

echo ""
echo "============================================================"
echo "Training complete!"
echo "============================================================"
echo ""
echo "Model weights:"
ls -lh "${MODEL_DIR}/"*/*.pkl 2>/dev/null || echo "  (none yet)"
echo ""
echo "Evaluation reports:"
ls "${EVAL_DIR}/"*/*.json 2>/dev/null || echo "  (none yet)"
echo ""
echo "Next steps:"
echo "  - Review eval/ plots for model performance"
echo "  - Use xgboost.pkl (best AUPRC) for inference"
echo "  - Run inference:"
echo "      python src/flair_boundary/infer/infer.py \\"
echo "          --candidates candidates.bed \\"
echo "          --genome ${GENOME} \\"
echo "          --tss-model ${MODEL_DIR}/tss/xgboost.pkl \\"
echo "          --tts-model ${MODEL_DIR}/tts/xgboost.pkl \\"
echo "          --bam aligned_reads.bam \\"
echo "          --type tss \\"
echo "          --output scored_boundaries.tsv \\"
echo "          --select"
