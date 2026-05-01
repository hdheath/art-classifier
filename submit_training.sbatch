#!/bin/bash
# =============================================================================
# flAIr Hybrid Boundary Classifier — SLURM Training Job
# Phoenix cluster | hdheath@ucsc.edu
#
# Runs TSS and TTS training sequentially in one job (two steps each):
#   Step 1: prepare_data.py — feature extraction from GRCh38 genome
#   Step 2: train.py — XGBoost / LogisticReg / RandomForest with CV
#   Step 3: evaluate.py — PR/ROC curves, feature importance plots
#
# Resources:
#   - Feature extraction: CPU-bound (pysam genome fetches), benefits from
#     multi-threading (--threads 16). ~200K sites × 2 types ≈ 30-60 min.
#   - Training: single-node, XGBoost uses all cores via n_jobs=-1.
#   - Total estimate: 3-5 hours on 16 cores.
#
# Submit:
#   sbatch submit_training.sbatch
#
# Monitor:
#   squeue --me
#   tail -f logs/flAIr_training_<jobid>.out
# =============================================================================

#SBATCH --job-name=flAIr_train
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/private/groups/brookslab/hdheath/projects/flAIr/logs/flAIr_training_%j.out
#SBATCH --error=/private/groups/brookslab/hdheath/projects/flAIr/logs/flAIr_training_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hdheath@ucsc.edu

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJ=/private/groups/brookslab/hdheath/projects/flAIr
GENOME=${PROJ}/data/ref/GRCh38.primary_assembly.genome.fa
GTF=${PROJ}/data/ref/gencode.v48.annotation.gtf

FANTOM_TSS=${PROJ}/TSS_db/FANTOM_TSS_human.bed
REFTSS=${PROJ}/TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt
POLYADB=${PROJ}/TTS_db/polyadb.hg38.weighted.bed
POLYADB_PAS=${PROJ}/TTS_db/polyAdb_human.PAS.txt

DATA_DIR=${PROJ}/data/derived
MODEL_DIR=${PROJ}/model_weights
EVAL_DIR=${PROJ}/data/derived/eval
SRC=${PROJ}/src/flair_boundary

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
MAX_SITES=200000     # Positive sites per type (200K TSS + 200K TTS)
NEG_RATIO=3          # 3x negatives per positive → 800K total per type
CV_FOLDS=5
HOLDOUT_CHROM=chr22
THREADS=${SLURM_CPUS_PER_TASK}   # Use all allocated CPUs for feature extraction
MODELS="xgboost logistic rf"

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
echo "============================================================"
echo "flAIr Hybrid Boundary Classifier Training"
echo "============================================================"
echo "Job ID:     ${SLURM_JOB_ID}"
echo "Node:       $(hostname)"
echo "CPUs:       ${SLURM_CPUS_PER_TASK}"
echo "Memory:     ${SLURM_MEM_PER_NODE}MB"
echo "Start:      $(date)"
echo "Genome:     ${GENOME}"
echo "Max sites:  ${MAX_SITES} pos per type"
echo "Neg ratio:  ${NEG_RATIO}x"
echo "Models:     ${MODELS}"
echo "Holdout:    ${HOLDOUT_CHROM}"
echo "============================================================"
echo ""

# Activate conda environment
source /private/home/hdheath/miniforge3/etc/profile.d/conda.sh
conda activate flair-ml

# Verify dependencies
python -c "import xgboost, pysam, sklearn, numpy, pandas; print('Dependencies OK')"

mkdir -p "${DATA_DIR}" "${MODEL_DIR}" "${EVAL_DIR}"

# ---------------------------------------------------------------------------
# Helper function: run one boundary type
# ---------------------------------------------------------------------------
run_boundary_type() {
    local BTYPE="$1"
    local T_START=$(date +%s)

    echo ""
    echo "============================================================"
    echo "  ${BTYPE^^} — Step 1/3: Feature extraction"
    echo "  $(date)"
    echo "============================================================"

    local FEATURES=${DATA_DIR}/${BTYPE}_features.npz

    if [ -f "${FEATURES}" ]; then
        echo "  Features already exist: ${FEATURES}"
        echo "  Delete to regenerate."
    else
        local PREPARE_ARGS=(
            --type "${BTYPE}"
            --genome "${GENOME}"
            --output "${FEATURES}"
            --max-sites "${MAX_SITES}"
            --neg-ratio "${NEG_RATIO}"
            --merge-window 25
            --neg-min-dist 500
            --threads "${THREADS}"
        )

        if [ "${BTYPE}" == "tss" ]; then
            PREPARE_ARGS+=(--fantom "${FANTOM_TSS}" --reftss "${REFTSS}")
        else
            PREPARE_ARGS+=(--polyadb "${POLYADB}" --polyadb-pas "${POLYADB_PAS}")
        fi

        # GTF is used for gene-distance context features (both TSS and TTS)
        [ -f "${GTF}" ] && PREPARE_ARGS+=(--gtf "${GTF}")

        PYTHONPATH=${PROJ}/src python ${SRC}/extract/prepare_data.py "${PREPARE_ARGS[@]}"
    fi

    echo ""
    echo "============================================================"
    echo "  ${BTYPE^^} — Step 2/3: Training models"
    echo "  $(date)"
    echo "============================================================"

    PYTHONPATH=${PROJ}/src python ${SRC}/train/train.py \
        --features "${FEATURES}" \
        --output "${MODEL_DIR}/${BTYPE}" \
        --type "${BTYPE}" \
        --model ${MODELS} \
        --cv-folds "${CV_FOLDS}" \
        --holdout-chrom "${HOLDOUT_CHROM}" \
        --feature-importance \
        --early-stopping 50

    echo ""
    echo "============================================================"
    echo "  ${BTYPE^^} — Step 3/3: Evaluation"
    echo "  $(date)"
    echo "============================================================"

    PYTHONPATH=${PROJ}/src python ${SRC}/eval/evaluate.py \
        --features "${FEATURES}" \
        --model-dir "${MODEL_DIR}/${BTYPE}" \
        --output "${EVAL_DIR}/${BTYPE}" \
        --type "${BTYPE}" \
        --holdout-chrom "${HOLDOUT_CHROM}"

    local T_END=$(date +%s)
    local ELAPSED=$(( (T_END - T_START) / 60 ))
    echo ""
    echo "  ${BTYPE^^} complete in ${ELAPSED} minutes."
    echo "  Models:   ${MODEL_DIR}/${BTYPE}/"
    echo "  Eval:     ${EVAL_DIR}/${BTYPE}/"
}

# ---------------------------------------------------------------------------
# Run TSS then TTS
# ---------------------------------------------------------------------------
JOB_START=$(date +%s)

run_boundary_type "tss"
run_boundary_type "tts"

JOB_END=$(date +%s)
TOTAL=$(( (JOB_END - JOB_START) / 60 ))

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "ALL DONE — Total time: ${TOTAL} minutes"
echo "$(date)"
echo "============================================================"
echo ""
echo "Output summary:"
echo ""
echo "Feature matrices:"
ls -lh ${DATA_DIR}/*.npz 2>/dev/null || echo "  (none)"
echo ""
echo "Trained models:"
ls -lh ${MODEL_DIR}/tss/*.pkl ${MODEL_DIR}/tts/*.pkl 2>/dev/null || echo "  (none)"
echo ""
echo "Evaluation plots:"
ls ${EVAL_DIR}/tss/*.png ${EVAL_DIR}/tts/*.png 2>/dev/null | wc -l | xargs echo "  PNG files:"
echo ""
echo "Performance summaries:"
for f in ${EVAL_DIR}/tss/tss_evaluation_summary.json ${EVAL_DIR}/tts/tts_evaluation_summary.json; do
    if [ -f "$f" ]; then
        echo "  $(basename $f):"
        python -c "
import json
with open('$f') as fp:
    d = json.load(fp)
for model, metrics in d.items():
    print(f'    {model}: AUROC={metrics.get(\"roc_auc\", 0):.4f}  AUPRC={metrics.get(\"avg_precision\", 0):.4f}  F1={metrics.get(\"f1_optimal\", 0):.4f}')
"
    fi
done
echo ""
echo "To run inference:"
echo "  conda activate flair-ml"
echo "  PYTHONPATH=${PROJ}/src python ${SRC}/infer/infer.py \\"
echo "      --candidates candidates.bed \\"
echo "      --genome ${GENOME} \\"
echo "      --tss-model ${MODEL_DIR}/tss/xgboost.pkl \\"
echo "      --tts-model ${MODEL_DIR}/tts/xgboost.pkl \\"
echo "      --bam aligned_reads.bam \\"
echo "      --type tss \\"
echo "      --output scored_boundaries.tsv \\"
echo "      --select"
