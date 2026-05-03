#!/bin/bash
# =============================================================================
# Artifact Classifier — SLURM Training Job (TSS + TTS, with width regressor)
# Phoenix cluster | hdheath@ucsc.edu
#
# This job runs the full training pipeline end-to-end for both boundary types:
#   Step 1: prepare_data.py — extract sequence features from GRCh38 around
#           positives (hg38-verified) + random/hard negatives.
#   Step 2: train.py — XGBoost / Logistic / RF classifier + (TSS only)
#           XGBoost width regressor on log1p(width).
#   Step 3: evaluate.py — PR/ROC curves, feature importance, calibration.
#
# IMPORTANT: this is the first re-training run after the hg19/hg38 build
# audit (commit f9e41ac). The previous TSS classifier was trained on
# hg19 positives evaluated against the hg38 genome — invalid. Any
# data/derived/*_features.npz files on disk before this run are stale
# and will be deleted by this script.
#
# Submit:
#   sbatch submit_training.sbatch
#
# Monitor:
#   squeue --me
#   tail -f logs/flAIr_training_<jobid>.out
# =============================================================================

#SBATCH --job-name=artclf_train
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --output=/private/home/hdheath/projects/artifact_classifier/logs/artclf_training_%j.out
#SBATCH --error=/private/home/hdheath/projects/artifact_classifier/logs/artclf_training_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hdheath@ucsc.edu

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths — point at the artifact_classifier repo, NOT the old flAIr tree.
# Reference databases live on brookslab; the repo links them via TSS_db/
# and TTS_db/ symlinks.
# ---------------------------------------------------------------------------
PROJ=/private/home/hdheath/projects/artifact_classifier
GENOME=${PROJ}/data/ref/GRCh38.primary_assembly.genome.fa
GTF=${PROJ}/data/ref/gencode.v49.basic.annotation.gtf

# TSS sources (all verified hg38 in commit f9e41ac):
FANTOM_TSS=${PROJ}/TSS_db/FANTOM_TSS_human.hg38.bed                          # positives source
FANTOM_FAIR=${PROJ}/TSS_db/fantom5_hg38/hg38_fair_CAGE_peaks_phase1and2.bed  # width target source

# TTS sources (verified byte-identical to fresh re-lift):
POLYADB=${PROJ}/TTS_db/polyadb.hg38.weighted.bed
POLYADB_PAS=${PROJ}/TTS_db/polyAdb_human.PAS.txt

DATA_DIR=${PROJ}/data/derived
MODEL_DIR=${PROJ}/model_weights
EVAL_DIR=${DATA_DIR}/eval
SRC=${PROJ}/src/flair_boundary

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
MAX_SITES=200000              # Max positive sites per type
NEG_RATIO=3                   # Negatives per positive
HARD_NEG_FRAC=0.5             # Half of negatives drawn from inside gene bodies
CV_FOLDS=5
HOLDOUT_CHROM=chr22
THREADS=${SLURM_CPUS_PER_TASK}
MODELS="xgboost logistic rf"

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
echo "============================================================"
echo "Artifact Classifier — Training Run"
echo "============================================================"
echo "Job ID:      ${SLURM_JOB_ID}"
echo "Node:        $(hostname)"
echo "CPUs:        ${SLURM_CPUS_PER_TASK}"
echo "Memory:      ${SLURM_MEM_PER_NODE} MB"
echo "Start:       $(date)"
echo "Repo:        ${PROJ}"
echo "Genome:      ${GENOME}"
echo "GTF:         ${GTF}"
echo "TSS pos:     ${FANTOM_TSS}"
echo "TSS width:   ${FANTOM_FAIR}"
echo "TTS:         ${POLYADB}"
echo "Max sites:   ${MAX_SITES} pos per type"
echo "Neg ratio:   ${NEG_RATIO}x  (hard fraction: ${HARD_NEG_FRAC})"
echo "Holdout:     ${HOLDOUT_CHROM}"
echo "Models:      ${MODELS}"
echo "============================================================"
echo

# Use the flair-ml conda env's interpreter directly. SLURM nodes have a
# system Python 3.10 (pandas 1.5) on PATH; `conda activate` inside a
# non-interactive job doesn't reliably shadow it. Calling the env's
# python by absolute path guarantees we get pandas 2.x / py3.11.
source /private/home/hdheath/miniforge3/etc/profile.d/conda.sh
conda activate flair-ml
PY=/private/home/hdheath/miniforge3/envs/flair-ml/bin/python

# Sanity-check we got the right interpreter, not the system one.
which_py=$($PY -c "import sys; print(sys.executable)")
py_ver=$($PY -c "import sys; print(sys.version_info[:2])")
pd_ver=$($PY -c "import pandas; print(pandas.__version__)")
echo "Interpreter: $which_py"
echo "Python:      $py_ver"
echo "Pandas:      $pd_ver"
if [[ "$which_py" != *"flair-ml"* ]]; then
    echo "ERROR: not running inside flair-ml env" >&2
    exit 1
fi
$PY -c "import xgboost, pysam, sklearn, numpy, pandas; print('Dependencies OK')"

mkdir -p "${DATA_DIR}" "${MODEL_DIR}" "${EVAL_DIR}" "${PROJ}/logs"

# ---------------------------------------------------------------------------
# Force regeneration: any cached features.npz from before the hg38 build
# fix is stale. Delete to be safe.
# ---------------------------------------------------------------------------
echo "Deleting any cached *_features.npz to force regeneration with hg38-clean data..."
for f in "${DATA_DIR}/tss_features.npz" "${DATA_DIR}/tts_features.npz"; do
    if [ -f "$f" ] || [ -L "$f" ]; then
        echo "  rm $f"
        rm -f "$f"
    fi
done

# ---------------------------------------------------------------------------
# Helper: run one boundary type
# ---------------------------------------------------------------------------
run_boundary_type() {
    local BTYPE="$1"
    local T_START=$(date +%s)

    echo
    echo "============================================================"
    echo "  ${BTYPE^^} — Step 1/3: Feature extraction"
    echo "  $(date)"
    echo "============================================================"

    local FEATURES=${DATA_DIR}/${BTYPE}_features.npz

    local PREPARE_ARGS=(
        --type "${BTYPE}"
        --genome "${GENOME}"
        --output "${FEATURES}"
        --max-sites "${MAX_SITES}"
        --neg-ratio "${NEG_RATIO}"
        --hard-neg-frac "${HARD_NEG_FRAC}"
        --merge-window 25
        --neg-min-dist 500
        --threads "${THREADS}"
    )

    if [ "${BTYPE}" == "tss" ]; then
        PREPARE_ARGS+=(--fantom "${FANTOM_TSS}")
        [ -f "${FANTOM_FAIR}" ] && PREPARE_ARGS+=(--fantom-fair "${FANTOM_FAIR}")
    else
        PREPARE_ARGS+=(--polyadb "${POLYADB}")
        [ -f "${POLYADB_PAS}" ] && PREPARE_ARGS+=(--polyadb-pas "${POLYADB_PAS}")
    fi

    # GTF needed for hard negatives (gene-body sampling) and proximity context.
    [ -f "${GTF}" ] && PREPARE_ARGS+=(--gtf "${GTF}")

    PYTHONPATH=${PROJ}/src $PY -m flair_boundary.extract.prepare_data "${PREPARE_ARGS[@]}"

    echo
    echo "============================================================"
    echo "  ${BTYPE^^} — Step 2/3: Training models"
    echo "  $(date)"
    echo "============================================================"

    # --width-model trains the regressor only for TSS (the train.py code
    # gates this internally — passing it for TTS is a no-op with a log
    # line saying so).
    PYTHONPATH=${PROJ}/src $PY -m flair_boundary.train.train \
        --features "${FEATURES}" \
        --output "${MODEL_DIR}/${BTYPE}" \
        --type "${BTYPE}" \
        --model ${MODELS} \
        --cv-folds "${CV_FOLDS}" \
        --holdout-chrom "${HOLDOUT_CHROM}" \
        --feature-importance \
        --early-stopping 50 \
        --width-model

    echo
    echo "============================================================"
    echo "  ${BTYPE^^} — Step 3/3: Evaluation"
    echo "  $(date)"
    echo "============================================================"

    PYTHONPATH=${PROJ}/src $PY -m flair_boundary.eval.evaluate \
        --features "${FEATURES}" \
        --model-dir "${MODEL_DIR}/${BTYPE}" \
        --output "${EVAL_DIR}/${BTYPE}" \
        --type "${BTYPE}" \
        --holdout-chrom "${HOLDOUT_CHROM}"

    local T_END=$(date +%s)
    local ELAPSED=$(( (T_END - T_START) / 60 ))
    echo
    echo "  ${BTYPE^^} complete in ${ELAPSED} minutes."
    echo "  Models:   ${MODEL_DIR}/${BTYPE}/"
    echo "  Eval:     ${EVAL_DIR}/${BTYPE}/"
}

# ---------------------------------------------------------------------------
# Run TSS, then TTS
# ---------------------------------------------------------------------------
JOB_START=$(date +%s)

run_boundary_type "tss"
run_boundary_type "tts"

JOB_END=$(date +%s)
TOTAL=$(( (JOB_END - JOB_START) / 60 ))

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo
echo "============================================================"
echo "ALL DONE — Total time: ${TOTAL} minutes"
echo "$(date)"
echo "============================================================"
echo
echo "Feature matrices:"
ls -lh ${DATA_DIR}/*.npz 2>/dev/null || echo "  (none)"
echo
echo "Trained models:"
ls -lh ${MODEL_DIR}/tss/*.pkl ${MODEL_DIR}/tts/*.pkl 2>/dev/null || echo "  (none)"
echo
echo "Evaluation summaries:"
for f in ${EVAL_DIR}/tss/tss_evaluation_summary.json ${EVAL_DIR}/tts/tts_evaluation_summary.json; do
    if [ -f "$f" ]; then
        echo "  $(basename $f):"
        $PY -c "
import json
with open('$f') as fp:
    d = json.load(fp)
for model, metrics in d.items():
    print(f'    {model:<24}  AUROC={metrics.get(\"roc_auc\", 0):.4f}  AUPRC={metrics.get(\"avg_precision\", 0):.4f}  F1={metrics.get(\"f1_optimal\", 0):.4f}')
"
    fi
done
echo
echo "Width regressor (TSS):"
if [ -f "${MODEL_DIR}/tss/training_metadata.json" ]; then
    $PY -c "
import json
with open('${MODEL_DIR}/tss/training_metadata.json') as fp:
    d = json.load(fp)
w = d.get('width_model')
if w:
    h = w.get('holdout_metrics', {})
    print(f\"  RMSE(log1p)={h.get('rmse_log1p', 0):.4f}  MAE(bp)={h.get('mae_bp', 0):.1f}  median|err|(bp)={h.get('median_ae_bp', 0):.1f}  Spearman r={h.get('spearman_r', 0):.4f}\")
else:
    print('  No width regressor metadata')
"
fi
