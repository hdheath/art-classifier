#!/bin/bash
# TTS-only training run. TSS already finished and validated; this job
# only does the TTS prepare/train/eval steps so we don't redo TSS work.

#SBATCH --job-name=artclf_tts
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/private/home/hdheath/projects/artifact_classifier/logs/artclf_tts_%j.out
#SBATCH --error=/private/home/hdheath/projects/artifact_classifier/logs/artclf_tts_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hdheath@ucsc.edu

set -euo pipefail

PROJ=/private/home/hdheath/projects/artifact_classifier
GENOME=${PROJ}/data/ref/GRCh38.primary_assembly.genome.fa
GTF=${PROJ}/data/ref/gencode.v49.basic.annotation.gtf
POLYADB=${PROJ}/TTS_db/polyadb.hg38.weighted.bed
POLYADB_PAS=${PROJ}/TTS_db/polyAdb_human.PAS.txt
DATA_DIR=${PROJ}/data/derived
MODEL_DIR=${PROJ}/model_weights
EVAL_DIR=${DATA_DIR}/eval
SRC=${PROJ}/src/flair_boundary
THREADS=${SLURM_CPUS_PER_TASK}

source /private/home/hdheath/miniforge3/etc/profile.d/conda.sh
conda activate flair-ml
PY=/private/home/hdheath/miniforge3/envs/flair-ml/bin/python
which_py=$($PY -c "import sys; print(sys.executable)")
if [[ "$which_py" != *"flair-ml"* ]]; then
    echo "ERROR: not running in flair-ml env (got $which_py)" >&2
    exit 1
fi
echo "Interpreter: $which_py"
$PY -c "import xgboost, pysam, sklearn, numpy, pandas; print('Dependencies OK')"

mkdir -p "${DATA_DIR}" "${MODEL_DIR}/tts" "${EVAL_DIR}/tts"

# Force regeneration of TTS features (in case stale).
rm -f ${DATA_DIR}/tts_features.npz

echo
echo "============================================================"
echo "TTS — Step 1/3: Feature extraction"
echo "$(date)"
echo "============================================================"
PYTHONPATH=${PROJ}/src $PY -m flair_boundary.extract.prepare_data \
    --type tts \
    --genome "${GENOME}" \
    --gtf "${GTF}" \
    --polyadb "${POLYADB}" \
    --polyadb-pas "${POLYADB_PAS}" \
    --output ${DATA_DIR}/tts_features.npz \
    --max-sites 200000 --neg-ratio 3 --hard-neg-frac 0.5 \
    --merge-window 25 --neg-min-dist 500 \
    --threads ${THREADS}

echo
echo "============================================================"
echo "TTS — Step 2/3: Training"
echo "$(date)"
echo "============================================================"
PYTHONPATH=${PROJ}/src $PY -m flair_boundary.train.train \
    --features ${DATA_DIR}/tts_features.npz \
    --output ${MODEL_DIR}/tts \
    --type tts \
    --model xgboost logistic rf \
    --cv-folds 5 --holdout-chrom chr22 \
    --feature-importance --early-stopping 50

echo
echo "============================================================"
echo "TTS — Step 3/3: Evaluation"
echo "$(date)"
echo "============================================================"
PYTHONPATH=${PROJ}/src $PY -m flair_boundary.eval.evaluate \
    --features ${DATA_DIR}/tts_features.npz \
    --model-dir ${MODEL_DIR}/tts \
    --output ${EVAL_DIR}/tts \
    --type tts \
    --holdout-chrom chr22

echo
echo "Done at $(date)"
