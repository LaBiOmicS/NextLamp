#!/bin/bash
# ==============================================================================
#  NextLAMP Systematized Pipeline — Full Babesia Dataset Run
# ==============================================================================
set -e

# Base directory setup
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DATA_DIR="${BASE_DIR}/data"
CONFIG_FILE="${BASE_DIR}/configs/babesia_full.yaml"
RESULTS_DIR="${BASE_DIR}/results/babesia_full"
LOGS_DIR="${BASE_DIR}/logs"

# Tools and Binaries
BOWTIE2_PATH="/home/fabiano.menegidio/miniforge3/envs/humann3_env/bin/bowtie2"
BOWTIE2_BUILD="/home/fabiano.menegidio/miniforge3/envs/humann3_env/bin/bowtie2-build"
PYTHON="/home/fabiano.menegidio/miniforge3/bin/python"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOGS_DIR}/nextlamp_babesia_full_${TIMESTAMP}.log"

mkdir -p "${RESULTS_DIR}"
mkdir -p "${LOGS_DIR}"

exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=========================================================================="
echo "  NextLAMP Whole-Genome LAMP Primer Design — Full Babesia Dataset"
echo "  Started at: $(date)"
echo "  Configuration file: ${CONFIG_FILE}"
echo "  Log file: ${LOG_FILE}"
echo "=========================================================================="

BT2_IDX="${DATA_DIR}/db_completo_bt2_idx"

# Step 1: Check / Build Bowtie 2 index for db_completo.fa
echo ""
echo "--- Step 1: Verifying Bowtie 2 Index for Complete Background Database ---"
if [ ! -f "${BT2_IDX}.1.bt2l" ] && [ ! -f "${BT2_IDX}.1.bt2" ]; then
    echo "Bowtie 2 index not found. Building large index for ${DATA_DIR}/db_completo.fa..."
    ${BOWTIE2_BUILD} --large-index "${DATA_DIR}/db_completo.fa" "${BT2_IDX}"
    echo "Index construction completed successfully."
else
    echo "Bowtie 2 index verified (${BT2_IDX})."
fi

# Step 2: Execute NextLAMP Pipeline with YAML configuration
echo ""
echo "--- Step 2: Running NextLAMP Pipeline ---"
export PYTHONPATH="${BASE_DIR}:${PYTHONPATH}"

${PYTHON} "${BASE_DIR}/scripts/nextlamp/run_nextlamp.py" \
    --config "${CONFIG_FILE}" \
    --bowtie2-path "${BOWTIE2_PATH}"

echo ""
echo "=========================================================================="
echo "  Pipeline Execution Finished Successfully!"
echo "  Finished at: $(date)"
echo "  Output Directory: ${RESULTS_DIR}"
echo "=========================================================================="
