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
if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
    conda activate ncbi_env 2>/dev/null || conda activate nextlamp 2>/dev/null || true
fi

BOWTIE2_PATH="$(which bowtie2 2>/dev/null || echo '/home/fabiano.menegidio/miniforge3/envs/ncbi_env/bin/bowtie2')"
BOWTIE2_BUILD="$(which bowtie2-build 2>/dev/null || echo '/home/fabiano.menegidio/miniforge3/envs/ncbi_env/bin/bowtie2-build')"
PYTHON="$(which python3 || echo '/home/fabiano.menegidio/miniforge3/bin/python')"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOGS_DIR}/nextlamp_babesia_full_${TIMESTAMP}.log"

mkdir -p "${RESULTS_DIR}"
mkdir -p "${LOGS_DIR}"
mkdir -p "${BASE_DIR}/results/babesia_paper_export/tables"
mkdir -p "${BASE_DIR}/results/babesia_paper_export/reports"
mkdir -p "${BASE_DIR}/results/babesia_paper_export/raw_data"

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

# Step 3: Organise Publication Package
echo ""
echo "--- Step 3: Structuring Publication Package ---"
PAPER_DIR="${BASE_DIR}/results/babesia_paper_export"

# Copy TSV tables and TXT reports
cp "${PAPER_DIR}/raw_data/nextlamp_babesia_results_primers.tsv" "${PAPER_DIR}/tables/Table_S1_Babesia_LAMP_Primers.tsv" 2>/dev/null || true
cp "${PAPER_DIR}/raw_data/nextlamp_babesia_results_report.txt" "${PAPER_DIR}/reports/NextLAMP_Execution_Summary.txt" 2>/dev/null || true

# Generate README for paper collaborators
cat << 'EOF' > "${PAPER_DIR}/README_PAPER_EXPORT.md"
# NextLAMP Results Package — Babesia Genome-Wide Primer Design

This directory contains the structured dataset and publication-ready tables for the scientific article.

## Directory Structure
- `tables/`: Contains formatted TSV/CSV tables for manuscript Supplementary Material.
  - `Table_S1_Babesia_LAMP_Primers.tsv`: Primer sequences, Tm, GC%, alignment specificity scores, and amplicon locations.
- `reports/`: Execution logs and summary reports for reproducibility.
- `raw_data/`: Complete JSON output bundle containing full metadata, parameter sets, and raw alignments.

## Reproducibility Metadata
- Target Sequence: `data/target_babesia_canis.fa`
- Background DB: `data/db_completo.fa` (501 Apicomplexa, 4 Dog, 1 Cat, 21 Tick genomes)
- Tool Version: NextLAMP v1.0.0
EOF

echo ""
echo "=========================================================================="
echo "  Pipeline Execution Finished Successfully!"
echo "  Finished at: $(date)"
echo "  Publication Package Directory: ${PAPER_DIR}"
echo "=========================================================================="
