#!/bin/bash
#SBATCH --job-name=nextlamp_run
#SBATCH --partition=fast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=nextlamp_run_%j.log
#SBATCH --error=nextlamp_run_%j.err

echo "Starting NextLAMP SLURM job..."
date

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${BASE_DIR}/data"

# Auto-detect Python and Bowtie 2
PYTHON_BIN="$(which python3 || which python)"
BOWTIE2_BUILD="$(which bowtie2-build || echo "/home/fabiano.menegidio/miniforge3/envs/humann3_env/bin/bowtie2-build")"

# 1. Build Bowtie 2 index if not present
if [ ! -f "${DATA_DIR}/db_completo_idx.1.bt2" ]; then
    echo "Building Bowtie 2 index..."
    ${BOWTIE2_BUILD} --threads 4 "${DATA_DIR}/db_completo.fa" "${DATA_DIR}/db_completo_idx"
else
    echo "Bowtie 2 index already exists."
fi

# 2. Run NextLAMP pipeline
echo "Running NextLAMP design..."
export PYTHONPATH="${BASE_DIR}:${PYTHONPATH}"

${PYTHON_BIN} "${BASE_DIR}/scripts/nextlamp/run_nextlamp.py" \
    --target-fasta "${DATA_DIR}/target_babesia_canis.fa" \
    --index-prefix "${DATA_DIR}/db_completo_idx" \
    --targets-list "${DATA_DIR}/targets_list.txt" \
    --background-list "${DATA_DIR}/background_list.txt" \
    --out "${BASE_DIR}/nextlamp_success.json" \
    --threads 4 \
    --max-sets 10

echo "NextLAMP job finished!"
date
