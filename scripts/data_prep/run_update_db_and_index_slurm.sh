#!/bin/bash
#SBATCH --job-name=update_db_index
#SBATCH --partition=medium
#SBATCH --nodelist=n010
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --output=/home/fabiano.menegidio/workdir/Omics/genomics/babesia/slurm_logs/update_db_index_%j.log
#SBATCH --error=/home/fabiano.menegidio/workdir/Omics/genomics/babesia/slurm_logs/update_db_index_%j.err

set -e

BASE_DIR="/home/fabiano.menegidio/workdir/Omics/genomics/babesia"
DATA_DIR="${BASE_DIR}/data"
# Activate environment if available
if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
    conda activate ncbi_env 2>/dev/null || conda activate humann3_env 2>/dev/null || true
fi

BOWTIE2_BUILD="$(which bowtie2-build 2>/dev/null || echo '/home/fabiano.menegidio/miniforge3/envs/ncbi_env/bin/bowtie2-build')"

INDEX_DIR="${DATA_DIR}/indexes"
mkdir -p "${INDEX_DIR}"

echo "--- Step 1: Sorting Segmented FASTA Datasets ---"
python3 scripts/data_prep/process_genomes.py

echo "--- Step 2: Building Segmented Bowtie2 Indexes ---"
echo "Using Bowtie2-build at: ${BOWTIE2_BUILD}"

THREADS=${SLURM_CPUS_PER_TASK:-8}

# Helper function to build index if fasta exists
build_idx() {
    local fa="$1"
    local prefix="$2"
    local name="$3"
    if [ -f "$fa" ]; then
        echo "[Index Build] Building index for ${name} (${fa})..."
        ${BOWTIE2_BUILD} --threads ${THREADS} "${fa}" "${prefix}"
    else
        echo "[Index Skip] ${fa} not found, skipping."
    fi
}

build_idx "${DATA_DIR}/host_dog.fa" "${INDEX_DIR}/idx_dog" "Dog Host"
build_idx "${DATA_DIR}/host_cat.fa" "${INDEX_DIR}/idx_cat" "Cat Host"
build_idx "${DATA_DIR}/host_human.fa" "${INDEX_DIR}/idx_human" "Human Host"
build_idx "${DATA_DIR}/vectors_ticks.fa" "${INDEX_DIR}/idx_ticks" "Tick Vectors"
build_idx "${DATA_DIR}/bg_apicomplexa.fa" "${INDEX_DIR}/idx_apicomplexa" "Apicomplexa Background"
build_idx "${DATA_DIR}/target_babesia_canis.fa" "${INDEX_DIR}/idx_targets" "Babesia Target"

echo "=========================================================="
echo "  Segmented Database preparation & Bowtie2 indexing completed!"
echo "  Finished at: $(date -Iseconds)"
echo "=========================================================="
