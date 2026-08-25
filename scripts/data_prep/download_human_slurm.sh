#!/bin/bash
#SBATCH --job-name=download_human
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/home/fabiano.menegidio/workdir/Omics/genomics/babesia/slurm_logs/download_human_%j.log
#SBATCH --error=/home/fabiano.menegidio/workdir/Omics/genomics/babesia/slurm_logs/download_human_%j.err

set -e

BASE_DIR="/home/fabiano.menegidio/workdir/Omics/genomics/babesia"
DATA_DIR="${BASE_DIR}/data"
HUMAN_ZIP="${DATA_DIR}/human.zip"
HUMAN_RAW="${DATA_DIR}/human_raw"

mkdir -p "${BASE_DIR}/slurm_logs"

DATASETS_BIN="/home/fabiano.menegidio/miniforge3/envs/ncbi_env/bin/datasets"
if [ ! -f "${DATASETS_BIN}" ]; then
    DATASETS_BIN="$(which datasets 2>/dev/null || echo 'datasets')"
fi

echo "=========================================================="
echo "  Downloading Human Genome (GRCh38.p14 - GCF_000001405.40)"
echo "  Started at: $(date -Iseconds)"
echo "=========================================================="

if [ ! -d "${HUMAN_RAW}" ]; then
    echo "[1/3] Downloading dehydrated Human genome assembly..."
    ${DATASETS_BIN} download genome accession GCF_000001405.40 --dehydrated --filename "${HUMAN_ZIP}"
    
    echo "[2/3] Extracting Human dataset structure..."
    unzip -q "${HUMAN_ZIP}" -d "${HUMAN_RAW}"
    rm -f "${HUMAN_ZIP}"
    
    echo "[3/3] Rehydrating Human genome (downloading FASTA files)..."
    ${DATASETS_BIN} rehydrate --directory "${HUMAN_RAW}"
    
    echo "[SUCCESS] Human genome downloaded & rehydrated into ${HUMAN_RAW}"
else
    echo "[SKIP] ${HUMAN_RAW} already exists."
fi

echo "=========================================================="
echo "  Finished at: $(date -Iseconds)"
echo "=========================================================="
