#!/bin/bash
#SBATCH --job-name=nextlamp_benchmark
#SBATCH --output=/home/fabiano.menegidio/workdir/Omics/genomics/babesia/slurm_logs/benchmark_%j.log
#SBATCH --error=/home/fabiano.menegidio/workdir/Omics/genomics/babesia/slurm_logs/benchmark_%j.err
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00

set -e

# Resolve script directory dynamically (FAIR & Location-Agnostic)
BASE_DIR="/home/fabiano.menegidio/workdir/Omics/genomics/babesia"

mkdir -p "${BASE_DIR}/slurm_logs" "${BASE_DIR}/tests/results"

echo "=========================================================="
echo "  SLURM NextLAMP vs GLAPD Benchmark Job"
echo "  Job ID: $SLURM_JOB_ID"
echo "  Node:   $SLURMD_NODENAME"
echo "  Date:   $(date -Iseconds)"
echo "=========================================================="

# Activate Conda environment if available
if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
    conda activate nextlamp 2>/dev/null || conda activate humann3_env 2>/dev/null || true
fi

export PYTHONPATH="$BASE_DIR:$PYTHONPATH"

# Run Python Benchmark Script
python scripts/benchmark_slurm.py

echo "=========================================================="
echo "  Benchmark finished successfully!"
echo "=========================================================="
