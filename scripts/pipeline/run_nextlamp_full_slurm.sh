#!/bin/bash
#SBATCH --job-name=nextlamp_full
#SBATCH --partition=fast,long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=/home/fabiano.menegidio/workdir/Omics/genomics/babesia/logs/slurm_nextlamp_full_%j.log
#SBATCH --error=/home/fabiano.menegidio/workdir/Omics/genomics/babesia/logs/slurm_nextlamp_full_%j.err

echo "=========================================================="
echo "  SLURM Job: NextLAMP Full Babesia Dataset Run"
echo "  Job ID: $SLURM_JOB_ID"
echo "  Node:   $SLURMD_NODENAME"
echo "  Date:   $(date -Iseconds)"
echo "=========================================================="

BASE_DIR="/home/fabiano.menegidio/workdir/Omics/genomics/babesia"
cd "${BASE_DIR}"

mkdir -p logs results/babesia_full

bash scripts/pipeline/run_nextlamp_full.sh

echo "=========================================================="
echo "  NextLAMP Full Run Completed!"
echo "=========================================================="
