#!/bin/bash
#SBATCH --job-name=parallel_loop_comp
#SBATCH --partition=fast,long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/home/fabiano.menegidio/workdir/Omics/genomics/babesia/logs/slurm_parallel_loop_%j.log
#SBATCH --error=/home/fabiano.menegidio/workdir/Omics/genomics/babesia/logs/slurm_parallel_loop_%j.err

echo "=========================================================="
echo "  SLURM Job: Parallel NextLAMP vs GLAPD Evaluation (Loop Primers)"
echo "  Job ID: $SLURM_JOB_ID"
echo "  Node:   $SLURMD_NODENAME"
echo "  Date:   $(date -Iseconds)"
echo "=========================================================="

BASE_DIR="/home/fabiano.menegidio/workdir/Omics/genomics/babesia"
cd "${BASE_DIR}"

export PYTHONPATH="${BASE_DIR}:${PYTHONPATH}"

# Execute Parallel Loop Comparison
python -m nextlamp.tests.glapd_comparison.run_parallel_loop_evaluation

echo "=========================================================="
echo "  Parallel Loop Comparison Job Finished Successfully!"
echo "=========================================================="
