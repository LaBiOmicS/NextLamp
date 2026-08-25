#!/bin/bash
#SBATCH --job-name=nextlamp_glapd_loop_comp
#SBATCH --partition=fast,long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/home/fabiano.menegidio/workdir/Omics/genomics/babesia/logs/slurm_loop_comparison_%j.log
#SBATCH --error=/home/fabiano.menegidio/workdir/Omics/genomics/babesia/logs/slurm_loop_comparison_%j.err

echo "=========================================================="
echo "  SLURM Job: NextLAMP vs GLAPD Comparison with Loop Primers"
echo "  Job ID: $SLURM_JOB_ID"
echo "  Node:   $SLURMD_NODENAME"
echo "  Date:   $(date -Iseconds)"
echo "=========================================================="

BASE_DIR="/home/fabiano.menegidio/workdir/Omics/genomics/babesia"
cd "${BASE_DIR}"

export PYTHONPATH="${BASE_DIR}:${PYTHONPATH}"

# Execute NextLAMP vs GLAPD Comparison (with Loop primers support)
python -m nextlamp.tests.glapd_comparison.run_comparison

echo "=========================================================="
echo "  Comparison Job Finished Successfully!"
echo "=========================================================="
