#!/bin/bash
#SBATCH --job-name=glapd_prepare
#SBATCH --partition=fast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=glapd_prepare_%j.log
#SBATCH --error=glapd_prepare_%j.err

echo "Starting GLAPD genome preparation pipeline..."
date

# Run python script with absolute path
/home/fabiano.menegidio/miniforge3/bin/python /home/fabiano.menegidio/workdir/Omics/genomics/babesia/scripts/data_prep/process_genomes.py

echo "Pipeline finished!"
date
