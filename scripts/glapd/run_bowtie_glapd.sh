#!/bin/bash
#SBATCH --job-name=glapd_run
#SBATCH --partition=fast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=glapd_run_%j.log
#SBATCH --error=glapd_run_%j.err

echo "Starting SLURM job for Bowtie index building and GLAPD..."
date

BASE_DIR="/home/fabiano.menegidio/workdir/Omics/genomics/babesia"
DATA_DIR="${BASE_DIR}/data"
GLAPD_DIR="${BASE_DIR}/GLAPD"

# 1. Build Bowtie index (CPU and Memory intensive)
echo "Building Bowtie index..."
${GLAPD_DIR}/bowtie/bowtie-build ${DATA_DIR}/db_completo.fa ${DATA_DIR}/db_completo_idx

# 2. Run GLAPD Single to find candidate single primers
echo "Running GLAPD Single..."
cd ${GLAPD_DIR}
./Single -in ${DATA_DIR}/target_babesia_canis.fa -out Bcanis

# 3. Run par.pl (Bowtie alignments against target & background)
echo "Running par.pl (Bowtie alignments)..."
# We increase threads to 4 for bowtie alignment
perl par.pl --in Bcanis --ref ${DATA_DIR}/target_babesia_canis.fa --common ${DATA_DIR}/targets_list.txt --specific ${DATA_DIR}/background_list.txt --bowtie ./bowtie/bowtie --index ${DATA_DIR}/db_completo_idx --threads 4

# 4. Run LAMP to combine primers and generate final specific set
echo "Running GLAPD LAMP..."
./LAMP -in Bcanis -ref ${DATA_DIR}/target_babesia_canis.fa -out ${BASE_DIR}/success_Bcanis.txt -common -specific

echo "GLAPD run finished!"
date
