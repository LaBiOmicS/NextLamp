#!/bin/bash
#SBATCH --job-name=glapd_specific
#SBATCH --partition=fast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=glapd_specific_%j.log
#SBATCH --error=glapd_specific_%j.err

echo "Starting GLAPD Specificity Run on Babesia canis..."
date

BASE_DIR="/home/fabiano.menegidio/workdir/Omics/genomics/babesia"
DATA_DIR="${BASE_DIR}/data"
GLAPD_DIR="${BASE_DIR}/GLAPD"

# 1. Build Bowtie index if not already built
if [ ! -f "${DATA_DIR}/db_completo_idx.1.ebwt" ]; then
    echo "Building Bowtie index..."
    ${GLAPD_DIR}/bowtie/bowtie-build ${DATA_DIR}/db_completo.fa ${DATA_DIR}/db_completo_idx
else
    echo "Bowtie index already exists."
fi

# 2. Run GLAPD Single
echo "Running GLAPD Single..."
cd ${GLAPD_DIR}
./Single -in ${DATA_DIR}/target_babesia_canis.fa -out Bcanis_glapd

# 3. Run par.pl with specificity only (resolves target contig commonality issue)
echo "Running par.pl for specificity checks..."
perl par.pl --in Bcanis_glapd --ref ${DATA_DIR}/target_babesia_canis.fa --specific ${DATA_DIR}/background_list.txt --bowtie ./bowtie/bowtie --index ${DATA_DIR}/db_completo_idx --threads 4

# 4. Run LAMP with specificity filter
echo "Running GLAPD LAMP..."
./LAMP -in Bcanis_glapd -ref ${DATA_DIR}/target_babesia_canis.fa -out ${BASE_DIR}/success_glapd.txt -specific

echo "GLAPD Specificity Run finished!"
date
