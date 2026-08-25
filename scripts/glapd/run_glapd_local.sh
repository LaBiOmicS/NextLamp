#!/bin/bash
# Exit immediately if a command exits with a non-zero status
set -e

echo "=========================================================="
echo "Starting local GLAPD analysis for Babesia spp. (all Babesia)..."
echo "Started at: $(date)"
echo "=========================================================="

BASE_DIR="/home/fabiano.menegidio/workdir/Omics/genomics/babesia"
DATA_DIR="${BASE_DIR}/data"
GLAPD_DIR="${BASE_DIR}/GLAPD"

# 1. Build Bowtie 1 indices for each partition if they don't exist
echo "--- Step 1: Checking and building Bowtie 1 indices for partitions ---"
for i in {1..5}; do
    PART_FA="${DATA_DIR}/db_part${i}.fa"
    IDX_PREFIX="${DATA_DIR}/db_part${i}_idx"
    if [ ! -f "${IDX_PREFIX}.1.ebwt" ]; then
        echo "Building index for part ${i}..."
        ${GLAPD_DIR}/bowtie/bowtie-build "${PART_FA}" "${IDX_PREFIX}"
        echo "Finished building index for part ${i}."
    else
        echo "Index for part ${i} already exists."
    fi
done

# 2. Run GLAPD Single
echo "--- Step 2: Running GLAPD Single ---"
cd "${GLAPD_DIR}"
./Single -in "${DATA_DIR}/target_babesia_canis.fa" -out Bcanis_glapd

# 3. Run par.pl with all partition indices passed as a comma-separated list
echo "--- Step 3: Running par.pl for alignment checks ---"
INDEX_LIST="${DATA_DIR}/db_part1_idx,${DATA_DIR}/db_part2_idx,${DATA_DIR}/db_part3_idx,${DATA_DIR}/db_part4_idx,${DATA_DIR}/db_part5_idx"
perl par.pl --in Bcanis_glapd \
            --ref "${DATA_DIR}/target_babesia_canis.fa" \
            --specific "${DATA_DIR}/background_list.txt" \
            --bowtie ./bowtie/bowtie \
            --index "${INDEX_LIST}" \
            --threads 6

# 4. Run GLAPD LAMP to compile final primers
echo "--- Step 4: Running GLAPD LAMP ---"
./LAMP -in Bcanis_glapd \
       -ref "${DATA_DIR}/target_babesia_canis.fa" \
       -out "${BASE_DIR}/success_glapd.txt" \
       -specific

echo "=========================================================="
echo "GLAPD local analysis completed successfully!"
echo "Ended at: $(date)"
echo "Results saved to ${BASE_DIR}/success_glapd.txt"
echo "=========================================================="
