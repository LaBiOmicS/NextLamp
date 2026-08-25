#!/bin/bash
#SBATCH --job-name=glapd_spp_cluster
#SBATCH --partition=fast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=cluster_run/glapd_cluster_%j.log
#SBATCH --error=cluster_run/glapd_cluster_%j.err

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${BASE_DIR}/data"
CLUSTER_DIR="${BASE_DIR}/cluster_run"

echo "=========================================================="
echo "Starting cluster GLAPD backup run for Babesia spp. (all Babesia)..."
echo "Started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "=========================================================="

# 1. Setup isolated cluster directory
mkdir -p "${CLUSTER_DIR}"
cd "${CLUSTER_DIR}"

# Create symlinks to GLAPD binaries and scripts to keep run isolated
ln -sf "${BASE_DIR}/GLAPD/Single" ./Single
ln -sf "${BASE_DIR}/GLAPD/LAMP" ./LAMP
ln -sf "${BASE_DIR}/GLAPD/par.pl" ./par.pl
ln -sfn "${BASE_DIR}/GLAPD/bowtie" ./bowtie
ln -sfn "${BASE_DIR}/GLAPD/Par" ./Par


# 2. Build Bowtie 1 indices for each partition if they don't exist
echo "--- Step 1: Checking and building Bowtie 1 indices for partitions ---"
for i in {1..5}; do
    PART_FA="${DATA_DIR}/db_part${i}.fa"
    IDX_PREFIX="${DATA_DIR}/db_part${i}_idx"
    if [ ! -f "${IDX_PREFIX}.1.ebwt" ]; then
        echo "Building index for part ${i}..."
        ./bowtie/bowtie-build "${PART_FA}" "${IDX_PREFIX}"
        echo "Finished building index for part ${i}."
    else
        echo "Index for part ${i} already exists."
    fi
done

# 3. Run GLAPD Single
echo "--- Step 2: Running GLAPD Single ---"
./Single -in "${DATA_DIR}/target_babesia_canis.fa" -out Bcanis_glapd_cluster

# 4. Run par.pl with all partition indices passed as a comma-separated list
echo "--- Step 3: Running par.pl for alignment checks ---"
INDEX_LIST="${DATA_DIR}/db_part1_idx,${DATA_DIR}/db_part2_idx,${DATA_DIR}/db_part3_idx,${DATA_DIR}/db_part4_idx,${DATA_DIR}/db_part5_idx"
perl par.pl --in Bcanis_glapd_cluster \
            --ref "${DATA_DIR}/target_babesia_canis.fa" \
            --specific "${DATA_DIR}/background_list.txt" \
            --bowtie ./bowtie/bowtie \
            --index "${INDEX_LIST}" \
            --threads 6

# 5. Run GLAPD LAMP to compile final primers
echo "--- Step 4: Running GLAPD LAMP ---"
./LAMP -in Bcanis_glapd_cluster \
       -ref "${DATA_DIR}/target_babesia_canis.fa" \
       -out "${CLUSTER_DIR}/success_glapd_cluster.txt" \
       -specific

echo "=========================================================="
echo "GLAPD cluster backup run completed successfully!"
echo "Ended at: $(date)"
echo "Results saved to ${CLUSTER_DIR}/success_glapd_cluster.txt"
echo "=========================================================="
