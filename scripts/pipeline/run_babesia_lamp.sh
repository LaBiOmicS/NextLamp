#!/bin/bash
#SBATCH --job-name=babesia_lamp
#SBATCH --partition=fast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=babesia_lamp_%j.log
#SBATCH --error=babesia_lamp_%j.err

set -e

echo "=========================================================="
echo "  Babesia spp. LAMP Primer Design — GLAPD + NextLAMP"
echo "  Started at: $(date)"
echo "  Job ID: $SLURM_JOB_ID"
echo "=========================================================="

BASE_DIR="/home/fabiano.menegidio/workdir/Omics/genomics/babesia"
DATA_DIR="${BASE_DIR}/data"
GLAPD_DIR="${BASE_DIR}/GLAPD"
RESULTS_DIR="${BASE_DIR}/results"

BOWTIE1_BUILD="${GLAPD_DIR}/bowtie/bowtie-build"
BOWTIE1="${GLAPD_DIR}/bowtie/bowtie"
BOWTIE2_BUILD="/home/fabiano.menegidio/miniforge3/envs/humann3_env/bin/bowtie2-build"
BOWTIE2="/home/fabiano.menegidio/miniforge3/envs/humann3_env/bin/bowtie2"
PYTHON="/home/fabiano.menegidio/miniforge3/bin/python"

TARGET_FA="${DATA_DIR}/target_babesia_canis.fa"
TARGETS_LIST="${DATA_DIR}/targets_list.txt"
BACKGROUND_LIST="${DATA_DIR}/background_list.txt"

mkdir -p "${RESULTS_DIR}"

##########################################################################
#  PART A: GLAPD (Bowtie 1 — partitioned indices)
##########################################################################
echo ""
echo "=========================================================="
echo "  PART A: GLAPD"
echo "=========================================================="

cd "${GLAPD_DIR}"

# A1. Build Bowtie 1 indices for each partition (if needed)
echo "--- A1: Checking Bowtie 1 indices for partitions ---"
for i in 1 2 3 4 5; do
    PART_FA="${DATA_DIR}/db_part${i}.fa"
    IDX_PREFIX="${DATA_DIR}/db_part${i}_idx"
    if [ ! -f "${IDX_PREFIX}.1.ebwt" ]; then
        echo "Building Bowtie 1 index for part ${i}..."
        ${BOWTIE1_BUILD} "${PART_FA}" "${IDX_PREFIX}"
        echo "Index for part ${i} built."
    else
        echo "Index for part ${i} already exists."
    fi
done

# A2. GLAPD Single
echo "--- A2: Running GLAPD Single ---"
./Single -in "${TARGET_FA}" -out Bcanis_full

# A3. par.pl with specificity (multi-index)
echo "--- A3: Running par.pl (specificity) ---"
INDEX_LIST="${DATA_DIR}/db_part1_idx,${DATA_DIR}/db_part2_idx,${DATA_DIR}/db_part3_idx,${DATA_DIR}/db_part4_idx,${DATA_DIR}/db_part5_idx"
perl par.pl --in Bcanis_full \
            --ref "${TARGET_FA}" \
            --specific "${BACKGROUND_LIST}" \
            --bowtie "${BOWTIE1}" \
            --index "${INDEX_LIST}" \
            --threads 6

# A4. GLAPD LAMP
echo "--- A4: Running GLAPD LAMP ---"
./LAMP -in Bcanis_full \
       -ref "${TARGET_FA}" \
       -out "${RESULTS_DIR}/glapd_babesia_success.txt" \
       -specific

echo "--- GLAPD finished at $(date) ---"

##########################################################################
#  PART B: NextLAMP (Bowtie 2 — large index)
##########################################################################
echo ""
echo "=========================================================="
echo "  PART B: NextLAMP"
echo "=========================================================="

BT2_IDX="${DATA_DIR}/db_completo_bt2_idx"

# B1. Build Bowtie 2 large index (if needed)
if [ ! -f "${BT2_IDX}.1.bt2l" ] && [ ! -f "${BT2_IDX}.1.bt2" ]; then
    echo "--- B1: Building Bowtie 2 large index ---"
    ${BOWTIE2_BUILD} --large-index "${DATA_DIR}/db_completo.fa" "${BT2_IDX}"
    echo "Bowtie 2 index built."
else
    echo "--- B1: Bowtie 2 index already exists ---"
fi

# B2. Run NextLAMP pipeline
echo "--- B2: Running NextLAMP pipeline ---"
export PYTHONPATH="${BASE_DIR}:${PYTHONPATH}"

${PYTHON} "${BASE_DIR}/scripts/nextlamp/run_nextlamp.py" \
    --target-fasta "${TARGET_FA}" \
    --bowtie2-path "${BOWTIE2}" \
    --index-prefix "${BT2_IDX}" \
    --targets-list "${TARGETS_LIST}" \
    --background-list "${BACKGROUND_LIST}" \
    --out "${RESULTS_DIR}/nextlamp_babesia_success.json" \
    --max-sets 10

echo "--- NextLAMP finished at $(date) ---"

##########################################################################
#  DONE
##########################################################################
echo ""
echo "=========================================================="
echo "  ALL DONE"
echo "  Ended at: $(date)"
echo ""
echo "  Results:"
echo "    GLAPD:    ${RESULTS_DIR}/glapd_babesia_success.txt"
echo "    NextLAMP: ${RESULTS_DIR}/nextlamp_babesia_success.json"
echo "=========================================================="
