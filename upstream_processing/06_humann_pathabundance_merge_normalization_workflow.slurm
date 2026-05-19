#!/bin/bash
#SBATCH --job-name=humann_path_merge
#SBATCH --account=project_2012806
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/project_2012191/chifang1/HUMAnN_pathabundance/log/humann_path_merge_%j.out
#SBATCH --error=/scratch/project_2012191/chifang1/HUMAnN_pathabundance/log/humann_path_merge_%j.err

# ============================================================
# HUMAnN pathway-abundance table collection, merging, and normalization
# ============================================================
# This script documents the workflow used to generate the merged
# HUMAnN pathway-abundance table for downstream functional analyses
# in the Parkinson's disease gut metagenomic meta-analysis.
#
# HUMAnN was run separately for each KneadData-filtered metagenomic
# sample. This script collects sample-level HUMAnN pathabundance files
# from multiple output directories, standardizes their filenames,
# merges them into a pathway-by-sample abundance matrix using R, and
# normalizes the merged table to relative abundance using
# humann_renorm_table.
#
# Main processing steps:
#   1. Collect sample-level HUMAnN pathabundance files
#   2. Merge sample-level pathway tables into a single matrix
#   3. Normalize the merged pathway-abundance table to relative abundance
#
# Main outputs:
#   - final_merged_pathabundance.tsv
#   - final_merged_pathabundance_relab.tsv
#
# The normalized pathway-abundance table was used for downstream
# functional profile analysis, pathway-level association testing,
# and figure generation in the manuscript.
# ============================================================

set -euo pipefail

# -------------------------
# 1. Define directories
# -------------------------
PROJECT_DIR="/scratch/project_2012191/chifang1"

WORK_DIR="${PROJECT_DIR}/HUMAnN_pathabundance"
COLLECTED_DIR="${WORK_DIR}/collected_pathabundance"
MERGED_DIR="${WORK_DIR}/merged_pathabundance"
SCRIPT_DIR="${WORK_DIR}/scripts"
LOG_DIR="${WORK_DIR}/log"

mkdir -p "${COLLECTED_DIR}" "${MERGED_DIR}" "${SCRIPT_DIR}" "${LOG_DIR}"

# -------------------------
# 2. Source directories containing HUMAnN pathabundance outputs
# -------------------------
# These directories correspond to HUMAnN outputs generated in different
# sample batches. File names may differ slightly across batches
# depending on whether the HUMAnN input was named "kneaddata" or
# "combined".

SOURCE_DIRS=(
    "/scratch/project_2012191/chifang1/HUMAnN/1_4318_4390"
    "/scratch/project_2012191/chifang1/HUMAnN/1_Results_5523_6030"
    "/scratch/project_2012191/chifang1/HUMAnN_4316_5522"
    "/scratch/project_2012191/chifang1/HUMAnN_path/8_pathabundance_5112_5522"
)

# -------------------------
# 3. Collect and standardize pathabundance files
# -------------------------
echo "Collecting HUMAnN pathabundance files..."
echo "Collection directory: ${COLLECTED_DIR}"

for SRC_DIR in "${SOURCE_DIRS[@]}"; do
    if [[ ! -d "${SRC_DIR}" ]]; then
        echo "Warning: source directory not found: ${SRC_DIR}"
        continue
    fi

    echo "Scanning source directory: ${SRC_DIR}"

    while IFS= read -r FILE; do
        BASENAME=$(basename "${FILE}")

        SAMPLE_ID=$(echo "${BASENAME}" | grep -oE 'SRR[0-9]+' | head -n 1)

        if [[ -z "${SAMPLE_ID}" ]]; then
            echo "Warning: could not extract sample ID from file: ${FILE}"
            continue
        fi

        DEST_FILE="${COLLECTED_DIR}/${SAMPLE_ID}_pathabundance.tsv"

        cp "${FILE}" "${DEST_FILE}"
        echo "Copied ${FILE} -> ${DEST_FILE}"

    done < <(find "${SRC_DIR}" -type f -name "*pathabundance.tsv")
done

N_FILES=$(find "${COLLECTED_DIR}" -type f -name "*_pathabundance.tsv" | wc -l)
echo "Number of collected pathabundance files: ${N_FILES}"

# -------------------------
# 4. Merge sample-level pathabundance files
# -------------------------
MERGE_R_SCRIPT="${SCRIPT_DIR}/05_merge_humann_pathabundance_tables.R"
MERGED_RAW="${MERGED_DIR}/final_merged_pathabundance.tsv"

echo "Merging sample-level HUMAnN pathabundance files..."
echo "Input directory: ${COLLECTED_DIR}"
echo "Merged output: ${MERGED_RAW}"

source /scratch/project_2012191/chifang1/00_Software/miniconda3/etc/profile.d/conda.sh
conda activate r_env

Rscript "${MERGE_R_SCRIPT}" "${COLLECTED_DIR}" "${MERGED_RAW}"

# -------------------------
# 5. Normalize merged pathway-abundance table
# -------------------------
MERGED_RELAB="${MERGED_DIR}/final_merged_pathabundance_relab.tsv"

echo "Normalizing merged pathway-abundance table to relative abundance..."
echo "Normalized output: ${MERGED_RELAB}"

conda activate humann4_env

humann_renorm_table \
    --input "${MERGED_RAW}" \
    --output "${MERGED_RELAB}" \
    --units relab

echo "HUMAnN pathway-abundance merging and normalization completed."
echo "Final merged table: ${MERGED_RAW}"
echo "Final relative-abundance table: ${MERGED_RELAB}"
echo "End time: $(date)"
