#!/bin/bash
#SBATCH --job-name=kneaddata_5301_5400
#SBATCH --account=project_2012806
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem=370G
#SBATCH --cpus-per-task=40
#SBATCH --output=kneaddata_5301_5400_%j.out
#SBATCH --error=kneaddata_5301_5400_%j.err

# ============================================================
# Representative KneadData quality-control and host-removal batch
# ============================================================
# This script shows one representative SLURM batch used for read-level
# quality control and host-read removal in the Parkinson's disease gut
# metagenome meta-analysis.
#
# The full dataset contained 1,609 shotgun metagenomic samples and was
# processed in multiple SLURM batches on Puhti for scheduling and
# resource-management purposes. All batches used the same KneadData
# workflow and software environment; only the sample ID range was changed
# between batches.
#
# This representative batch processes samples SRR19065301 to SRR19065400.
#
# Input:
#   Raw paired-end FASTQ files:
#   - ${SAMPLE_ID}_1.fastq.gz
#   - ${SAMPLE_ID}_2.fastq.gz
#
# Main processing steps:
#   - Read-level quality control using KneadData
#   - Host-read removal using a human genome reference database
#   - Removal of intermediate files
#
# Main outputs:
#   - KneadData-filtered paired-end FASTQ files:
#     ${SAMPLE_ID}_1_kneaddata_paired_1.fastq
#     ${SAMPLE_ID}_1_kneaddata_paired_2.fastq
#
# These filtered reads were subsequently used for downstream taxonomic
# profiling with MetaPhlAn and functional profiling with HUMAnN.
# ============================================================

set -euo pipefail

# -------------------------
# 1. Software environment
# -------------------------
source /scratch/project_2012191/chifang1/00_Software/miniconda3/etc/profile.d/conda.sh
conda activate kneaddata_env

# -------------------------
# 2. Input, output, log, and database paths
# -------------------------
INPUT_DIR="/scratch/project_2012191/chifang1/raw_data"
OUTPUT_BASE_DIR="/scratch/project_2012191/chifang1/tmp"
QC_DEHOST_DIR="${OUTPUT_BASE_DIR}/QC_deHost"
LOG_DIR="${OUTPUT_BASE_DIR}/Log"
DATABASE="/scratch/project_2012191/chifang1/00_Databases/Human_genome"
TRIMMOMATIC_PATH="/scratch/project_2012191/chifang1/00_Software/miniconda3/envs/kneaddata_env/share/trimmomatic-0.39-2"

mkdir -p "${QC_DEHOST_DIR}" "${LOG_DIR}"

# -------------------------
# 3. Define sample range for this representative batch
# -------------------------
# This batch processes samples SRR19065301 to SRR19065400.
# Other batches used the same workflow but different sample ID ranges.

SAMPLES=()
for i in $(seq 5301 5400); do
    SAMPLES+=("SRR1906${i}")
done

# -------------------------
# 4. Process each sample
# -------------------------
for SAMPLE_ID in "${SAMPLES[@]}"; do

    echo "============================================================"
    echo "Processing sample: ${SAMPLE_ID}"
    echo "Start time: $(date)"
    echo "============================================================"

    # Skip samples that have already been processed
    if [[ -f "${QC_DEHOST_DIR}/${SAMPLE_ID}_1_kneaddata_paired_1.fastq" ]]; then
        echo "Sample ${SAMPLE_ID} has already been processed. Skipping."
        continue
    fi

    # Define sample-specific temporary output directory
    OUTPUT_DIR="${OUTPUT_BASE_DIR}/${SAMPLE_ID}"
    mkdir -p "${OUTPUT_DIR}"

    # Check input files
    FORWARD_READ="${INPUT_DIR}/${SAMPLE_ID}_1.fastq.gz"
    REVERSE_READ="${INPUT_DIR}/${SAMPLE_ID}_2.fastq.gz"

    if [[ ! -f "${FORWARD_READ}" ]]; then
        echo "Forward read file not found: ${FORWARD_READ}"
        continue
    fi

    if [[ ! -f "${REVERSE_READ}" ]]; then
        echo "Reverse read file not found: ${REVERSE_READ}"
        continue
    fi

    # Run KneadData
    kneaddata \
        --input1 "${FORWARD_READ}" \
        --input2 "${REVERSE_READ}" \
        --reference-db "${DATABASE}" \
        --output "${OUTPUT_DIR}" \
        --trimmomatic "${TRIMMOMATIC_PATH}" \
        --threads 40 \
        --log "${LOG_DIR}/kneaddata_${SAMPLE_ID}.log" \
        --remove-intermediate-output

    # Check whether KneadData successfully generated paired output files
    PAIRED_1="${OUTPUT_DIR}/${SAMPLE_ID}_1_kneaddata_paired_1.fastq"
    PAIRED_2="${OUTPUT_DIR}/${SAMPLE_ID}_1_kneaddata_paired_2.fastq"

    if [[ -f "${PAIRED_1}" && -f "${PAIRED_2}" ]]; then
        echo "Sample ${SAMPLE_ID} processed successfully. Moving output files..."

        mv "${PAIRED_1}" "${QC_DEHOST_DIR}/"
        mv "${PAIRED_2}" "${QC_DEHOST_DIR}/"

        # Remove temporary sample directory
        rm -rf "${OUTPUT_DIR}"

        echo "Sample ${SAMPLE_ID} completed successfully."
    else
        echo "Sample ${SAMPLE_ID} failed. Expected KneadData output files were not found."
        echo "Please check the log file: ${LOG_DIR}/kneaddata_${SAMPLE_ID}.log"
    fi

    echo "End time for ${SAMPLE_ID}: $(date)"
done

echo "All samples in this representative batch have been processed."
echo "Script finished at: $(date)"
