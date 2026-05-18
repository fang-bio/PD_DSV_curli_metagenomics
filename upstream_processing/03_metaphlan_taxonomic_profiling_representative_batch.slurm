#!/bin/bash
#SBATCH --job-name=metaphlan_5101_5200
#SBATCH --output=/scratch/project_2012191/chifang1/metaphlan/log/metaphlan_5101_5200_%j.out
#SBATCH --error=/scratch/project_2012191/chifang1/metaphlan/log/metaphlan_5101_5200_%j.err
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=372G
#SBATCH --partition=small
#SBATCH --account=project_2012806

# ============================================================
# Representative MetaPhlAn taxonomic-profiling batch
# ============================================================
# This script shows one representative SLURM batch used for
# taxonomic profiling of KneadData-filtered shotgun metagenomic
# samples in the Parkinson's disease gut metagenome meta-analysis.
#
# The final analyzed dataset contained 1,609 shotgun metagenomic
# samples after standardized quality control and host-read removal.
# Taxonomic profiling was performed in multiple SLURM batches on
# Puhti for scheduling and resource-management purposes. All batches
# used the same MetaPhlAn workflow, software environment, database
# path, and processing parameters; only the sample ID range was
# changed between batches.
#
# This representative batch processes samples SRR19065101 to
# SRR19065200.
#
# Input:
#   KneadData-filtered paired-end FASTQ files:
#   - SRR${sample_id}_1_kneaddata_paired_1.fastq
#   - SRR${sample_id}_1_kneaddata_paired_2.fastq
#
# Main processing steps:
#   - Taxonomic profiling using MetaPhlAn
#   - Bowtie2 alignment against the MetaPhlAn marker-gene database
#   - Generation of sample-level taxonomic profiles
#
# Main outputs:
#   - MetaPhlAn taxonomic profile files
#   - Bowtie2 alignment output files
#
# These taxonomic profiles were subsequently merged and used for
# downstream community-structure analysis, differential abundance
# testing, prevalence analysis, and figure generation in the manuscript.
# ============================================================

set -euo pipefail

# -------------------------
# 1. Software environment
# -------------------------
source /scratch/project_2012191/chifang1/00_Software/miniconda3/etc/profile.d/conda.sh
conda activate mpa_env

# -------------------------
# 2. Input, output, log, and database paths
# -------------------------
INPUT_DIR="/scratch/project_2012191/chifang1/tmp/QC_deHost"
OUTPUT_TXT_DIR="/scratch/project_2012191/chifang1/metaphlan/final_txt"
OUTPUT_SAM_DIR="/scratch/project_2012191/chifang1/metaphlan/final_sam"
LOG_DIR="/scratch/project_2012191/chifang1/metaphlan/log"
BOWTIE2_DB="/scratch/project_2012191/chifang1/00_Databases/metaphlan_db"

mkdir -p "${OUTPUT_TXT_DIR}" "${OUTPUT_SAM_DIR}" "${LOG_DIR}"

# -------------------------
# 3. Define sample range for this representative batch
# -------------------------
# This batch processes samples SRR19065101 to SRR19065200.
# Other batches used the same workflow but different sample ID ranges.

START_ID=19065101
END_ID=19065200

# -------------------------
# 4. Process each sample
# -------------------------
for ((sample_id=${START_ID}; sample_id<=${END_ID}; sample_id++)); do

    SAMPLE="SRR${sample_id}"

    INPUT_FASTQ1="${INPUT_DIR}/${SAMPLE}_1_kneaddata_paired_1.fastq"
    INPUT_FASTQ2="${INPUT_DIR}/${SAMPLE}_1_kneaddata_paired_2.fastq"
    OUTPUT_TXT="${OUTPUT_TXT_DIR}/${SAMPLE}_metaphlan_output.txt"
    OUTPUT_SAM="${OUTPUT_SAM_DIR}/${SAMPLE}_bowtie2_output.sam"
    SAMPLE_LOG="${LOG_DIR}/${SAMPLE}_metaphlan.log"

    echo "============================================================" | tee -a "${SAMPLE_LOG}"
    echo "Processing sample: ${SAMPLE}" | tee -a "${SAMPLE_LOG}"
    echo "Start time: $(date)" | tee -a "${SAMPLE_LOG}"
    echo "Forward read: ${INPUT_FASTQ1}" | tee -a "${SAMPLE_LOG}"
    echo "Reverse read: ${INPUT_FASTQ2}" | tee -a "${SAMPLE_LOG}"
    echo "============================================================" | tee -a "${SAMPLE_LOG}"

    # Skip samples that have already been profiled
    if [[ -f "${OUTPUT_TXT}" ]]; then
        echo "MetaPhlAn output already exists for ${SAMPLE}. Skipping." | tee -a "${SAMPLE_LOG}"
        continue
    fi

    # Check input files
    if [[ ! -f "${INPUT_FASTQ1}" ]]; then
        echo "Forward read file not found: ${INPUT_FASTQ1}. Skipping ${SAMPLE}." | tee -a "${SAMPLE_LOG}"
        continue
    fi

    if [[ ! -f "${INPUT_FASTQ2}" ]]; then
        echo "Reverse read file not found: ${INPUT_FASTQ2}. Skipping ${SAMPLE}." | tee -a "${SAMPLE_LOG}"
        continue
    fi

    # Run MetaPhlAn
    echo "Running MetaPhlAn for ${SAMPLE}..." | tee -a "${SAMPLE_LOG}"

    if metaphlan "${INPUT_FASTQ1},${INPUT_FASTQ2}" \
        --input_type fastq \
        --output_file "${OUTPUT_TXT}" \
        --bowtie2db "${BOWTIE2_DB}" \
        --bowtie2out "${OUTPUT_SAM}" \
        --nproc "${SLURM_CPUS_PER_TASK}" \
        >> "${SAMPLE_LOG}" 2>&1; then

        echo "MetaPhlAn completed successfully for ${SAMPLE}." | tee -a "${SAMPLE_LOG}"
    else
        echo "MetaPhlAn failed for ${SAMPLE}. Please check the log file: ${SAMPLE_LOG}" | tee -a "${SAMPLE_LOG}"
        continue
    fi

    echo "End time for ${SAMPLE}: $(date)" | tee -a "${SAMPLE_LOG}"

done

echo "All samples in this representative MetaPhlAn batch have been processed."
echo "Script finished at: $(date)"
