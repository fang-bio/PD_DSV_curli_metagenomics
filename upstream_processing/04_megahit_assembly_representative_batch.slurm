#!/bin/bash
#SBATCH --job-name=megahit_5401_5522
#SBATCH --account=project_2012806
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=372G
#SBATCH --partition=small
#SBATCH --output=/scratch/project_2012191/chifang1/contigs/log/megahit_5401_5522_%j.out
#SBATCH --error=/scratch/project_2012191/chifang1/contigs/log/megahit_5401_5522_%j.err

# ============================================================
# Representative MEGAHIT metagenomic assembly batch
# ============================================================
# This script shows one representative SLURM batch used for
# metagenomic assembly of KneadData-filtered shotgun metagenomic
# reads using MEGAHIT.
#
# The samples were processed in multiple SLURM batches on Puhti for
# scheduling and resource-management purposes. All batches used the
# same MEGAHIT workflow, software environment, input structure, and
# assembly parameters; only the sample ID range was changed between
# batches.
#
# This representative batch processes samples SRR19065401 to
# SRR19065522.
#
# Input:
#   KneadData-filtered paired-end FASTQ files:
#   - ${SAMPLE}_1_kneaddata_paired_1.fastq
#   - ${SAMPLE}_1_kneaddata_paired_2.fastq
#
# Main processing steps:
#   - Per-sample metagenomic assembly using MEGAHIT
#   - Collection of final.contigs.fa files
#   - Removal of intermediate assembly directories
#
# Main outputs:
#   - Per-sample assembled contigs:
#     ${SAMPLE}_final.contigs.fa
#
# Note:
#   Include this script in the manuscript repository only if
#   assembly-derived contigs were used in the final analyses.
# ============================================================

set -euo pipefail

# -------------------------
# 1. Software environment
# -------------------------
module load biokit

# -------------------------
# 2. Input and output directories
# -------------------------
INPUT_DIR="/scratch/project_2012191/chifang1/tmp/QC_deHost"
OUTPUT_DIR="/scratch/project_2012191/chifang1/contigs"
FINAL_DIR="${OUTPUT_DIR}/final_contigs"
LOG_DIR="${OUTPUT_DIR}/log"

mkdir -p "${FINAL_DIR}" "${LOG_DIR}"

# -------------------------
# 3. Define sample range for this representative batch
# -------------------------
# This batch processes samples SRR19065401 to SRR19065522.
# Other batches used the same workflow but different sample ID ranges.

SAMPLES=()
for i in $(seq 5401 5522); do
    SAMPLES+=("SRR1906${i}")
done

# -------------------------
# 4. Process each sample
# -------------------------
for SAMPLE in "${SAMPLES[@]}"; do

    echo "============================================================"
    echo "Processing sample: ${SAMPLE}"
    echo "Start time: $(date)"
    echo "============================================================"

    FORWARD_READ="${INPUT_DIR}/${SAMPLE}_1_kneaddata_paired_1.fastq"
    REVERSE_READ="${INPUT_DIR}/${SAMPLE}_1_kneaddata_paired_2.fastq"
    SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE}_output"
    SAMPLE_LOG="${LOG_DIR}/${SAMPLE}_megahit.log"
    FINAL_CONTIGS="${FINAL_DIR}/${SAMPLE}_final.contigs.fa"

    # Skip samples that have already been assembled
    if [[ -f "${FINAL_CONTIGS}" ]]; then
        echo "Sample ${SAMPLE} has already been assembled. Skipping."
        continue
    fi

    # Check input files
    if [[ ! -f "${FORWARD_READ}" ]]; then
        echo "Forward read file not found: ${FORWARD_READ}" | tee -a "${SAMPLE_LOG}"
        continue
    fi

    if [[ ! -f "${REVERSE_READ}" ]]; then
        echo "Reverse read file not found: ${REVERSE_READ}" | tee -a "${SAMPLE_LOG}"
        continue
    fi

    # Remove existing incomplete output directory if present
    if [[ -d "${SAMPLE_OUTPUT_DIR}" ]]; then
        echo "Removing existing output directory: ${SAMPLE_OUTPUT_DIR}" | tee -a "${SAMPLE_LOG}"
        rm -rf "${SAMPLE_OUTPUT_DIR}"
    fi

    # Run MEGAHIT assembly
    echo "Running MEGAHIT assembly for sample ${SAMPLE}..." | tee -a "${SAMPLE_LOG}"

    srun megahit \
        -1 "${FORWARD_READ}" \
        -2 "${REVERSE_READ}" \
        -t "${SLURM_CPUS_PER_TASK}" \
        --mem-flag 2 \
        -o "${SAMPLE_OUTPUT_DIR}" \
        >> "${SAMPLE_LOG}" 2>&1

    # Check MEGAHIT output
    if [[ ! -f "${SAMPLE_OUTPUT_DIR}/final.contigs.fa" ]]; then
        echo "Error: MEGAHIT failed for sample ${SAMPLE}; final.contigs.fa was not found." | tee -a "${SAMPLE_LOG}"
        continue
    fi

    # Move and rename final contigs
    mv "${SAMPLE_OUTPUT_DIR}/final.contigs.fa" "${FINAL_CONTIGS}"

    # Remove intermediate assembly directory
    rm -rf "${SAMPLE_OUTPUT_DIR}"

    echo "Sample ${SAMPLE} assembly completed successfully." | tee -a "${SAMPLE_LOG}"
    echo "End time: $(date)" | tee -a "${SAMPLE_LOG}"

done

echo "All samples in this representative MEGAHIT batch have been processed."
echo "Script finished at: $(date)"
