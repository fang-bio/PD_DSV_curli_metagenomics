#!/bin/bash

#SBATCH --account=project_2012806
#SBATCH --partition=small
#SBATCH --time=1-00:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=download_PRJDB15956
#SBATCH --output=/scratch/project_2012191/chifang1/Meta_analysis/9_PRJDB15956/logs/PRJDB15956_download_%j.out
#SBATCH --error=/scratch/project_2012191/chifang1/Meta_analysis/9_PRJDB15956/logs/PRJDB15956_download_%j.err
#SBATCH --ntasks=1

# ============================================================
# Raw FASTQ download for public shotgun metagenomic datasets
# ============================================================
# This SLURM script submits the raw-read download step for one
# publicly available fecal shotgun metagenomic cohort included in
# the Parkinson's disease gut microbiome meta-analysis.
#
# The actual download commands are stored in the accompanying
# project-specific FTP script:
#
#   PRJDB15956-fastq_ftp.sh
#
# This script was used to retrieve raw FASTQ files from public
# sequencing repositories before standardized quality control,
# host-read removal, taxonomic profiling, and functional profiling.
#
# Raw sequencing files are not included in the GitHub repository
# because of their large file size. Public accession numbers are
# provided in the manuscript Supplementary Table 1.
# ============================================================

set -euo pipefail

# -------------------------
# 1. Define working directory
# -------------------------
WORK_DIR="/scratch/project_2012191/chifang1/Meta_analysis/9_PRJDB15956"
LOG_DIR="${WORK_DIR}/logs"

mkdir -p "${LOG_DIR}"

cd "${WORK_DIR}"

echo "Starting raw FASTQ download for PRJDB15956"
echo "Working directory: ${WORK_DIR}"
echo "Start time: $(date)"

# -------------------------
# 2. Run project-specific FTP download script
# -------------------------
bash "${WORK_DIR}/PRJDB15956-fastq_ftp.sh"

echo "Finished raw FASTQ download for PRJDB15956"
echo "End time: $(date)"
