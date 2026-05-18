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
# Representative raw FASTQ download script for one public cohort
# ============================================================
# This SLURM script documents one representative cohort-specific
# raw-read download step used in the Parkinson's disease gut
# metagenome meta-analysis.
#
# The full study included multiple publicly available fecal shotgun
# metagenomic cohorts. Raw FASTQ files were downloaded using
# cohort-specific FTP scripts. Across cohorts, the same overall
# download workflow was used; only the project accession, working
# directory, and accession-specific FTP links differed.
#
# This representative script was used for project PRJDB15956.
# The actual FTP download commands are stored in the accompanying
# project-specific script:
#
#   PRJDB15956-fastq_ftp.sh
#
# Downloaded raw FASTQ files were subsequently processed using the
# standardized upstream workflow, including read-level quality control,
# host-read removal, taxonomic profiling, and functional profiling.
#
# Raw sequencing files are not included in this GitHub repository
# because of their large file size and public availability. Complete
# accession numbers and cohort-level information are provided in
# Supplementary Table 1 of the manuscript.
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
