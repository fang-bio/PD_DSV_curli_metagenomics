#!/bin/bash
#SBATCH --job-name=curli_extract
#SBATCH --account=project_2012806
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=76G
#SBATCH --output=curli_extract_%j.out
#SBATCH --error=curli_extract_%j.err

############################################################
# E. coli curli-gene extraction from HUMAnN genefamilies
#
# Aim:
# Extract curli-gene abundance from HUMAnN genefamilies output
# using curated UniRef90 identifiers corresponding to Escherichia
# coli csgA-csgG curli genes.
#
# Rationale:
# Experimental evidence linking bacterial curli amyloid fibers to
# alpha-synuclein aggregation has primarily focused on E. coli curli.
# Therefore, this workflow uses curated UniRef90 identifiers for
# E. coli csgA-csgG genes to construct a targeted curli-gene
# abundance feature for downstream analyses.
#
# Main steps:
# 1. Combine and deduplicate curated UniRef90 identifiers for
#    E. coli csgA-csgG genes.
# 2. Normalize line endings in the UniRef90 ID list.
# 3. Extract matching UniRef90 gene families from HUMAnN
#    genefamilies output across samples.
# 4. Use unstratified HUMAnN gene-family rows to avoid double
#    counting taxonomically stratified entries.
# 5. Generate sample-gene-level and sample-level curli abundance
#    tables for downstream statistical analyses.
#
# Input:
# 1. Curated UniRef90 ID CSV files for E. coli csgA-csgG genes.
# 2. HUMAnN *_combined_genefamilies.tsv files.
#
# Output:
# 1. curli_UniRef90_IDs_cleaned.csv
# 2. curli_gene_level_abundance.csv
# 3. curli_sample_level_abundance.csv
# 4. curli_extraction_summary.txt
############################################################

set -euo pipefail

# ==========================================================
# 1. Path configuration
# ==========================================================

BASE_DIR="/scratch/project_2012191/chifang1/R_studio/curli"
ID_DIR="${BASE_DIR}/csgA_G"

GENEFAMILY_DIR="/scratch/project_2012191/chifang1/HUMAnN/1206_genefamilies_all"

OUT_DIR="${BASE_DIR}/1206/curli_extracted_data"
TMP_DIR="${BASE_DIR}/1206/tmp_curli_extraction"

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

RAW_ID_FILE="${OUT_DIR}/curli_UniRef90_IDs_raw.csv"
CLEAN_ID_FILE="${OUT_DIR}/curli_UniRef90_IDs_cleaned.csv"

GENE_LEVEL_RAW="${TMP_DIR}/curli_gene_level_raw.csv"
GENE_LEVEL_OUT="${OUT_DIR}/curli_gene_level_abundance.csv"
SAMPLE_LEVEL_OUT="${OUT_DIR}/curli_sample_level_abundance.csv"
GENE_SUMMARY_OUT="${OUT_DIR}/curli_gene_summary.csv"
SUMMARY_OUT="${OUT_DIR}/curli_extraction_summary.txt"

echo "Start time: $(date)"
echo "ID directory: ${ID_DIR}"
echo "HUMAnN genefamily directory: ${GENEFAMILY_DIR}"
echo "Output directory: ${OUT_DIR}"

# ==========================================================
# 2. Validate input directories
# ==========================================================

if [ ! -d "${ID_DIR}" ]; then
  echo "ERROR: UniRef90 ID directory not found: ${ID_DIR}" >&2
  exit 1
fi

if [ ! -d "${GENEFAMILY_DIR}" ]; then
  echo "ERROR: HUMAnN genefamily directory not found: ${GENEFAMILY_DIR}" >&2
  exit 1
fi

# ==========================================================
# 3. Combine and clean curated UniRef90 IDs
# ==========================================================

echo "Combining curated UniRef90 ID files..."

find "${ID_DIR}" \
  -maxdepth 1 \
  -type f \
  -name "*.csv" \
  ! -name "curli_UniRef90_IDs_raw.csv" \
  ! -name "curli_UniRef90_IDs_cleaned.csv" \
  -print0 | xargs -0 cat > "${RAW_ID_FILE}"

if [ ! -s "${RAW_ID_FILE}" ]; then
  echo "ERROR: No input UniRef90 ID content found in ${ID_DIR}" >&2
  exit 1
fi

grep -o "UniRef90_[A-Za-z0-9_]*" "${RAW_ID_FILE}" | \
  tr -d "\r" | \
  sort -u > "${CLEAN_ID_FILE}"

if [ ! -s "${CLEAN_ID_FILE}" ]; then
  echo "ERROR: Cleaned UniRef90 ID list is empty." >&2
  exit 1
fi

N_IDS=$(wc -l < "${CLEAN_ID_FILE}")
echo "Curated UniRef90 IDs retained: ${N_IDS}"

# ==========================================================
# 4. Extract curli-related gene families across samples
# ==========================================================

echo "Extracting curli-related UniRef90 gene families from HUMAnN outputs..."

echo "Sample,Gene,Abundance" > "${GENE_LEVEL_RAW}"

shopt -s nullglob
GENEFAMILY_FILES=("${GENEFAMILY_DIR}"/*_combined_genefamilies.tsv)

if [ "${#GENEFAMILY_FILES[@]}" -eq 0 ]; then
  echo "ERROR: No *_combined_genefamilies.tsv files found in ${GENEFAMILY_DIR}" >&2
  exit 1
fi

for file in "${GENEFAMILY_FILES[@]}"; do

  sample=$(basename "${file}" "_combined_genefamilies.tsv")

  awk -v sample="${sample}" '
    BEGIN {
      FS = "\t"
      OFS = ","
    }

    NR == FNR {
      ids[$1] = 1
      next
    }

    FNR == 1 {
      next
    }

    {
      gene = $1

      # Use unstratified HUMAnN rows only.
      # Stratified rows contain "|" and are skipped to avoid double counting.
      if (gene ~ /\|/) {
        next
      }

      if (gene in ids) {
        abundance = $2
        print sample, gene, abundance
      }
    }
  ' "${CLEAN_ID_FILE}" "${file}" >> "${GENE_LEVEL_RAW}"

done

N_FILES="${#GENEFAMILY_FILES[@]}"
echo "Processed HUMAnN genefamily files: ${N_FILES}"

# ==========================================================
# 5. Collapse duplicate sample-gene rows if present
# ==========================================================

echo "Collapsing duplicate sample-gene rows if present..."

{
  echo "Sample,Gene,Abundance"
  awk -F',' '
    NR == 1 {
      next
    }

    {
      key = $1 "," $2
      abundance[key] += $3
    }

    END {
      for (key in abundance) {
        print key "," abundance[key]
      }
    }
  ' "${GENE_LEVEL_RAW}" | sort
} > "${GENE_LEVEL_OUT}"

N_GENE_ROWS=$(( $(wc -l < "${GENE_LEVEL_OUT}") - 1 ))
echo "Sample-gene abundance rows retained: ${N_GENE_ROWS}" 

# ==========================================================
# 6. Generate sample-level curli abundance
# ==========================================================

echo "Generating sample-level curli abundance table..."

{
  echo "Sample,Curli_gene_abundance,Number_of_curli_gene_hits"
  awk -F',' '
    NR == 1 {
      next
    }

    {
      sample = $1
      abundance[sample] += $3
      hits[sample] += 1
    }

    END {
      for (sample in abundance) {
        print sample "," abundance[sample] "," hits[sample]
      }
    }
  ' "${GENE_LEVEL_OUT}" | sort
} > "${SAMPLE_LEVEL_OUT}"

N_SAMPLES_WITH_CURli=$(( $(wc -l < "${SAMPLE_LEVEL_OUT}") - 1 ))
echo "Samples with at least one curli-gene hit: ${N_SAMPLES_WITH_CURli}"

# ==========================================================
# 7. Generate gene-level summary
# ==========================================================

echo "Generating gene-level curli summary..."

{
  echo "Gene,Total_abundance,Number_of_samples_detected"
  awk -F',' '
    NR == 1 {
      next
    }

    {
      gene = $2
      total[gene] += $3
      samples[gene "," $1] = 1
    }

    END {
      for (key in samples) {
        split(key, arr, ",")
        gene = arr[1]
        sample_count[gene] += 1
      }

      for (gene in total) {
        print gene "," total[gene] "," sample_count[gene]
      }
    }
  ' "${GENE_LEVEL_OUT}" | sort
} > "${GENE_SUMMARY_OUT}"

# ==========================================================
# 8. Export extraction summary
# ==========================================================

{
  echo "Curli-gene extraction summary"
  echo "Date: $(date)"
  echo ""
  echo "Curated UniRef90 IDs: ${N_IDS}"
  echo "Processed HUMAnN genefamily files: ${N_FILES}"
  echo "Sample-gene abundance rows: ${N_GENE_ROWS}"
  echo "Samples with curli-gene hits: ${N_SAMPLES_WITH_CURli}"
  echo ""
  echo "Input UniRef90 ID directory:"
  echo "${ID_DIR}"
  echo ""
  echo "Input HUMAnN genefamily directory:"
  echo "${GENEFAMILY_DIR}"
  echo ""
  echo "Output files:"
  echo "${CLEAN_ID_FILE}"
  echo "${GENE_LEVEL_OUT}"
  echo "${SAMPLE_LEVEL_OUT}"
  echo "${GENE_SUMMARY_OUT}"
} > "${SUMMARY_OUT}"

echo "Curli-gene extraction completed successfully."
echo "Clean UniRef90 ID list: ${CLEAN_ID_FILE}"
echo "Sample-gene-level output: ${GENE_LEVEL_OUT}"
echo "Sample-level output: ${SAMPLE_LEVEL_OUT}"
echo "Gene-level summary: ${GENE_SUMMARY_OUT}"
echo "Summary: ${SUMMARY_OUT}"
echo "End time: $(date)"
