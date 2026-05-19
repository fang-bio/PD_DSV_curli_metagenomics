#!/bin/bash
#SBATCH --job-name=extract_curli_genes
#SBATCH --output=/scratch/project_2012191/chifang1/R_studio/curli/1206_relab/1206_curli_UniRef90_%j.out
#SBATCH --error=/scratch/project_2012191/chifang1/R_studio/curli/1206_relab/1206_curli_UniRef90_%j.err
#SBATCH --account=project_2012806
#SBATCH --partition=small
#SBATCH --time=5:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=20

############################################################
# Extract E. coli curli-related UniRef90 gene families
#
# Aim:
# Extract E. coli curli-gene abundance from HUMAnN gene-family
# relative-abundance profiles using a curated UniRef90 ID list
# corresponding to csgA-csgG curli genes.
#
# Rationale:
# Experimental evidence linking bacterial curli amyloid fibers to
# alpha-synuclein aggregation has primarily focused on E. coli curli.
# Therefore, this workflow constructs a targeted E. coli curli-gene
# abundance feature for downstream analyses.
#
# Main steps:
# 1. Clean the curated UniRef90 ID list.
# 2. Iterate over sample-level HUMAnN gene-family relative-abundance files.
# 3. Extract exact UniRef90 matches for curli-related genes.
# 4. Export sample-gene-level curli abundance.
#
# Input:
# 1. *_rel_abundance.tsv files
# 2. curli_UniRef90_IDs_cleaned.csv
#
# Output:
# 1206_relab_curli_all_samples.csv
############################################################

set -euo pipefail

# ----------------------------------------------------------
# 1. Path configuration
# ----------------------------------------------------------

OUTPUT_DIR="/scratch/project_2012191/chifang1/R_studio/curli/1206_relab"
INPUT_DIR="/scratch/project_2012191/chifang1/R_studio/curli/genefamilies_relative_abundance"

CUR_LI_ID_FILE="$OUTPUT_DIR/curli_UniRef90_IDs_cleaned.csv"
OUTPUT_FILE="$OUTPUT_DIR/1206_relab_curli_all_samples.csv"
LOG_FILE="$OUTPUT_DIR/curli_UniRef90_log.txt"
SUMMARY_FILE="$OUTPUT_DIR/curli_UniRef90_extraction_summary.txt"

mkdir -p "$OUTPUT_DIR"

echo "Start curli-gene extraction: $(date)" > "$LOG_FILE"
echo "Input directory: $INPUT_DIR" | tee -a "$LOG_FILE"
echo "Curli UniRef90 ID file: $CUR_LI_ID_FILE" | tee -a "$LOG_FILE"
echo "Output file: $OUTPUT_FILE" | tee -a "$LOG_FILE"

# ----------------------------------------------------------
# 2. Validate inputs
# ----------------------------------------------------------

if [ ! -d "$INPUT_DIR" ]; then
  echo "ERROR: Input directory not found: $INPUT_DIR" | tee -a "$LOG_FILE" >&2
  exit 1
fi

if [ ! -s "$CUR_LI_ID_FILE" ]; then
  echo "ERROR: Curli UniRef90 ID file not found or empty: $CUR_LI_ID_FILE" | tee -a "$LOG_FILE" >&2
  exit 1
fi

# Remove Windows CRLF line endings and empty lines
sed -i 's/\r$//' "$CUR_LI_ID_FILE"
sed -i '/^[[:space:]]*$/d' "$CUR_LI_ID_FILE"

N_IDS=$(wc -l < "$CUR_LI_ID_FILE")
echo "Curated curli UniRef90 IDs: $N_IDS" | tee -a "$LOG_FILE"

shopt -s nullglob
input_files=("$INPUT_DIR"/*_rel_abundance.tsv)

if [ "${#input_files[@]}" -eq 0 ]; then
  echo "ERROR: No *_rel_abundance.tsv files found in $INPUT_DIR" | tee -a "$LOG_FILE" >&2
  exit 1
fi

echo "Number of relative-abundance files: ${#input_files[@]}" | tee -a "$LOG_FILE"

# ----------------------------------------------------------
# 3. Extract exact UniRef90 matches
# ----------------------------------------------------------

echo "Sample,Gene,Abundance" > "$OUTPUT_FILE"

processed_count=0
matched_sample_count=0
total_match_count=0

for file in "${input_files[@]}"; do

  sample=$(basename "$file" | sed 's/_rel_abundance.tsv//')
  processed_count=$((processed_count + 1))

  match_count=$(awk '
    BEGIN {
      FS = "\t"
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

      # Exact UniRef90 match only.
      if (gene in ids) {
        count += 1
      }
    }

    END {
      print count + 0
    }
  ' "$CUR_LI_ID_FILE" "$file")

  if [ "$match_count" -gt 0 ]; then

    awk -v sample="$sample" '
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
        abundance = $2

        if (gene in ids) {
          print sample, gene, abundance
        }
      }
    ' "$CUR_LI_ID_FILE" "$file" >> "$OUTPUT_FILE"

    matched_sample_count=$((matched_sample_count + 1))
    total_match_count=$((total_match_count + match_count))

    echo "Processed: $sample; matched curli genes: $match_count" | tee -a "$LOG_FILE"

  else

    echo "Processed: $sample; no matching curli genes detected" | tee -a "$LOG_FILE"

  fi

done

# ----------------------------------------------------------
# 4. Export extraction summary
# ----------------------------------------------------------

{
  echo "Curli UniRef90 extraction summary"
  echo "Date: $(date)"
  echo ""
  echo "Curated UniRef90 IDs: $N_IDS"
  echo "Relative-abundance files processed: $processed_count"
  echo "Samples with curli-gene matches: $matched_sample_count"
  echo "Total sample-gene matches: $total_match_count"
  echo ""
  echo "Input directory:"
  echo "$INPUT_DIR"
  echo ""
  echo "Curli UniRef90 ID file:"
  echo "$CUR_LI_ID_FILE"
  echo ""
  echo "Output file:"
  echo "$OUTPUT_FILE"
} > "$SUMMARY_FILE"

echo "Curli-gene extraction completed: $(date)" | tee -a "$LOG_FILE"
echo "Summary written to: $SUMMARY_FILE" | tee -a "$LOG_FILE"
