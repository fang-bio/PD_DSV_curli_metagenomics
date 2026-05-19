#!/bin/bash
#SBATCH --job-name=convert_rpk_to_rel_abund
#SBATCH --output=/scratch/project_2012191/chifang1/R_studio/curli/1206_relab/convert_rpk_to_rel_abund_%j.out
#SBATCH --error=/scratch/project_2012191/chifang1/R_studio/curli/1206_relab/convert_rpk_to_rel_abund_%j.err
#SBATCH --account=project_2012806
#SBATCH --partition=small
#SBATCH --time=3:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=20

############################################################
# Convert HUMAnN genefamilies RPK values to relative abundance
#
# Aim:
# Convert HUMAnN genefamilies output from RPK values to
# sample-wise relative abundance values before targeted
# extraction of E. coli curli-related gene families.
#
# Main steps:
# 1. Iterate over cleaned HUMAnN genefamilies tables.
# 2. Exclude unmapped and unintegrated features from the total.
# 3. Calculate the total mapped RPK value for each sample.
# 4. Convert each gene-family RPK value to relative abundance.
#
# Input:
# *_genefamilies_cleaned.tsv files
#
# Output:
# *_rel_abundance.tsv files
############################################################

set -euo pipefail

# ----------------------------------------------------------
# 1. Path configuration
# ----------------------------------------------------------

INPUT_DIR="/scratch/project_2012191/chifang1/R_studio/curli/genefamilies_cleaned"
OUTPUT_DIR="/scratch/project_2012191/chifang1/R_studio/curli/genefamilies_relative_abundance"

mkdir -p "$OUTPUT_DIR"

LOG_FILE="$OUTPUT_DIR/conversion_log.txt"

echo "Start RPK-to-relative-abundance conversion: $(date)" > "$LOG_FILE"
echo "Input directory: $INPUT_DIR" | tee -a "$LOG_FILE"
echo "Output directory: $OUTPUT_DIR" | tee -a "$LOG_FILE"

# ----------------------------------------------------------
# 2. Validate input directory
# ----------------------------------------------------------

if [ ! -d "$INPUT_DIR" ]; then
  echo "ERROR: Input directory not found: $INPUT_DIR" | tee -a "$LOG_FILE" >&2
  exit 1
fi

shopt -s nullglob
input_files=("$INPUT_DIR"/*_genefamilies_cleaned.tsv)

if [ "${#input_files[@]}" -eq 0 ]; then
  echo "ERROR: No *_genefamilies_cleaned.tsv files found in $INPUT_DIR" | tee -a "$LOG_FILE" >&2
  exit 1
fi

echo "Number of input files: ${#input_files[@]}" | tee -a "$LOG_FILE"

# ----------------------------------------------------------
# 3. Convert each sample
# ----------------------------------------------------------

for file in "${input_files[@]}"; do

  sample=$(basename "$file" "_genefamilies_cleaned.tsv")
  output_file="$OUTPUT_DIR/${sample}_rel_abundance.tsv"

  total_rpk=$(awk '
    NR > 1 &&
    $1 != "UNMAPPED" &&
    $1 != "UNINTEGRATED" &&
    $2 ~ /^[0-9.eE+-]+$/ {
      sum += $2
    }
    END {
      print sum
    }
  ' "$file")

  if [ -z "$total_rpk" ] || [ "$total_rpk" = "0" ]; then
    echo "Skipping $sample: total mapped RPK is zero or missing." | tee -a "$LOG_FILE"
    continue
  fi

  awk -v total="$total_rpk" '
    BEGIN {
      FS = "\t"
      OFS = "\t"
    }

    NR == 1 {
      print $1, "Relative_abundance"
      next
    }

    $1 != "UNMAPPED" &&
    $1 != "UNINTEGRATED" &&
    $2 ~ /^[0-9.eE+-]+$/ {
      print $1, $2 / total
    }
  ' "$file" > "$output_file"

  echo "Converted: $sample" | tee -a "$LOG_FILE"

done

echo "RPK-to-relative-abundance conversion completed: $(date)" | tee -a "$LOG_FILE"
