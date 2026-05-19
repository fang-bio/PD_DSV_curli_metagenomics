############################################################
# Figure 1C: Targeted MaAsLin3 analysis
#
# Aim:
# Test whether Desulfovibrio, Escherichia coli, and bacterial
# curli genes differ between Parkinson's disease patients and
# controls after accounting for cohort-level heterogeneity.
#
# Method:
# MaAsLin3 mixed-effects model
#
# Model:
# Feature ~ Group + (1 | Cohort)
#
# Input:
# 1. 1C_3Features_long_correct.xlsx
#    - rows = features
#    - columns = samples
#
# 2. metadata_924_685.xlsx
#    - rows = samples
#    - required columns: SampleID, Group, Cohort
#
# Output:
# MaAsLin3 result files for Figure 1C.
############################################################

suppressWarnings(suppressMessages({
  library(maaslin3)
  library(readxl)
  library(dplyr)
  library(tibble)
}))

set.seed(2025)

# ==========================================================
# 1. File paths
# ==========================================================

feature_file  <- "C:/Users/chifang/OneDrive - University of Helsinki/Manuscripts/Tommi_Vatanen/3Figures_share/1C_3Features_long_correct.xlsx"
metadata_file <- "C:/Users/chifang/OneDrive - University of Helsinki/Manuscripts/Tommi_Vatanen/3Figures_share/metadata_924_685.xlsx"

output_dir <- "C:/Users/chifang/OneDrive - University of Helsinki/Manuscripts/Tommi_Vatanen/3Figures_share/1C_MaAsLin3_results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================================
# 2. Load feature table
#
# Original format:
# rows = features, columns = samples
#
# MaAsLin3 required format:
# rows = samples, columns = features
# ==========================================================

feat_raw <- read_excel(feature_file)

# The first column contains feature names; the remaining columns are SampleIDs
names(feat_raw)[1] <- "feature"

feat_mat <- feat_raw %>%
  column_to_rownames("feature") %>%
  as.matrix()

# Transpose to Sample × Feature
input_data <- as.data.frame(t(feat_mat))

# Convert feature values to numeric
input_data[] <- lapply(input_data, as.numeric)

# Check for missing values introduced during numeric conversion
if (anyNA(input_data)) {
  stop("Feature table contains NA values after numeric conversion. Please check the input file.")
}

# ==========================================================
# 3. Load metadata
# ==========================================================

input_metadata <- read_excel(metadata_file) %>%
  mutate(
    SampleID = as.character(SampleID),
    Group    = factor(Group, levels = c("Control", "PD")),
    Cohort   = factor(Cohort)
  ) %>%
  distinct(SampleID, .keep_all = TRUE) %>%
  column_to_rownames("SampleID")

# ==========================================================
# 4. Align samples between feature table and metadata
# ==========================================================

common_samples <- intersect(rownames(input_data), rownames(input_metadata))

cat("Matched samples:", length(common_samples), "\n")

input_data     <- input_data[common_samples, , drop = FALSE]
input_metadata <- input_metadata[common_samples, , drop = FALSE]

stopifnot(identical(rownames(input_data), rownames(input_metadata)))

# ==========================================================
# 5. Sanity checks
# ==========================================================

cat("\nFeature names:\n")
print(colnames(input_data))

cat("\nFeature value range:\n")
print(summary(as.vector(as.matrix(input_data))))

cat("\nZero proportion:\n")
print(mean(input_data == 0, na.rm = TRUE))

cat("\nGroup counts:\n")
print(table(input_metadata$Group))

cat("\nCohort counts:\n")
print(table(input_metadata$Cohort))

# ==========================================================
# 6. Run MaAsLin3
# ==========================================================

fit_out <- maaslin3(
  input_data     = input_data,
  input_metadata = input_metadata,
  output         = output_dir,

  fixed_effects  = c("Group"),
  random_effects = c("Cohort"),

  normalization  = "NONE",   # feature values have already been preprocessed
  transform      = "LOG",
  correction     = "BH",

  max_significance = 1.0,    # retain complete output
  standardize      = FALSE,
  cores            = 1
)

cat("\nMaAsLin3 finished successfully.\n")
cat("Results saved to:", output_dir, "\n")

# Save session information for reproducibility
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)
