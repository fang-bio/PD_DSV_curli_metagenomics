############################################################
# Figure 2A: Preprocessing of targeted microbial features
#
# Aim:
# Prepare sample-level metadata and targeted microbial features
# for Figure 2 analyses, including Desulfovibrio, Escherichia coli,
# and bacterial curli genes.
#
# Input:
# 1. metadata_924_685.xlsx
#    - required columns: SampleID, Group, Cohort
#
# 2. 1C_3Features_wide_correct.xlsx
#    - required columns: SampleID, Desulfovibrio_relab,
#      Ecoli_relab, Ecoli_curligene_relab
#
# Output:
# Preprocessed data table for downstream Figure 2 analyses.
############################################################

suppressWarnings(suppressMessages({
  library(tidyverse)
  library(readxl)
  library(writexl)
}))

set.seed(2025)

# ----------------------------------------------------------
# 1. Path configuration
# ----------------------------------------------------------

dir_base <- "C:/Users/chifang/OneDrive - University of Helsinki/Manuscripts/Tommi_Vatanen/3Figures_share"

file_meta <- file.path(dir_base, "metadata_924_685.xlsx")
file_feat <- file.path(dir_base, "1C_3Features_wide_correct.xlsx")

dir_out <- file.path(dir_base, "Figure2_Outputs", "Step0_Preprocessed")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(file_meta))
stopifnot(file.exists(file_feat))

cat("Metadata file:", file_meta, "\n")
cat("Feature file :", file_feat, "\n")

# ----------------------------------------------------------
# 2. Helper function
# ----------------------------------------------------------

clean_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}

get_half_min_pseudocount <- function(x) {
  nz <- x[x > 0 & is.finite(x)]
  if (length(nz) == 0) return(1e-12)
  min(nz) / 2
}

# ----------------------------------------------------------
# 3. Load metadata
# ----------------------------------------------------------

meta <- read_xlsx(file_meta)

required_meta <- c("SampleID", "Group", "Cohort")
missing_meta <- setdiff(required_meta, colnames(meta))

if (length(missing_meta) > 0) {
  stop(
    "Metadata missing required columns: ",
    paste(missing_meta, collapse = ", "),
    "\nDetected columns: ",
    paste(colnames(meta), collapse = ", ")
  )
}

meta_keep <- meta %>%
  mutate(
    SampleID = as.character(SampleID),
    Group = case_when(
      Group %in% c("PD", "Patient", "PD Patient") ~ "PD",
      Group %in% c("Control", "HC", "Healthy", "Population Control") ~ "Control",
      TRUE ~ as.character(Group)
    ),
    Group = factor(Group, levels = c("Control", "PD")),
    Cohort = factor(Cohort),
    Age = if ("Age" %in% colnames(.)) as.numeric(Age) else NA_real_,
    Sex = if ("Sex" %in% colnames(.)) as.factor(Sex) else factor(NA),
    BMI = if ("BMI" %in% colnames(.)) as.numeric(BMI) else NA_real_
  ) %>%
  select(SampleID, Group, Cohort, Age, Sex, BMI) %>%
  distinct(SampleID, .keep_all = TRUE)

if (any(is.na(meta_keep$Group))) {
  stop("Group contains NA after label normalization. Please check Group labels.")
}

# ----------------------------------------------------------
# 4. Load targeted feature table
# ----------------------------------------------------------

feat <- read_xlsx(file_feat)

required_feat <- c(
  "SampleID",
  "Desulfovibrio_relab",
  "Ecoli_relab",
  "Ecoli_curligene_relab"
)

missing_feat <- setdiff(required_feat, colnames(feat))

if (length(missing_feat) > 0) {
  stop(
    "Feature table missing required columns: ",
    paste(missing_feat, collapse = ", "),
    "\nDetected columns: ",
    paste(colnames(feat), collapse = ", ")
  )
}

feat_keep <- feat %>%
  mutate(
    SampleID  = as.character(SampleID),
    DSV_raw   = clean_numeric(Desulfovibrio_relab),
    Ecoli_raw = clean_numeric(Ecoli_relab),
    Curli_raw = clean_numeric(Ecoli_curligene_relab)
  ) %>%
  select(SampleID, DSV_raw, Ecoli_raw, Curli_raw)

if (anyNA(feat_keep$DSV_raw) ||
    anyNA(feat_keep$Ecoli_raw) ||
    anyNA(feat_keep$Curli_raw)) {
  stop("Feature table contains NA values after numeric conversion. Please check the input file.")
}

# ----------------------------------------------------------
# 5. Merge metadata and features
# ----------------------------------------------------------

df0 <- meta_keep %>%
  inner_join(feat_keep, by = "SampleID")

cat("Merged sample size:", nrow(df0), "\n")
cat("Group counts:\n")
print(table(df0$Group))
cat("Cohort counts:\n")
print(table(df0$Cohort))

# ----------------------------------------------------------
# 6. Log transformation and standardization
# ----------------------------------------------------------

pc_dsv   <- get_half_min_pseudocount(df0$DSV_raw)
pc_ecoli <- get_half_min_pseudocount(df0$Ecoli_raw)
pc_curli <- get_half_min_pseudocount(df0$Curli_raw)

cat("Pseudocount DSV  :", pc_dsv, "\n")
cat("Pseudocount Ecoli:", pc_ecoli, "\n")
cat("Pseudocount Curli:", pc_curli, "\n")

df0 <- df0 %>%
  mutate(
    DSV_detected   = DSV_raw > 0,
    Ecoli_detected = Ecoli_raw > 0,
    Curli_detected = Curli_raw > 0,

    DSV_pc   = DSV_raw + pc_dsv,
    Ecoli_pc = Ecoli_raw + pc_ecoli,
    Curli_pc = Curli_raw + pc_curli,

    log_DSV   = log10(DSV_pc),
    log_Ecoli = log10(Ecoli_pc),
    log_Curli = log10(Curli_pc),

    DSV_z   = as.numeric(scale(log_DSV)),
    Ecoli_z = as.numeric(scale(log_Ecoli)),
    Curli_z = as.numeric(scale(log_Curli))
  )

# ----------------------------------------------------------
# 7. Save outputs
# ----------------------------------------------------------

saveRDS(
  df0,
  file.path(dir_out, "df0_preprocessed.rds")
)

write_xlsx(
  df0,
  file.path(dir_out, "df0_preprocessed.xlsx")
)

summary_tbl <- df0 %>%
  summarise(
    n_total = n(),
    n_PD = sum(Group == "PD"),
    n_Control = sum(Group == "Control"),

    DSV_zero_rate = mean(DSV_raw == 0, na.rm = TRUE),
    Ecoli_zero_rate = mean(Ecoli_raw == 0, na.rm = TRUE),
    Curli_zero_rate = mean(Curli_raw == 0, na.rm = TRUE),

    DSV_detection_PD = mean(DSV_detected[Group == "PD"], na.rm = TRUE),
    DSV_detection_Control = mean(DSV_detected[Group == "Control"], na.rm = TRUE),

    Ecoli_detection_PD = mean(Ecoli_detected[Group == "PD"], na.rm = TRUE),
    Ecoli_detection_Control = mean(Ecoli_detected[Group == "Control"], na.rm = TRUE),

    Curli_detection_PD = mean(Curli_detected[Group == "PD"], na.rm = TRUE),
    Curli_detection_Control = mean(Curli_detected[Group == "Control"], na.rm = TRUE)
  )

write_xlsx(
  summary_tbl,
  file.path(dir_out, "Step0_summary.xlsx")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(dir_out, "sessionInfo.txt")
)

cat("Figure 2 preprocessing completed successfully.\n")
cat("Outputs written to:", dir_out, "\n")
