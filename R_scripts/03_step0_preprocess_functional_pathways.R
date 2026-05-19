############################################################
# Figure 3A: Preprocessing of functional pathway profiles
#
# Aim:
# Prepare pathway-level functional profiles for Figure 3 analyses
# by integrating sample metadata, Desulfovibrio and curli-gene
# abundance, and HUMAnN pathway abundance profiles.
#
# Main steps:
# 1. Load sample metadata and targeted microbial features.
# 2. Define DSV-curli combined exposure groups within each cohort.
# 3. Retain unstratified pathway profiles and remove unmapped or
#    unintegrated features.
# 4. Align pathway profiles with sample metadata.
# 5. Apply a pathway prevalence filter.
#
# Output:
# Preprocessed metadata and pathway abundance tables for downstream
# functional profiling analyses in Figure 3.
############################################################

suppressWarnings(suppressMessages({
  library(tidyverse)
  library(readxl)
  library(data.table)
  library(openxlsx)
}))

set.seed(2025)

# ==========================================================
# 1. Path configuration
# ==========================================================

dir_base <- "C:/Users/chifang/OneDrive - University of Helsinki/Manuscripts/Tommi_Vatanen/3Figures_share/Figure3_Outputs"
dir_out  <- "C:/Users/chifang/OneDrive - University of Helsinki/Manuscripts/Tommi_Vatanen/3Figures_share/Figure3_Outputs"

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

file_meta    <- "metadata_924_685.xlsx"
file_feat    <- "1C_3Features_wide_correct.xlsx"
file_pathway <- "1_pathabundance_4316_6030_relab.tsv"

file_meta_path    <- file.path(dir_base, file_meta)
file_feat_path    <- file.path(dir_base, file_feat)
file_pathway_path <- file.path(dir_base, file_pathway)

stopifnot(file.exists(file_meta_path))
stopifnot(file.exists(file_feat_path))
stopifnot(file.exists(file_pathway_path))

cat("Start time:", format(Sys.time(), "%H:%M:%S"), "\n")
cat("Metadata file:", file_meta_path, "\n")
cat("Feature file:", file_feat_path, "\n")
cat("Pathway file:", file_pathway_path, "\n")

# ==========================================================
# 2. Helper function
# ==========================================================

clean_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}

# ==========================================================
# 3. Load metadata and targeted microbial features
# ==========================================================

cat("Preparing metadata and targeted microbial features...\n")

meta <- read_xlsx(file_meta_path) %>%
  select(SampleID, Group, Cohort) %>%
  mutate(
    SampleID = as.character(SampleID),
    Group = factor(Group, levels = c("Control", "PD")),
    Cohort = factor(Cohort)
  ) %>%
  distinct(SampleID, .keep_all = TRUE)

if (anyNA(meta$Group)) {
  stop("Group contains NA values after factor conversion. Please check group labels.")
}

feat_raw <- read_xlsx(
  file_feat_path,
  col_types = "text"
)

required_feat <- c(
  "SampleID",
  "Desulfovibrio_relab",
  "Ecoli_curligene_relab"
)

missing_feat <- setdiff(required_feat, colnames(feat_raw))

if (length(missing_feat) > 0) {
  stop(
    "Feature table missing required columns: ",
    paste(missing_feat, collapse = ", "),
    "\nDetected columns: ",
    paste(colnames(feat_raw), collapse = ", ")
  )
}

feat_clean <- feat_raw %>%
  filter(SampleID != "SampleID") %>%
  transmute(
    SampleID = as.character(SampleID),
    DSV_relative_abundance = clean_numeric(Desulfovibrio_relab),
    Curli_relative_abundance = clean_numeric(Ecoli_curligene_relab)
  ) %>%
  filter(
    !is.na(DSV_relative_abundance),
    !is.na(Curli_relative_abundance)
  )

meta_feat <- meta %>%
  inner_join(feat_clean, by = "SampleID")

cat("Merged metadata and feature samples:", nrow(meta_feat), "\n")
cat("Group counts:\n")
print(table(meta_feat$Group))

cat("Cohort counts:\n")
print(table(meta_feat$Cohort))

# ==========================================================
# 4. Define DSV-curli combined exposure groups
# ==========================================================

cat("Defining DSV-curli combined exposure groups...\n")

meta_processed <- meta_feat %>%
  group_by(Cohort) %>%
  mutate(
    DSV_log = log1p(DSV_relative_abundance),
    Curli_log = log1p(Curli_relative_abundance),

    DSV_z = as.numeric(scale(DSV_log)),
    Curli_z = as.numeric(scale(Curli_log)),
    DSVxCurli = DSV_z * Curli_z,

    DSV_median = median(DSV_relative_abundance, na.rm = TRUE),
    Curli_median = median(Curli_relative_abundance, na.rm = TRUE),

    DSV_status = case_when(
      DSV_median == 0 & DSV_relative_abundance > 0 ~ "High",
      DSV_median == 0 & DSV_relative_abundance == 0 ~ "Low",
      DSV_relative_abundance >= DSV_median ~ "High",
      TRUE ~ "Low"
    ),

    Curli_status = case_when(
      Curli_median == 0 & Curli_relative_abundance > 0 ~ "High",
      Curli_median == 0 & Curli_relative_abundance == 0 ~ "Low",
      Curli_relative_abundance >= Curli_median ~ "High",
      TRUE ~ "Low"
    ),

    SynergyGroup = factor(
      paste0(
        ifelse(DSV_status == "High", "H", "L"),
        ifelse(Curli_status == "High", "H", "L")
      ),
      levels = c("LL", "HL", "LH", "HH")
    )
  ) %>%
  ungroup()

cat("Synergy group distribution before pathway alignment:\n")
print(table(meta_processed$SynergyGroup))

# ==========================================================
# 5. Load and clean pathway abundance table
# ==========================================================

cat("Loading and cleaning pathway abundance table...\n")

path_raw <- fread(
  file_pathway_path,
  data.table = FALSE
)

colnames(path_raw)[1] <- "Pathway"

path_clean <- path_raw %>%
  filter(!grepl("\\|", Pathway)) %>%
  filter(!grepl("UNMAPPED|UNINTEGRATED", Pathway)) %>%
  column_to_rownames("Pathway")

cat("Unstratified pathways before prevalence filtering:", nrow(path_clean), "\n")

path_unstratified_export <- path_clean %>%
  as.data.frame() %>%
  rownames_to_column("Pathway")

write.xlsx(
  path_unstratified_export,
  file.path(dir_out, "figure3_unstratified_pathways_before_prevalence_filter.xlsx")
)

fwrite(
  path_unstratified_export,
  file.path(dir_out, "figure3_unstratified_pathways_before_prevalence_filter.tsv"),
  sep = "\t"
)

# ==========================================================
# 6. Align pathway profiles with metadata
# ==========================================================

cat("Aligning pathway profiles with metadata...\n")

common_samples <- intersect(
  colnames(path_clean),
  meta_processed$SampleID
)

if (length(common_samples) == 0) {
  stop("No matching SampleIDs found between pathway table and metadata.")
}

meta_final <- meta_processed %>%
  filter(SampleID %in% common_samples) %>%
  arrange(SampleID) %>%
  as.data.frame()

rownames(meta_final) <- meta_final$SampleID

path_final <- path_clean[, meta_final$SampleID, drop = FALSE]

if (!all(colnames(path_final) == rownames(meta_final))) {
  stop("Sample alignment mismatch between pathway matrix and metadata.")
}

cat("Aligned samples:", ncol(path_final), "\n")

# ==========================================================
# 7. Apply pathway prevalence filter
# ==========================================================

cat("Applying pathway prevalence filter...\n")

prevalence <- rowSums(path_final > 0, na.rm = TRUE) / ncol(path_final)

prevalence_table <- tibble(
  Pathway = rownames(path_final),
  Prevalence = prevalence
)

path_final <- path_final[prevalence >= 0.10, , drop = FALSE]

cat("Pathways retained after prevalence >= 10% filter:", nrow(path_final), "\n")
cat("Final synergy group distribution after alignment:\n")
print(table(meta_final$SynergyGroup))

write.xlsx(
  prevalence_table,
  file.path(dir_out, "figure3_pathway_prevalence_summary.xlsx")
)

# ==========================================================
# 8. Export preprocessed data
# ==========================================================

cat("Exporting preprocessed Figure 3 input files...\n")

write.xlsx(
  meta_final,
  file.path(dir_out, "figure3_metadata_final.xlsx")
)

fwrite(
  meta_final,
  file.path(dir_out, "figure3_metadata_final.csv")
)

path_export <- path_final %>%
  rownames_to_column("Pathway")

write.xlsx(
  path_export,
  file.path(dir_out, "figure3_pathabundance_final.xlsx")
)

fwrite(
  path_export,
  file.path(dir_out, "figure3_pathabundance_final.tsv"),
  sep = "\t"
)

saveRDS(
  list(
    metadata = meta_final,
    pathway_abundance = path_final,
    pathway_prevalence = prevalence_table
  ),
  file.path(dir_out, "figure3_preprocessed_objects.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(dir_out, "sessionInfo.txt")
)

cat("Figure 3 functional pathway preprocessing completed successfully.\n")
cat("Outputs written to:", dir_out, "\n")
cat("End time:", format(Sys.time(), "%H:%M:%S"), "\n")
