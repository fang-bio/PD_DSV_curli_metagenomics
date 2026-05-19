############################################################
# Figure 3C-D: Preparation of synergy-associated pathway heatmaps
#
# Aim:
# Prepare group-averaged pathway abundance matrices for the
# Figure 3C-D heatmaps based on pathway-level association results
# from the Desulfovibrio-curli interaction analysis.
#
# Main steps:
# 1. Load metadata, pathway abundance profiles, and MaAsLin3 results.
# 2. Select positively associated pathways for Figure 3C.
# 3. Select negatively associated pathways for Figure 3D.
# 4. Calculate log-transformed row-scaled pathway abundance.
# 5. Average pathway abundance across DSV-curli combined exposure groups.
#
# Input:
# 1. Fig3_Metadata_Final.xlsx
# 2. Fig3_Pathabundance_Unstratified_588.xlsx
# 3. Fig3B_Final_DSVxCurli_abundance_all.csv
#
# Output:
# Pathway lists and group-averaged heatmap matrices for Figure 3C-D.
############################################################

suppressWarnings(suppressMessages({
  library(tidyverse)
  library(readxl)
  library(data.table)
  library(writexl)
}))

set.seed(2025)

# ==========================================================
# 1. Path configuration
# ==========================================================

dir_data <- "C:/Users/chifang/OneDrive - University of Helsinki/Manuscripts/Tommi_Vatanen/3Figures_share/Figure3_Outputs"
dir_out  <- file.path(dir_data, "Fig3CD_heatmap_preprocessing")

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

file_meta <- file.path(dir_data, "Fig3_Metadata_Final.xlsx")
file_pathway <- file.path(dir_data, "Fig3_Pathabundance_Unstratified_588.xlsx")
file_results <- file.path(dir_data, "Fig3B_Final_DSVxCurli_abundance_all.csv")

stopifnot(file.exists(file_meta))
stopifnot(file.exists(file_pathway))
stopifnot(file.exists(file_results))

cat("Start Figure 3C-D heatmap preprocessing\n")
cat("Metadata file:", file_meta, "\n")
cat("Pathway file :", file_pathway, "\n")
cat("Result file  :", file_results, "\n")
cat("Output dir   :", dir_out, "\n")

# ==========================================================
# 2. Load metadata
# ==========================================================

meta <- read_xlsx(file_meta)

if (!"SampleID" %in% colnames(meta)) {
  colnames(meta)[1] <- "SampleID"
}

required_meta <- c("SampleID", "SynergyGroup")
missing_meta <- setdiff(required_meta, colnames(meta))

if (length(missing_meta) > 0) {
  stop(
    "Metadata missing required columns: ",
    paste(missing_meta, collapse = ", "),
    "\nDetected columns: ",
    paste(colnames(meta), collapse = ", ")
  )
}

meta <- meta %>%
  mutate(
    SampleID = as.character(SampleID),
    SynergyGroup = as.character(SynergyGroup)
  )

# Convert compact group labels to publication-friendly labels if needed
meta <- meta %>%
  mutate(
    SynergyGroup = recode(
      SynergyGroup,
      "LL" = "Both Low",
      "LH" = "Curli Only",
      "HL" = "DSV Only",
      "HH" = "Both High",
      .default = SynergyGroup
    ),
    SynergyGroup = factor(
      SynergyGroup,
      levels = c(
        "Both Low",
        "Curli Only",
        "DSV Only",
        "Both High"
      )
    )
  ) %>%
  distinct(SampleID, .keep_all = TRUE) %>%
  column_to_rownames("SampleID")

if (anyNA(meta$SynergyGroup)) {
  stop("SynergyGroup contains NA values after relabeling. Please check group labels.")
}

cat("Synergy group counts:\n")
print(table(meta$SynergyGroup))

# ==========================================================
# 3. Load pathway abundance matrix
# ==========================================================

path_raw <- read_xlsx(file_pathway)

if (!"Pathway" %in% colnames(path_raw)) {
  colnames(path_raw)[1] <- "Pathway"
}

path_mat <- path_raw %>%
  column_to_rownames("Pathway") %>%
  as.matrix()

mode(path_mat) <- "numeric"

if (anyNA(path_mat)) {
  stop("Pathway abundance matrix contains NA values after numeric conversion.")
}

cat("Pathway matrix:", nrow(path_mat), "pathways x", ncol(path_mat), "samples\n")

# ==========================================================
# 4. Load MaAsLin3 pathway association results
# ==========================================================

res_3b <- fread(file_results, data.table = FALSE) %>%
  as_tibble()

required_res <- c("feature", "coef")
missing_res <- setdiff(required_res, colnames(res_3b))

if (length(missing_res) > 0) {
  stop(
    "Result table missing required columns: ",
    paste(missing_res, collapse = ", "),
    "\nDetected columns: ",
    paste(colnames(res_3b), collapse = ", ")
  )
}

# Identify q-value column
if ("qval_joint" %in% colnames(res_3b)) {
  res_3b$target_q <- res_3b$qval_joint
} else if ("qval" %in% colnames(res_3b)) {
  res_3b$target_q <- res_3b$qval
} else if ("qval_individual" %in% colnames(res_3b)) {
  res_3b$target_q <- res_3b$qval_individual
} else {
  stop("No q-value column found. Expected one of: qval_joint, qval, qval_individual.")
}

res_3b <- res_3b %>%
  mutate(
    coef = as.numeric(coef),
    target_q = as.numeric(target_q)
  ) %>%
  filter(
    !is.na(coef),
    !is.na(target_q)
  )

cat("Loaded pathway association results:", nrow(res_3b), "\n")

# ==========================================================
# 5. Parameters
# ==========================================================

q_threshold <- 0.05
top_n <- 20

cat("Q-value threshold:", q_threshold, "\n")
cat("Top pathways per direction:", top_n, "\n")

# ==========================================================
# 6. Select positive and negative synergy-associated pathways
# ==========================================================

positive_paths <- res_3b %>%
  filter(
    target_q < q_threshold,
    coef > 0
  ) %>%
  arrange(desc(abs(coef))) %>%
  slice_head(n = top_n) %>%
  transmute(
    Pathway = feature,
    Beta_Coefficient = coef,
    Abs_Effect_Size = abs(coef),
    Q_value = target_q,
    Direction = "Positive association"
  )

negative_paths <- res_3b %>%
  filter(
    target_q < q_threshold,
    coef < 0
  ) %>%
  arrange(desc(abs(coef))) %>%
  slice_head(n = top_n) %>%
  transmute(
    Pathway = feature,
    Beta_Coefficient = coef,
    Abs_Effect_Size = abs(coef),
    Q_value = target_q,
    Direction = "Negative association"
  )

cat("Positive pathways selected:", nrow(positive_paths), "\n")
cat("Negative pathways selected:", nrow(negative_paths), "\n")

if (nrow(positive_paths) == 0) {
  warning("No positive pathways passed the selected q-value threshold.")
}

if (nrow(negative_paths) == 0) {
  warning("No negative pathways passed the selected q-value threshold.")
}

write_xlsx(
  positive_paths,
  file.path(dir_out, "figure3c_positive_synergy_pathways.xlsx")
)

write_xlsx(
  negative_paths,
  file.path(dir_out, "figure3d_negative_synergy_pathways.xlsx")
)

# ==========================================================
# 7. Build group-averaged heatmap matrices
# ==========================================================

get_half_min_pseudocount <- function(x) {
  nz <- x[x > 0 & is.finite(x)]

  if (length(nz) == 0) {
    return(1e-12)
  }

  min(nz) / 2
}

build_group_mean_matrix <- function(pathway_df, pathway_matrix, metadata) {

  if (nrow(pathway_df) == 0) {
    return(matrix(nrow = 0, ncol = 0))
  }

  missing_pathways <- setdiff(pathway_df$Pathway, rownames(pathway_matrix))

  if (length(missing_pathways) > 0) {
    warning(
      "The following selected pathways were not found in the pathway matrix: ",
      paste(missing_pathways, collapse = "; ")
    )
  }

  available_pathways <- intersect(pathway_df$Pathway, rownames(pathway_matrix))

  mat <- pathway_matrix[available_pathways, , drop = FALSE]

  common_samples <- intersect(colnames(mat), rownames(metadata))

  if (length(common_samples) == 0) {
    stop("No matching samples between pathway matrix and metadata.")
  }

  mat <- mat[, common_samples, drop = FALSE]
  meta_sub <- metadata[common_samples, , drop = FALSE]

  pseudocount <- get_half_min_pseudocount(as.vector(mat))

  mat_log <- log10(mat + pseudocount)

  mat_z <- t(scale(t(mat_log)))
  mat_z[is.na(mat_z)] <- 0

  mat_group <- mat_z %>%
    as.data.frame() %>%
    rownames_to_column("Pathway") %>%
    pivot_longer(
      cols = -Pathway,
      names_to = "SampleID",
      values_to = "Z"
    ) %>%
    left_join(
      meta_sub %>%
        rownames_to_column("SampleID") %>%
        select(SampleID, SynergyGroup),
      by = "SampleID"
    ) %>%
    filter(!is.na(SynergyGroup)) %>%
    group_by(Pathway, SynergyGroup) %>%
    summarise(
      MeanZ = mean(Z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = SynergyGroup,
      values_from = MeanZ
    ) %>%
    column_to_rownames("Pathway") %>%
    as.matrix()

  group_order <- c(
    "Both Low",
    "Curli Only",
    "DSV Only",
    "Both High"
  )

  missing_groups <- setdiff(group_order, colnames(mat_group))

  if (length(missing_groups) > 0) {
    for (grp in missing_groups) {
      mat_group <- cbind(mat_group, setNames(data.frame(NA_real_), grp))
    }
  }

  mat_group <- mat_group[
    available_pathways,
    group_order,
    drop = FALSE
  ]

  mat_group
}

mat_positive <- build_group_mean_matrix(
  pathway_df = positive_paths,
  pathway_matrix = path_mat,
  metadata = meta
)

mat_negative <- build_group_mean_matrix(
  pathway_df = negative_paths,
  pathway_matrix = path_mat,
  metadata = meta
)

# ==========================================================
# 8. Export heatmap matrices and source tables
# ==========================================================

write_xlsx(
  as.data.frame(mat_positive) %>%
    rownames_to_column("Pathway"),
  file.path(dir_out, "figure3c_positive_synergy_heatmap_matrix.xlsx")
)

write_xlsx(
  as.data.frame(mat_negative) %>%
    rownames_to_column("Pathway"),
  file.path(dir_out, "figure3d_negative_synergy_heatmap_matrix.xlsx")
)

write.csv(
  as.data.frame(mat_positive) %>%
    rownames_to_column("Pathway"),
  file.path(dir_out, "figure3c_positive_synergy_heatmap_matrix.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(mat_negative) %>%
    rownames_to_column("Pathway"),
  file.path(dir_out, "figure3d_negative_synergy_heatmap_matrix.csv"),
  row.names = FALSE
)

summary_tbl <- tibble(
  Metric = c(
    "Q-value threshold",
    "Top pathways per direction",
    "Positive pathways selected",
    "Negative pathways selected",
    "Samples with metadata",
    "Pathways in abundance matrix"
  ),
  Value = c(
    q_threshold,
    top_n,
    nrow(positive_paths),
    nrow(negative_paths),
    nrow(meta),
    nrow(path_mat)
  )
)

write_xlsx(
  summary_tbl,
  file.path(dir_out, "figure3cd_heatmap_preprocessing_summary.xlsx")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(dir_out, "sessionInfo.txt")
)

cat("Figure 3C-D heatmap preprocessing completed successfully.\n")
cat("Outputs written to:", dir_out, "\n")
