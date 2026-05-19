############################################################
# Figure 3C-D: Synergy-associated pathway heatmaps
#
# Aim:
# Generate group-averaged heatmaps for pathways positively and
# negatively associated with the Desulfovibrio-curli interaction.
#
# Input:
# 1. figure3c_positive_synergy_heatmap_matrix.xlsx
# 2. figure3d_negative_synergy_heatmap_matrix.xlsx
# 3. figure3c_positive_synergy_pathways.xlsx
# 4. figure3d_negative_synergy_pathways.xlsx
#
# Output:
# Figure 3C and Figure 3D heatmaps, together with source data and
# pathway statistics used for plotting.
############################################################

suppressWarnings(suppressMessages({
  library(tidyverse)
  library(readxl)
  library(pheatmap)
  library(RColorBrewer)
  library(writexl)
}))

set.seed(2025)

# ==========================================================
# 1. Path configuration
# ==========================================================

dir_data <- "C:/Users/chifang/OneDrive - University of Helsinki/Manuscripts/Tommi_Vatanen/3Figures_share/Figure3_Outputs/Fig3CD_heatmap_preprocessing"

dir_out <- file.path(
  "C:/Users/chifang/OneDrive - University of Helsinki/Manuscripts/Tommi_Vatanen/3Figures_share/Figure3_Outputs",
  "Fig3CD_heatmaps"
)

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

file_mat_pos <- file.path(dir_data, "figure3c_positive_synergy_heatmap_matrix.xlsx")
file_mat_neg <- file.path(dir_data, "figure3d_negative_synergy_heatmap_matrix.xlsx")
file_pos_paths <- file.path(dir_data, "figure3c_positive_synergy_pathways.xlsx")
file_neg_paths <- file.path(dir_data, "figure3d_negative_synergy_pathways.xlsx")

stopifnot(file.exists(file_mat_pos))
stopifnot(file.exists(file_mat_neg))
stopifnot(file.exists(file_pos_paths))
stopifnot(file.exists(file_neg_paths))

cat("Start Figure 3C-D heatmap plotting\n")
cat("Positive matrix:", file_mat_pos, "\n")
cat("Negative matrix:", file_mat_neg, "\n")
cat("Output dir:", dir_out, "\n")

# ==========================================================
# 2. Load heatmap matrices and pathway tables
# ==========================================================

mat_pos <- read_xlsx(file_mat_pos) %>%
  column_to_rownames("Pathway") %>%
  as.matrix()

mat_neg <- read_xlsx(file_mat_neg) %>%
  column_to_rownames("Pathway") %>%
  as.matrix()

mode(mat_pos) <- "numeric"
mode(mat_neg) <- "numeric"

if (anyNA(mat_pos)) {
  stop("Positive heatmap matrix contains NA values.")
}

if (anyNA(mat_neg)) {
  stop("Negative heatmap matrix contains NA values.")
}

pos_paths <- read_xlsx(file_pos_paths)
neg_paths <- read_xlsx(file_neg_paths)

cat("Positive heatmap matrix:", nrow(mat_pos), "pathways x", ncol(mat_pos), "groups\n")
cat("Negative heatmap matrix:", nrow(mat_neg), "pathways x", ncol(mat_neg), "groups\n")

# ==========================================================
# 3. Ensure consistent column order
# ==========================================================

group_order <- c(
  "Both Low",
  "Curli Only",
  "DSV Only",
  "Both High"
)

missing_pos_groups <- setdiff(group_order, colnames(mat_pos))
missing_neg_groups <- setdiff(group_order, colnames(mat_neg))

if (length(missing_pos_groups) > 0) {
  stop(
    "Positive matrix missing group columns: ",
    paste(missing_pos_groups, collapse = ", ")
  )
}

if (length(missing_neg_groups) > 0) {
  stop(
    "Negative matrix missing group columns: ",
    paste(missing_neg_groups, collapse = ", ")
  )
}

mat_pos <- mat_pos[, group_order, drop = FALSE]
mat_neg <- mat_neg[, group_order, drop = FALSE]

# ==========================================================
# 4. Heatmap color and annotations
# ==========================================================

heatmap_colors <- colorRampPalette(
  rev(brewer.pal(11, "RdBu"))
)(100)

ann_col <- data.frame(
  SynergyGroup = factor(
    group_order,
    levels = group_order
  )
)

rownames(ann_col) <- group_order

ann_colors <- list(
  SynergyGroup = c(
    "Both Low" = "#4DAF4A",
    "Curli Only" = "#FFB000",
    "DSV Only" = "#00BFC4",
    "Both High" = "#C77CFF"
  )
)

# ==========================================================
# 5. Plot Figure 3C: positive association heatmap
# ==========================================================

pheatmap(
  mat_pos,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 9,
  fontsize_col = 11,
  border_color = NA,
  color = heatmap_colors,
  filename = file.path(dir_out, "figure3c_positive_synergy_groupmean_heatmap.pdf"),
  width = 9,
  height = 8
)

pheatmap(
  mat_pos,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 9,
  fontsize_col = 11,
  border_color = NA,
  color = heatmap_colors,
  filename = file.path(dir_out, "figure3c_positive_synergy_groupmean_heatmap.png"),
  width = 9,
  height = 8
)

# ==========================================================
# 6. Plot Figure 3D: negative association heatmap
# ==========================================================

pheatmap(
  mat_neg,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 9,
  fontsize_col = 11,
  border_color = NA,
  color = heatmap_colors,
  filename = file.path(dir_out, "figure3d_negative_synergy_groupmean_heatmap.pdf"),
  width = 9,
  height = 8
)

pheatmap(
  mat_neg,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 9,
  fontsize_col = 11,
  border_color = NA,
  color = heatmap_colors,
  filename = file.path(dir_out, "figure3d_negative_synergy_groupmean_heatmap.png"),
  width = 9,
  height = 8
)

# ==========================================================
# 7. Export source data and plotting statistics
# ==========================================================

cat("Exporting Figure 3C-D source data and statistics...\n")

mat_pos_df <- as.data.frame(mat_pos) %>%
  rownames_to_column("Pathway")

mat_neg_df <- as.data.frame(mat_neg) %>%
  rownames_to_column("Pathway")

pos_table <- pos_paths %>%
  mutate(
    Direction = "Positive association",
    Rank_by_abs_beta = row_number()
  )

neg_table <- neg_paths %>%
  mutate(
    Direction = "Negative association",
    Rank_by_abs_beta = row_number()
  )

summary_tbl <- tibble(
  Metric = c(
    "Positive pathways shown in Figure 3C",
    "Negative pathways shown in Figure 3D",
    "Groups shown in heatmaps"
  ),
  Value = c(
    nrow(mat_pos),
    nrow(mat_neg),
    ncol(mat_pos)
  )
)

write_xlsx(
  list(
    Figure3C_positive_pathways = pos_table,
    Figure3D_negative_pathways = neg_table,
    Figure3C_groupmean_matrix = mat_pos_df,
    Figure3D_groupmean_matrix = mat_neg_df,
    Summary_statistics = summary_tbl
  ),
  file.path(dir_out, "figure3cd_heatmap_source_data_and_statistics.xlsx")
)

write.csv(
  mat_pos_df,
  file.path(dir_out, "figure3c_groupmean_heatmap_matrix.csv"),
  row.names = FALSE
)

write.csv(
  mat_neg_df,
  file.path(dir_out, "figure3d_groupmean_heatmap_matrix.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(dir_out, "sessionInfo.txt")
)

print(summary_tbl)

cat("Figure 3C-D heatmap plotting completed successfully.\n")
cat("Outputs written to:", dir_out, "\n")
