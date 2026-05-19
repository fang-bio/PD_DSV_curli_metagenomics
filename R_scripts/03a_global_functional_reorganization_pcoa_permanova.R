############################################################
# Figure 3A: Global functional reorganization
#
# Aim:
# Evaluate whether combined Desulfovibrio and curli-gene elevation
# is associated with global shifts in microbiome functional pathway
# profiles.
#
# Main analyses:
# 1. CLR transformation of pathway abundance profiles.
# 2. Aitchison distance calculation.
# 3. Principal coordinates analysis.
# 4. PERMANOVA to test cohort, Desulfovibrio, curli-gene, and
#    interaction effects.
# 5. Variance partitioning based on PERMANOVA R2 values.
#
# Computing environment:
# The PERMANOVA analysis uses 9,999 permutations and was run on
# the Puhti high-performance computing environment because the
# pathway-level Aitchison-distance analysis is computationally
# intensive.
#
# Input:
# 1. Fig3_Metadata_Final.csv
# 2. Fig3_Pathabundance_Final.tsv
#
# Output:
# Figure 3A PCoA plot, cohort-colored supplementary PCoA,
# variance-explained plot, and PERMANOVA result tables.
############################################################

suppressWarnings(suppressMessages({
  library(tidyverse)
  library(data.table)
  library(vegan)
  library(compositions)
  library(openxlsx)
  library(RColorBrewer)
}))

set.seed(2025)

# ==========================================================
# 1. Path configuration
# ==========================================================

dir_data <- "/users/chifang1/3Figures/Figure3_puhti"
dir_out  <- file.path(dir_data, "Fig3A_global_functional_reorganization")

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

meta_path <- file.path(dir_data, "Fig3_Metadata_Final.csv")
path_path <- file.path(dir_data, "Fig3_Pathabundance_Final.tsv")

stopifnot(file.exists(meta_path))
stopifnot(file.exists(path_path))

cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Metadata file:", meta_path, "\n")
cat("Pathway file :", path_path, "\n")
cat("Output dir   :", dir_out, "\n")

# Use available SLURM cores if running on Puhti
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(n_cores) || n_cores < 1) {
  n_cores <- 1
}

cat("PERMANOVA cores:", n_cores, "\n")

# ==========================================================
# 2. Load metadata
# ==========================================================

cat("\nLoading metadata...\n")

meta <- fread(meta_path, data.table = FALSE)

# Ensure the first column is SampleID if exported with altered column names
if (!"SampleID" %in% colnames(meta)) {
  colnames(meta)[1] <- "SampleID"
}

required_meta <- c(
  "SampleID",
  "Group",
  "Cohort",
  "SynergyGroup",
  "DSV_z",
  "Curli_z"
)

missing_meta <- setdiff(required_meta, colnames(meta))

if (length(missing_meta) > 0) {
  stop(
    "Metadata missing required columns: ",
    paste(missing_meta, collapse = ", "),
    "\nPlease check the Figure 3 preprocessing script."
  )
}

meta <- meta %>%
  mutate(
    SampleID = as.character(SampleID),
    Group = factor(Group, levels = c("Control", "PD")),
    Cohort = factor(Cohort),
    SynergyGroup = factor(
      SynergyGroup,
      levels = c("LL", "HL", "LH", "HH")
    ),
    DSV_z = as.numeric(DSV_z),
    Curli_z = as.numeric(Curli_z),
    DSVxCurli = DSV_z * Curli_z
  ) %>%
  distinct(SampleID, .keep_all = TRUE)

if (anyNA(meta$Group)) {
  stop("Group contains NA values. Please check group labels.")
}

if (anyNA(meta$SynergyGroup)) {
  stop("SynergyGroup contains NA values. Please check group definitions.")
}

cat("Metadata samples:", nrow(meta), "\n")
cat("Group counts:\n")
print(table(meta$Group))

cat("Synergy group counts:\n")
print(table(meta$SynergyGroup))

cat("Cohort counts:\n")
print(table(meta$Cohort))

# ==========================================================
# 3. Load pathway abundance matrix
# ==========================================================

cat("\nLoading pathway abundance matrix...\n")

pathway <- fread(path_path, data.table = FALSE)

if (!"Pathway" %in% colnames(pathway)) {
  colnames(pathway)[1] <- "Pathway"
}

pathway <- pathway %>%
  column_to_rownames("Pathway")

# Convert pathway table to numeric matrix
pathway_mat <- as.matrix(pathway)
mode(pathway_mat) <- "numeric"

if (anyNA(pathway_mat)) {
  stop("Pathway matrix contains NA values after numeric conversion.")
}

cat("Pathway matrix:", nrow(pathway_mat), "pathways x", ncol(pathway_mat), "samples\n")

# ==========================================================
# 4. Align pathway matrix and metadata
# ==========================================================

cat("\nAligning samples...\n")

common_samples <- sort(intersect(colnames(pathway_mat), meta$SampleID))

if (length(common_samples) == 0) {
  cat("Pathway sample IDs, first entries:\n")
  print(head(colnames(pathway_mat)))
  cat("Metadata sample IDs, first entries:\n")
  print(head(meta$SampleID))
  stop("No matching SampleIDs found between pathway matrix and metadata.")
}

meta_aligned <- meta %>%
  filter(SampleID %in% common_samples) %>%
  arrange(SampleID) %>%
  as.data.frame()

rownames(meta_aligned) <- meta_aligned$SampleID

pathway_aligned <- pathway_mat[, meta_aligned$SampleID, drop = FALSE]

if (!all(colnames(pathway_aligned) == rownames(meta_aligned))) {
  stop("Sample alignment mismatch between pathway matrix and metadata.")
}

cat("Aligned pathway matrix:",
    nrow(pathway_aligned), "pathways x",
    ncol(pathway_aligned), "samples\n")

cat("Aligned synergy group counts:\n")
print(table(meta_aligned$SynergyGroup))

# ==========================================================
# 5. CLR transformation and Aitchison distance
# ==========================================================

cat("\nPerforming CLR transformation and calculating Aitchison distance...\n")

pseudocount <- 1e-6

# Samples must be rows and pathways must be columns for CLR transformation
pathway_for_clr <- t(pathway_aligned + pseudocount)

pathway_clr <- compositions::clr(pathway_for_clr)
pathway_clr <- as.matrix(pathway_clr)

dist_aitch <- dist(pathway_clr, method = "euclidean")

# ==========================================================
# 6. PCoA analysis
# ==========================================================

cat("\nRunning PCoA...\n")

pcoa_res <- cmdscale(dist_aitch, k = 2, eig = TRUE)

pcoa_df <- as.data.frame(pcoa_res$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")

pcoa_df$SampleID <- rownames(pcoa_df)

pcoa_df <- pcoa_df %>%
  left_join(
    meta_aligned %>% rownames_to_column("SampleID_row"),
    by = "SampleID"
  ) %>%
  select(-SampleID_row)

positive_eigs <- pcoa_res$eig[pcoa_res$eig > 0]
eig_pct <- round(pcoa_res$eig[1:2] / sum(positive_eigs) * 100, 1)

cat("PCoA variance explained:\n")
cat("PCoA1:", eig_pct[1], "%\n")
cat("PCoA2:", eig_pct[2], "%\n")

# ==========================================================
# 7. PERMANOVA analysis
# ==========================================================

cat("\nRunning PERMANOVA with 9,999 permutations...\n")

n_perm <- 9999

adonis_formula <- dist_aitch ~ Cohort + DSV_z + Curli_z + DSV_z:Curli_z

adonis_terms <- adonis2(
  formula = adonis_formula,
  data = meta_aligned,
  permutations = n_perm,
  by = "terms",
  parallel = n_cores
)

adonis_marginal <- adonis2(
  formula = adonis_formula,
  data = meta_aligned,
  permutations = n_perm,
  by = "margin",
  parallel = n_cores
)

adonis_terms_df <- as.data.frame(adonis_terms)
adonis_terms_df$Term <- rownames(adonis_terms_df)

adonis_marginal_df <- as.data.frame(adonis_marginal)
adonis_marginal_df$Term <- rownames(adonis_marginal_df)

interaction_row <- adonis_terms_df %>%
  filter(Term == "DSV_z:Curli_z")

if (nrow(interaction_row) == 0) {
  interaction_p <- NA_real_
  interaction_r2 <- NA_real_
} else {
  interaction_p <- interaction_row$`Pr(>F)`[1]
  interaction_r2 <- interaction_row$R2[1]
}

p_str <- case_when(
  is.na(interaction_p) ~ "P = NA",
  interaction_p < 0.001 ~ "P < 0.001",
  TRUE ~ paste0("P = ", sprintf("%.3f", interaction_p))
)

r2_str <- case_when(
  is.na(interaction_r2) ~ "R2 = NA",
  TRUE ~ paste0("R2 = ", sprintf("%.2f", interaction_r2 * 100), "%")
)

cat("Interaction term:", r2_str, ",", p_str, "\n")

# ==========================================================
# 8. Variance partitioning table
# ==========================================================

cat("\nPreparing variance partitioning table...\n")

var_df <- adonis_terms_df %>%
  filter(!Term %in% c("Residual", "Total")) %>%
  mutate(
    R2_Pct = R2 * 100,
    Term_Label = case_when(
      Term == "Cohort" ~ "Cohort effect",
      Term == "DSV_z" ~ "Desulfovibrio abundance",
      Term == "Curli_z" ~ "Curli-gene abundance",
      Term == "DSV_z:Curli_z" ~ "DSV x Curli interaction",
      TRUE ~ Term
    ),
    P_label = case_when(
      `Pr(>F)` < 0.001 ~ "< 0.001",
      TRUE ~ sprintf("%.3f", `Pr(>F)`)
    )
  )

print(var_df[, c("Term_Label", "R2_Pct", "P_label")])

# ==========================================================
# 9. Figure 3A: PCoA colored by DSV-curli group
# ==========================================================

cat("\nGenerating Figure 3A PCoA plot...\n")

synergy_colors <- c(
  "LL" = "grey70",
  "HL" = "#4575B4",
  "LH" = "#E69F00",
  "HH" = "#D73027"
)

p_pcoa <- ggplot(
  pcoa_df,
  aes(x = PCoA1, y = PCoA2, color = SynergyGroup)
) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_color_manual(
    values = synergy_colors,
    name = "DSV-curli group",
    drop = FALSE
  ) +
  labs(
    title = "Global functional pathway profiles",
    subtitle = paste0(
      "DSV x Curli interaction: ",
      r2_str,
      ", ",
      p_str,
      "; PERMANOVA with cohort adjustment"
    ),
    x = paste0("PCoA1 (", eig_pct[1], "%)"),
    y = paste0("PCoA2 (", eig_pct[2], "%)")
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(dir_out, "figure3a_global_functional_pcoa.pdf"),
  plot = p_pcoa,
  width = 6.5,
  height = 5,
  device = "pdf"
)

ggsave(
  filename = file.path(dir_out, "figure3a_global_functional_pcoa.tiff"),
  plot = p_pcoa,
  width = 6.5,
  height = 5,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# ==========================================================
# 10. Variance-explained plot
# ==========================================================

cat("\nGenerating variance-explained plot...\n")

variance_colors <- c(
  "Cohort" = "#999999",
  "DSV_z" = "#66C2A5",
  "Curli_z" = "#FC8D62",
  "DSV_z:Curli_z" = "#D53E4F"
)

p_variance <- ggplot(
  var_df,
  aes(
    x = reorder(Term_Label, R2_Pct),
    y = R2_Pct,
    fill = Term
  )
) +
  geom_col(width = 0.7, color = "black", alpha = 0.85) +
  geom_text(
    aes(label = paste0(sprintf("%.2f", R2_Pct), "%")),
    hjust = -0.15,
    fontface = "bold",
    size = 3.8
  ) +
  scale_fill_manual(values = variance_colors) +
  scale_y_continuous(
    limits = c(0, max(var_df$R2_Pct, na.rm = TRUE) * 1.25)
  ) +
  coord_flip() +
  labs(
    title = "Variance sources in functional pathway profiles",
    subtitle = "Sequential PERMANOVA R2 values",
    x = NULL,
    y = "Variance explained (%)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(dir_out, "figure3a_variance_explained.pdf"),
  plot = p_variance,
  width = 6,
  height = 4,
  device = "pdf"
)

ggsave(
  filename = file.path(dir_out, "figure3a_variance_explained.tiff"),
  plot = p_variance,
  width = 6,
  height = 4,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# ==========================================================
# 11. Supplementary cohort-colored PCoA
# ==========================================================

cat("\nGenerating cohort-colored supplementary PCoA plot...\n")

n_cohorts <- length(unique(pcoa_df$Cohort))

if (n_cohorts <= 12) {
  cohort_colors <- brewer.pal(max(3, n_cohorts), "Set3")[seq_len(n_cohorts)]
} else {
  cohort_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_cohorts)
}

names(cohort_colors) <- levels(factor(pcoa_df$Cohort))

p_cohort <- ggplot(
  pcoa_df,
  aes(x = PCoA1, y = PCoA2, color = Cohort)
) +
  geom_point(size = 1.4, alpha = 0.70) +
  scale_color_manual(values = cohort_colors) +
  labs(
    title = "Cohort structure in functional pathway profiles",
    subtitle = "PCoA based on Aitchison distance",
    x = paste0("PCoA1 (", eig_pct[1], "%)"),
    y = paste0("PCoA2 (", eig_pct[2], "%)"),
    color = "Cohort"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 8),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(dir_out, "suppfig3a_cohort_colored_functional_pcoa.pdf"),
  plot = p_cohort,
  width = 9,
  height = 6.5,
  device = "pdf"
)

ggsave(
  filename = file.path(dir_out, "suppfig3a_cohort_colored_functional_pcoa.tiff"),
  plot = p_cohort,
  width = 9,
  height = 6.5,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# ==========================================================
# 12. Export result tables
# ==========================================================

cat("\nExporting result tables...\n")

wb <- createWorkbook()

addWorksheet(wb, "PCoA_coordinates")
writeData(wb, "PCoA_coordinates", pcoa_df)

addWorksheet(wb, "PERMANOVA_terms")
writeData(wb, "PERMANOVA_terms", adonis_terms_df)

addWorksheet(wb, "PERMANOVA_marginal")
writeData(wb, "PERMANOVA_marginal", adonis_marginal_df)

addWorksheet(wb, "Variance_R2")
writeData(wb, "Variance_R2", var_df)

saveWorkbook(
  wb,
  file.path(dir_out, "figure3a_pcoa_permanova_outputs.xlsx"),
  overwrite = TRUE
)

saveRDS(
  list(
    metadata = meta_aligned,
    pathway_matrix = pathway_aligned,
    pathway_clr = pathway_clr,
    distance = dist_aitch,
    pcoa = pcoa_res,
    pcoa_coordinates = pcoa_df,
    adonis_terms = adonis_terms,
    adonis_marginal = adonis_marginal,
    variance_table = var_df,
    p_pcoa = p_pcoa,
    p_variance = p_variance,
    p_cohort = p_cohort
  ),
  file.path(dir_out, "figure3a_analysis_objects.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(dir_out, "sessionInfo.txt")
)

cat("\nFigure 3A global functional reorganization analysis completed successfully.\n")
cat("Outputs written to:", dir_out, "\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
