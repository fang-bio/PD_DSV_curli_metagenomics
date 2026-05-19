############################################################
# Figure 1B: Global microbiome structure by disease status
#
# Aim:
# Evaluate global species-level microbiome structure between
# Parkinson's disease patients and controls using Aitchison
# distance and PCoA ordination.
#
# Main analyses:
# 1. Load and align sample metadata with species-level abundance profiles.
# 2. Apply CLR transformation to species-level relative abundance data.
# 3. Compute Aitchison distance.
# 4. Perform PCoA ordination.
# 5. Test the disease-status effect using cohort-stratified PERMANOVA.
# 6. Generate Figure 1B and export source data.
#
# Computing environment:
# The PERMANOVA analysis used 9,999 permutations and was run on
# the Puhti high-performance computing environment because the
# Aitchison-distance analysis of species-level profiles was
# computationally intensive.
#
# Input:
# 1. metadata_924_685.xlsx
# 2. species_collapsed_fraction.xlsx
#
# Output:
# Figure 1B PCoA plot, PCoA coordinates, Aitchison distance matrix,
# PERMANOVA result tables, and source-data tables.
############################################################

suppressWarnings(suppressMessages({
  library(tidyverse)
  library(readxl)
  library(openxlsx)
  library(vegan)
  library(ape)
}))

set.seed(2025)

# ----------------------------------------------------------
# 1. Path configuration
# ----------------------------------------------------------

base_dir <- "/users/chifang1/CorrectMetada/3Figures_share"

meta_file <- file.path(base_dir, "metadata_924_685.xlsx")
species_file <- file.path(base_dir, "species_collapsed_fraction.xlsx")

outdir <- file.path(base_dir, "Figure1B_global_microbiome_structure")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(meta_file))
stopifnot(file.exists(species_file))

cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Metadata file:", meta_file, "\n")
cat("Species file :", species_file, "\n")
cat("Output dir   :", outdir, "\n")

n_perm <- 9999

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(n_cores) || n_cores < 1) {
  n_cores <- 1
}

cat("PERMANOVA permutations:", n_perm, "\n")
cat("PERMANOVA cores:", n_cores, "\n")

# ----------------------------------------------------------
# 2. Load and harmonize metadata
# ----------------------------------------------------------

cat("\nLoading metadata...\n")

meta <- read_xlsx(meta_file)

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

meta <- meta %>%
  mutate(
    SampleID = as.character(SampleID),
    Group = case_when(
      Group %in% c("PD", "Patient", "PD Patient") ~ "PD",
      Group %in% c("Control", "HC", "Healthy", "Population Control") ~ "Control",
      TRUE ~ as.character(Group)
    ),
    Group = factor(Group, levels = c("Control", "PD")),
    Cohort = factor(Cohort),
    Continent = if ("Continent" %in% colnames(.)) factor(Continent) else factor(NA),
    Platform = if ("Platform" %in% colnames(.)) factor(Platform) else factor(NA)
  ) %>%
  distinct(SampleID, .keep_all = TRUE)

if (anyNA(meta$Group)) {
  stop("Group contains NA values after label harmonization. Please check metadata.")
}

cat("Metadata samples:", nrow(meta), "\n")
cat("Group counts:\n")
print(table(meta$Group))
cat("Cohort count:", nlevels(meta$Cohort), "\n")

# ----------------------------------------------------------
# 3. Load and align species abundance table
# ----------------------------------------------------------

cat("\nLoading species-level abundance table...\n")

abun_raw <- read_xlsx(species_file)

if (!"Species" %in% colnames(abun_raw)) {
  colnames(abun_raw)[1] <- "Species"
}

sample_cols <- colnames(abun_raw)[colnames(abun_raw) != "Species"]

cat("Species:", nrow(abun_raw), "\n")
cat("Sample columns in species table:", length(sample_cols), "\n")

common_samples <- intersect(meta$SampleID, sample_cols)

if (length(common_samples) == 0) {
  stop("No matching sample IDs between metadata and species table.")
}

cat("Matched samples:", length(common_samples), "\n")

meta <- meta %>%
  filter(SampleID %in% common_samples) %>%
  arrange(SampleID) %>%
  as.data.frame()

rownames(meta) <- meta$SampleID

mat <- abun_raw %>%
  select(Species, all_of(meta$SampleID)) %>%
  column_to_rownames("Species") %>%
  as.matrix() %>%
  t()

mode(mat) <- "numeric"

if (anyNA(mat)) {
  stop("Species abundance matrix contains NA values after numeric conversion.")
}

if (!all(rownames(mat) == meta$SampleID)) {
  stop("Sample alignment mismatch between species matrix and metadata.")
}

cat("Aligned species matrix:", nrow(mat), "samples x", ncol(mat), "species\n")

# ----------------------------------------------------------
# 4. Relative abundance scaling, CLR transformation, and Aitchison distance
# ----------------------------------------------------------

cat("\nPreparing compositional data...\n")

row_sums <- rowSums(mat, na.rm = TRUE)

if (any(row_sums == 0)) {
  stop("At least one sample has zero total abundance.")
}

if (all(abs(row_sums - 1) < 0.01)) {
  cat("Species table detected as fractional relative abundance.\n")
  mat_rel <- mat
} else if (all(abs(row_sums - 100) < 1)) {
  cat("Species table detected as percentage abundance. Converting to fractions.\n")
  mat_rel <- mat / 100
} else {
  cat("Species table detected as non-normalized abundance. Converting to relative abundance.\n")
  mat_rel <- sweep(mat, 1, row_sums, "/")
}

zero_percentage <- 100 * mean(mat_rel == 0, na.rm = TRUE)

cat(
  "Relative abundance row-sum range:",
  paste(round(range(rowSums(mat_rel)), 4), collapse = " - "),
  "\n"
)
cat("Zero percentage:", round(zero_percentage, 2), "%\n")

pseudocount <- 1e-6

log_mat <- log(mat_rel + pseudocount)
mat_clr <- sweep(log_mat, 1, rowMeans(log_mat), "-")

dist_ait <- dist(mat_clr, method = "euclidean")

saveRDS(
  mat_rel,
  file.path(outdir, "figure1b_species_relative_abundance.rds")
)

saveRDS(
  mat_clr,
  file.path(outdir, "figure1b_species_clr_matrix.rds")
)

saveRDS(
  dist_ait,
  file.path(outdir, "figure1b_aitchison_distance.rds")
)

abundance_qc <- data.frame(
  Metric = c(
    "Samples analyzed",
    "Species analyzed",
    "Minimum row sum before normalization",
    "Maximum row sum before normalization",
    "Median row sum before normalization",
    "Zero percentage",
    "Pseudocount"
  ),
  Value = c(
    nrow(mat),
    ncol(mat),
    min(row_sums),
    max(row_sums),
    median(row_sums),
    zero_percentage,
    pseudocount
  )
)

write.xlsx(
  abundance_qc,
  file.path(outdir, "figure1b_abundance_qc_summary.xlsx"),
  rowNames = FALSE
)

# ----------------------------------------------------------
# 5. PCoA ordination
# ----------------------------------------------------------

cat("\nRunning PCoA...\n")

pcoa_res <- ape::pcoa(dist_ait)

var_exp <- round(100 * pcoa_res$values$Relative_eig[1:5], 4)

cat("PCoA variance explained:\n")
for (i in seq_along(var_exp)) {
  cat(sprintf("PCoA%d: %.4f%%\n", i, var_exp[i]))
}

pcoa_df <- as.data.frame(pcoa_res$vectors[, 1:5, drop = FALSE])
colnames(pcoa_df) <- paste0("PCoA", seq_len(ncol(pcoa_df)))
pcoa_df$SampleID <- rownames(pcoa_df)

pcoa_df <- pcoa_df %>%
  left_join(
    meta %>%
      rownames_to_column("SampleID_row") %>%
      select(-SampleID_row),
    by = "SampleID"
  )

group_centers <- pcoa_df %>%
  group_by(Group) %>%
  summarise(
    PCoA1 = mean(PCoA1, na.rm = TRUE),
    PCoA2 = mean(PCoA2, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

pcoa_variance <- data.frame(
  Axis = paste0("PCoA", seq_along(var_exp)),
  Variance_Explained_Percent = var_exp
)

# ----------------------------------------------------------
# 6. PERMANOVA analyses
# ----------------------------------------------------------

cat("\nRunning PERMANOVA analyses...\n")

adon_stratified <- vegan::adonis2(
  dist_ait ~ Group,
  data = meta,
  permutations = n_perm,
  strata = meta$Cohort,
  parallel = n_cores
)

adon_adjusted <- vegan::adonis2(
  dist_ait ~ Cohort + Group,
  data = meta,
  permutations = n_perm,
  by = "terms",
  parallel = n_cores
)

adon_cohort <- vegan::adonis2(
  dist_ait ~ Cohort,
  data = meta,
  permutations = n_perm,
  parallel = n_cores
)

stratified_df <- as.data.frame(adon_stratified)
stratified_df$Term <- rownames(stratified_df)

adjusted_df <- as.data.frame(adon_adjusted)
adjusted_df$Term <- rownames(adjusted_df)

cohort_df <- as.data.frame(adon_cohort)
cohort_df$Term <- rownames(cohort_df)

R2_stratified <- stratified_df$R2[1]
P_stratified <- stratified_df$`Pr(>F)`[1]

R2_adjusted_group <- adjusted_df %>%
  filter(Term == "Group") %>%
  pull(R2)

P_adjusted_group <- adjusted_df %>%
  filter(Term == "Group") %>%
  pull(`Pr(>F)`)

R2_cohort <- cohort_df$R2[1]
P_cohort <- cohort_df$`Pr(>F)`[1]

format_p <- function(p) {
  case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3g", p)
  )
}

permanova_summary <- tibble(
  Model = c(
    "Disease effect stratified by cohort",
    "Disease effect adjusted for cohort",
    "Cohort effect"
  ),
  Formula = c(
    "Aitchison distance ~ Group | strata(Cohort)",
    "Aitchison distance ~ Cohort + Group",
    "Aitchison distance ~ Cohort"
  ),
  Term = c(
    "Group",
    "Group",
    "Cohort"
  ),
  R2 = c(
    R2_stratified,
    R2_adjusted_group,
    R2_cohort
  ),
  R2_percent = 100 * R2,
  P_value = c(
    P_stratified,
    P_adjusted_group,
    P_cohort
  ),
  P_label = format_p(P_value)
)

write.xlsx(
  list(
    PERMANOVA_summary = permanova_summary,
    Stratified_PERMANOVA = stratified_df,
    Adjusted_PERMANOVA = adjusted_df,
    Cohort_PERMANOVA = cohort_df
  ),
  file.path(outdir, "figure1b_permanova_results.xlsx"),
  overwrite = TRUE
)

# ----------------------------------------------------------
# 7. Generate Figure 1B
# ----------------------------------------------------------

cat("\nGenerating Figure 1B plot...\n")

group_colors <- c(
  "Control" = "#4DAF4A",
  "PD" = "#F8766D"
)

annotation_label <- paste0(
  "Stratified PERMANOVA\n",
  "P = ", format_p(P_stratified), "\n",
  "R2 = ", sprintf("%.3f", R2_stratified)
)

p_figure1b <- ggplot(
  pcoa_df,
  aes(x = PCoA1, y = PCoA2, color = Group)
) +
  geom_point(
    size = 1.3,
    alpha = 0.75
  ) +
  stat_ellipse(
    aes(group = Group),
    level = 0.95,
    type = "norm",
    linetype = "dashed",
    linewidth = 0.55,
    alpha = 0.75,
    show.legend = FALSE
  ) +
  geom_point(
    data = group_centers,
    aes(x = PCoA1, y = PCoA2, color = Group),
    shape = 18,
    size = 4,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = annotation_label,
    hjust = 1.08,
    vjust = 1.25,
    size = 3.2
  ) +
  scale_color_manual(
    values = group_colors,
    name = "Group"
  ) +
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
  )

ggsave(
  filename = file.path(outdir, "figure1b_global_microbiome_pcoa.pdf"),
  plot = p_figure1b,
  width = 6.5,
  height = 5,
  device = "pdf"
)

ggsave(
  filename = file.path(outdir, "figure1b_global_microbiome_pcoa.tiff"),
  plot = p_figure1b,
  width = 6.5,
  height = 5,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(outdir, "figure1b_global_microbiome_pcoa.png"),
  plot = p_figure1b,
  width = 6.5,
  height = 5,
  dpi = 600,
  bg = "white"
)

# ----------------------------------------------------------
# 8. Export source data and analysis objects
# ----------------------------------------------------------

cat("\nExporting source data and analysis objects...\n")

write.xlsx(
  list(
    PCoA_coordinates = pcoa_df,
    Group_centers = group_centers,
    PCoA_variance = pcoa_variance,
    PERMANOVA_summary = permanova_summary,
    Abundance_QC = abundance_qc
  ),
  file.path(outdir, "figure1b_source_data_and_statistics.xlsx"),
  overwrite = TRUE
)

write.xlsx(
  pcoa_df,
  file.path(outdir, "figure1b_pcoa_coordinates.xlsx"),
  overwrite = TRUE
)

saveRDS(
  list(
    metadata = meta,
    species_relative_abundance = mat_rel,
    species_clr_matrix = mat_clr,
    aitchison_distance = dist_ait,
    pcoa = pcoa_res,
    pcoa_coordinates = pcoa_df,
    group_centers = group_centers,
    pcoa_variance = pcoa_variance,
    adon_stratified = adon_stratified,
    adon_adjusted = adon_adjusted,
    adon_cohort = adon_cohort,
    permanova_summary = permanova_summary,
    figure1b_plot = p_figure1b
  ),
  file.path(outdir, "figure1b_analysis_objects.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "sessionInfo.txt")
)

cat("\nFigure 1B global microbiome structure analysis completed successfully.\n")
cat("Outputs written to:", outdir, "\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
