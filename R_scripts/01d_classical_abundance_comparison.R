############################################################
# Figure 1D: Abundance validation using violin plots
#
# Aim:
# Compare the relative abundance distributions of Desulfovibrio,
# Escherichia coli, and bacterial curli genes between Parkinson's
# disease patients and controls using classical group-wise tests.
#
# Method:
# Wilcoxon rank-sum test
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
# Violin plots and summary statistics for Figure 1D.
############################################################

suppressWarnings(suppressMessages({
  library(tidyverse)
  library(readxl)
  library(cowplot)
}))

set.seed(2025)

# ==========================================================
# 1. File paths and parameters
# ==========================================================

feature_file  <- "C:/Users/chifang/OneDrive - University of Helsinki/Manuscripts/Tommi_Vatanen/3Figures_share/1C_3Features_long_correct.xlsx"
metadata_file <- "C:/Users/chifang/OneDrive - University of Helsinki/Manuscripts/Tommi_Vatanen/3Figures_share/metadata_924_685.xlsx"

output_dir <- "C:/Users/chifang/OneDrive - University of Helsinki/Manuscripts/Tommi_Vatanen/3Figures_share/1D_abundance_validation_results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pseudocount <- 1e-6

colors <- c(
  "Control" = "#009E73",
  "PD"      = "#E64B35"
)

# ==========================================================
# 2. Load and prepare feature data
# ==========================================================

cat("Loading feature table...\n")

feat_raw <- read_excel(feature_file)

names(feat_raw)[1] <- "feature"

feat_long <- feat_raw %>%
  pivot_longer(
    cols = -feature,
    names_to = "SampleID",
    values_to = "Abundance"
  ) %>%
  mutate(
    SampleID  = as.character(SampleID),
    Abundance = as.numeric(Abundance),
    Feature   = factor(feature)
  )

if (anyNA(feat_long$Abundance)) {
  stop("Feature table contains NA values after numeric conversion. Please check the input file.")
}

# ==========================================================
# 3. Load metadata
# ==========================================================

cat("Loading metadata...\n")

meta <- read_excel(metadata_file) %>%
  mutate(
    SampleID = as.character(SampleID),
    Group    = factor(Group, levels = c("Control", "PD")),
    Cohort   = factor(Cohort)
  ) %>%
  distinct(SampleID, .keep_all = TRUE)

# ==========================================================
# 4. Merge feature table and metadata
# ==========================================================

plot_df <- feat_long %>%
  inner_join(meta, by = "SampleID") %>%
  mutate(
    Abundance_log10 = log10(Abundance + pseudocount),
    Detected        = Abundance > 0
  )

cat("Matched samples:", length(unique(plot_df$SampleID)), "\n")
cat("Features:\n")
print(unique(plot_df$Feature))

write_csv(
  plot_df,
  file.path(output_dir, "figure1d_merged_plotting_data.csv")
)

# ==========================================================
# 5. Statistical analysis
# ==========================================================

cat("Performing Wilcoxon rank-sum tests...\n")

stat_df <- plot_df %>%
  group_by(Feature) %>%
  summarise(
    Wilcoxon_p = wilcox.test(Abundance ~ Group)$p.value,
    Wilcoxon_statistic = as.numeric(wilcox.test(Abundance ~ Group)$statistic),

    Mean_Control = mean(Abundance[Group == "Control"], na.rm = TRUE),
    Mean_PD      = mean(Abundance[Group == "PD"], na.rm = TRUE),

    Median_Control = median(Abundance[Group == "Control"], na.rm = TRUE),
    Median_PD      = median(Abundance[Group == "PD"], na.rm = TRUE),

    Detection_Control = mean(Abundance[Group == "Control"] > 0, na.rm = TRUE) * 100,
    Detection_PD      = mean(Abundance[Group == "PD"] > 0, na.rm = TRUE) * 100,

    n_Control = sum(Group == "Control"),
    n_PD      = sum(Group == "PD"),

    .groups = "drop"
  ) %>%
  mutate(
    Wilcoxon_adj = p.adjust(Wilcoxon_p, method = "BH"),

    Fold_Change = (Mean_PD + pseudocount) / (Mean_Control + pseudocount),
    Log2_Fold_Change = log2(Fold_Change),

    Detection_Diff = Detection_PD - Detection_Control,

    Significance = case_when(
      Wilcoxon_adj < 0.001 ~ "***",
      Wilcoxon_adj < 0.01  ~ "**",
      Wilcoxon_adj < 0.05  ~ "*",
      TRUE ~ "ns"
    ),

    Direction = ifelse(Mean_PD > Mean_Control, "Increased in PD", "Decreased in PD")
  )

write_csv(
  stat_df,
  file.path(output_dir, "figure1d_wilcoxon_statistics.csv")
)

# ==========================================================
# 6. Create violin plots
# ==========================================================

cat("Creating violin plots...\n")

feature_labels <- c(
  "Desulfovibrio_relab"    = expression(italic("Desulfovibrio")),
  "Ecoli_relab"            = expression(italic("E. coli")),
  "Ecoli_curligene_relab"  = expression(italic("E. coli")~"curli genes")
)

create_violin_plot <- function(data, feature_name, stat_info) {

  feature_data  <- data %>% filter(Feature == feature_name)
  feature_stats <- stat_info %>% filter(Feature == feature_name)

  ymax <- max(feature_data$Abundance_log10, na.rm = TRUE)

  ggplot(feature_data, aes(x = Group, y = Abundance_log10, fill = Group)) +

    geom_violin(
      trim = FALSE,
      scale = "width",
      color = "grey30",
      alpha = 0.85,
      adjust = 1.5
    ) +

    geom_boxplot(
      width = 0.14,
      outlier.shape = NA,
      color = "black",
      fill = "white",
      alpha = 0.45,
      linewidth = 0.45
    ) +

    geom_jitter(
      width = 0.10,
      height = 0,
      alpha = 0.18,
      size = 0.45,
      color = "grey20"
    ) +

    annotate(
      "text",
      x = 1.5,
      y = ymax + 0.25,
      label = sprintf("Wilcoxon p = %.2e", feature_stats$Wilcoxon_p),
      size = 3.4,
      hjust = 0.5
    ) +

    scale_fill_manual(values = colors) +

    scale_y_continuous(
      expand = expansion(mult = c(0.04, 0.12))
    ) +

    labs(
      x = NULL,
      y = expression(log[10]~"(relative abundance + 10"^{-6}*")"),
      title = feature_labels[[as.character(feature_name)]]
    ) +

    theme_classic(base_size = 12) +

    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      axis.line = element_line(linewidth = 0.45),
      axis.ticks = element_line(linewidth = 0.45)
    )
}

plots <- list()

for (feature in unique(plot_df$Feature)) {
  plots[[as.character(feature)]] <- create_violin_plot(plot_df, feature, stat_df)
}

combined_plot <- plot_grid(
  plotlist = plots,
  ncol = 3,
  align = "hv"
)

ggsave(
  file.path(output_dir, "figure1d_abundance_violin_plots.pdf"),
  combined_plot,
  width = 15,
  height = 4.8,
  device = "pdf"
)

ggsave(
  file.path(output_dir, "figure1d_abundance_violin_plots.tiff"),
  combined_plot,
  width = 15,
  height = 4.8,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

cat("Figure 1D analysis finished successfully.\n")
cat("Results saved to:", output_dir, "\n")
