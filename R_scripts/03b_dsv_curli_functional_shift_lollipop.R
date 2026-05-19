############################################################
# Figure 3B: DSV-curli-associated functional pathway shifts
#
# Aim:
# Visualize the top functional pathways associated with the
# Desulfovibrio-curli interaction effect using a diverging
# lollipop plot.
#
# Input:
# Fig3B_Final_DSVxCurli_abundance_all.csv
#
# Main steps:
# 1. Load pathway-level differential association results.
# 2. Identify the q-value column used for multiple-testing correction.
# 3. Retain pathways with q < 0.1.
# 4. Select the top 15 positively and top 15 negatively associated
#    pathways based on absolute coefficient size.
# 5. Generate a diverging lollipop plot for Figure 3B.
#
# Output:
# Figure 3B lollipop plot and the selected pathway table.
############################################################

suppressWarnings(suppressMessages({
  library(ggplot2)
  library(stringr)
}))

set.seed(2025)

# ==========================================================
# 1. File paths
# ==========================================================

file_path <- "C:/Users/chifang/OneDrive - University of Helsinki/Manuscripts/Tommi_Vatanen/3Figures_share/Figure3_Outputs/Fig3B/Fig3B_Final_DSVxCurli_abundance_all.csv"

dir_out <- dirname(file_path)

if (!file.exists(file_path)) {
  stop("Input file not found: ", file_path)
}

cat("Input file:", file_path, "\n")
cat("Output directory:", dir_out, "\n")

# ==========================================================
# 2. Load differential pathway results
# ==========================================================

df <- read.csv(file_path, stringsAsFactors = FALSE)

cat("Loaded rows:", nrow(df), "\n")

required_cols <- c("feature", "coef")
missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    "Input table missing required columns: ",
    paste(missing_cols, collapse = ", "),
    "\nDetected columns: ",
    paste(colnames(df), collapse = ", ")
  )
}

# ==========================================================
# 3. Identify q-value column
# ==========================================================

if ("qval_joint" %in% colnames(df)) {
  df$target_q <- df$qval_joint
} else if ("qval" %in% colnames(df)) {
  df$target_q <- df$qval
} else if ("qval_individual" %in% colnames(df)) {
  df$target_q <- df$qval_individual
} else {
  stop("No q-value column found. Expected one of: qval_joint, qval, qval_individual.")
}

df$coef <- as.numeric(df$coef)
df$target_q <- as.numeric(df$target_q)

df <- df[!is.na(df$coef) & !is.na(df$target_q), ]

cat("Rows after removing missing coef/q values:", nrow(df), "\n")

# ==========================================================
# 4. Filter significant pathways
# ==========================================================

q_threshold <- 0.1

df_sig <- df[df$target_q < q_threshold, ]

cat("Significant pathways with q <", q_threshold, ":", nrow(df_sig), "\n")

if (nrow(df_sig) == 0) {
  stop("No pathways passed the q-value threshold.")
}

# ==========================================================
# 5. Select top positive and negative pathways
# ==========================================================

n_top <- 15

up_rows <- df_sig[df_sig$coef > 0, ]
down_rows <- df_sig[df_sig$coef < 0, ]

up_rows <- up_rows[order(abs(up_rows$coef), decreasing = TRUE), ]
down_rows <- down_rows[order(abs(down_rows$coef), decreasing = TRUE), ]

top_up <- head(up_rows, n_top)
top_down <- head(down_rows, n_top)

plot_data <- rbind(top_up, top_down)

cat("Selected pathways for plotting:", nrow(plot_data), "\n")
cat("Positive pathways:", nrow(top_up), "\n")
cat("Negative pathways:", nrow(top_down), "\n")

if (nrow(plot_data) == 0) {
  stop("No positive or negative pathways available for plotting.")
}

# ==========================================================
# 6. Prepare pathway labels
# ==========================================================

plot_data$Direction <- ifelse(
  plot_data$coef > 0,
  "Positive association",
  "Negative association"
)

plot_data$CleanName <- gsub("^[^:]+: ", "", plot_data$feature)
plot_data$CleanName <- stringr::str_wrap(plot_data$CleanName, width = 60)

plot_data$CleanName <- factor(
  plot_data$CleanName,
  levels = plot_data$CleanName[order(plot_data$coef)]
)

plot_data$minus_log10_q <- -log10(plot_data$target_q)

# ==========================================================
# 7. Save selected pathway table
# ==========================================================

write.csv(
  plot_data,
  file.path(dir_out, "figure3b_lollipop_selected_pathways.csv"),
  row.names = FALSE
)

# ==========================================================
# 8. Generate Figure 3B lollipop plot
# ==========================================================

p_lollipop <- ggplot(plot_data, aes(x = CleanName, y = coef)) +
  geom_segment(
    aes(
      x = CleanName,
      xend = CleanName,
      y = 0,
      yend = coef,
      color = Direction
    ),
    linewidth = 1.2,
    alpha = 0.85
  ) +
  geom_point(
    aes(
      color = Direction,
      size = minus_log10_q
    ),
    alpha = 1
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.4
  ) +
  scale_color_manual(
    values = c(
      "Positive association" = "#D73027",
      "Negative association" = "#4575B4"
    )
  ) +
  coord_flip() +
  labs(
    title = "DSV-curli-associated functional pathway shifts",
    subtitle = "Top positively and negatively associated pathways",
    x = NULL,
    y = "Effect size coefficient",
    color = NULL,
    size = expression(-log[10]~"(q-value)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(
      size = 9,
      color = "black",
      lineheight = 0.85
    ),
    axis.text.x = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10,
      unit = "pt"
    )
  )

# ==========================================================
# 9. Save outputs
# ==========================================================

ggsave(
  filename = file.path(dir_out, "figure3b_dsv_curli_functional_shift_lollipop.pdf"),
  plot = p_lollipop,
  width = 12,
  height = 10,
  device = "pdf"
)

ggsave(
  filename = file.path(dir_out, "figure3b_dsv_curli_functional_shift_lollipop.tiff"),
  plot = p_lollipop,
  width = 12,
  height = 10,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(dir_out, "sessionInfo_figure3b_lollipop.txt")
)

cat("Figure 3B lollipop plot completed successfully.\n")
cat("Outputs written to:", dir_out, "\n")

# ==========================================================
# 10. Export Figure 3B source data and summary statistics
# ==========================================================

cat("Exporting Figure 3B source data and summary statistics...\n")

fig3b_plot_data <- plot_data %>%
  transmute(
    Feature_ID = feature,
    Pathway_Name = as.character(CleanName),
    Direction = Direction,
    Beta_Coefficient = coef,
    Abs_Effect_Size = abs(coef),
    Q_value = target_q,
    NegLog10_Q = -log10(target_q)
  ) %>%
  arrange(desc(Abs_Effect_Size))

fig3b_all_significant <- df_sig %>%
  mutate(
    Direction = ifelse(
      coef > 0,
      "Positive association",
      "Negative association"
    )
  ) %>%
  transmute(
    Feature_ID = feature,
    Beta_Coefficient = coef,
    Abs_Effect_Size = abs(coef),
    Q_value = target_q,
    Direction = Direction
  ) %>%
  arrange(desc(Abs_Effect_Size))

fig3b_summary <- tibble(
  Metric = c(
    "Total pathways tested",
    "Significant pathways (q < 0.1)",
    "Positive pathways",
    "Negative pathways",
    "Pathways shown in Figure 3B"
  ),
  Value = c(
    nrow(df),
    nrow(df_sig),
    sum(df_sig$coef > 0, na.rm = TRUE),
    sum(df_sig$coef < 0, na.rm = TRUE),
    nrow(plot_data)
  )
)

write_xlsx(
  list(
    Figure3B_plot_data = fig3b_plot_data,
    Significant_pathways_q0.1 = fig3b_all_significant,
    Summary_statistics = fig3b_summary
  ),
  file.path(dir_out, "figure3b_source_data_and_statistics.xlsx")
)

write.csv(
  fig3b_plot_data,
  file.path(dir_out, "figure3b_lollipop_plot_data.csv"),
  row.names = FALSE
)

write.csv(
  fig3b_all_significant,
  file.path(dir_out, "figure3b_all_significant_pathways_q0.1.csv"),
  row.names = FALSE
)

write.csv(
  fig3b_summary,
  file.path(dir_out, "figure3b_summary_statistics.csv"),
  row.names = FALSE
)

print(fig3b_summary)

cat("Figure 3B source data exported successfully.\n")
