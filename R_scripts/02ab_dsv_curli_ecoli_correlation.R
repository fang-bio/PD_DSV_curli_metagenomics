############################################################
# Figure 2A-B: Specific ecological coupling analysis
#
# Aim:
# Evaluate whether Desulfovibrio abundance shows stronger
# disease-specific coupling with bacterial curli genes than with
# total Escherichia coli abundance.
#
# Analysis:
# 1. Desulfovibrio vs. curli genes among samples positive for both
# 2. Desulfovibrio vs. Escherichia coli among samples positive for both
#
# Method:
# Spearman correlation, stratified by disease group.
#
# Input:
# df0_preprocessed.rds from Figure 2 preprocessing.
#
# Output:
# Correlation statistics and scatter plots for Figure 2A-B.
############################################################

suppressWarnings(suppressMessages({
  library(tidyverse)
  library(patchwork)
  library(writexl)
}))

set.seed(2025)

# ----------------------------------------------------------
# 1. Path configuration
# ----------------------------------------------------------

dir_base <- "C:/Users/chifang/OneDrive - University of Helsinki/Manuscripts/Tommi_Vatanen/3Figures_share"

dir_in <- file.path(dir_base, "Figure2_Outputs", "Step0_Preprocessed")
dir_out <- file.path(dir_base, "Figure2_Outputs", "Fig2AB_specific_coupling")

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

file_input <- file.path(dir_in, "df0_preprocessed.rds")

stopifnot(file.exists(file_input))

cat("Input file:", file_input, "\n")
cat("Output directory:", dir_out, "\n")

# ----------------------------------------------------------
# 2. Load preprocessed data
# ----------------------------------------------------------

df0 <- readRDS(file_input)

required_cols <- c(
  "SampleID",
  "Group",
  "Cohort",
  "DSV_raw",
  "Ecoli_raw",
  "Curli_raw"
)

missing_cols <- setdiff(required_cols, colnames(df0))

if (length(missing_cols) > 0) {
  stop(
    "Input data missing required columns: ",
    paste(missing_cols, collapse = ", "),
    "\nPlease check the Figure 2 preprocessing script."
  )
}

df0 <- df0 %>%
  mutate(
    SampleID = as.character(SampleID),
    Group = factor(Group, levels = c("Control", "PD")),
    Cohort = factor(Cohort),
    DSV_raw = as.numeric(DSV_raw),
    Ecoli_raw = as.numeric(Ecoli_raw),
    Curli_raw = as.numeric(Curli_raw)
  )

if (anyNA(df0$Group)) {
  stop("Group contains NA values. Please check group labels.")
}

# ----------------------------------------------------------
# 3. Prepare carrier-only subsets
# ----------------------------------------------------------

df_curli <- df0 %>%
  filter(DSV_raw > 0, Curli_raw > 0) %>%
  mutate(
    log_DSV = log10(DSV_raw),
    log_Curli = log10(Curli_raw)
  )

df_ecoli <- df0 %>%
  filter(DSV_raw > 0, Ecoli_raw > 0) %>%
  mutate(
    log_DSV = log10(DSV_raw),
    log_Ecoli = log10(Ecoli_raw)
  )

cat("Samples for DSV vs curli:", nrow(df_curli), "\n")
cat("Samples for DSV vs E. coli:", nrow(df_ecoli), "\n")

cat("Group counts for DSV vs curli:\n")
print(table(df_curli$Group))

cat("Group counts for DSV vs E. coli:\n")
print(table(df_ecoli$Group))

# ----------------------------------------------------------
# 4. Spearman correlation analysis
# ----------------------------------------------------------

calculate_spearman <- function(data, x_col, y_col, comparison_name) {

  data %>%
    group_by(Group) %>%
    summarise(
      Comparison = comparison_name,
      N = n(),
      Rho = ifelse(
        N >= 3 &&
          sd(.data[[x_col]], na.rm = TRUE) > 0 &&
          sd(.data[[y_col]], na.rm = TRUE) > 0,
        cor(.data[[x_col]], .data[[y_col]], method = "spearman"),
        NA_real_
      ),
      P_value = ifelse(
        N >= 3 &&
          sd(.data[[x_col]], na.rm = TRUE) > 0 &&
          sd(.data[[y_col]], na.rm = TRUE) > 0,
        cor.test(
          .data[[x_col]],
          .data[[y_col]],
          method = "spearman",
          exact = FALSE
        )$p.value,
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      P_fmt = case_when(
        is.na(P_value) ~ "NA",
        P_value < 0.001 ~ "< 0.001",
        TRUE ~ sprintf("%.3f", P_value)
      ),
      Rho_fmt = ifelse(is.na(Rho), "NA", sprintf("%.2f", Rho)),
      Label = paste0(
        Group,
        ": rho = ",
        Rho_fmt,
        ", P = ",
        P_fmt,
        ", N = ",
        N
      ),
      Significance = case_when(
        is.na(P_value) ~ "NA",
        P_value < 0.001 ~ "***",
        P_value < 0.01 ~ "**",
        P_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
}

stats_curli <- calculate_spearman(
  data = df_curli,
  x_col = "log_DSV",
  y_col = "log_Curli",
  comparison_name = "Desulfovibrio vs curli genes"
)

stats_ecoli <- calculate_spearman(
  data = df_ecoli,
  x_col = "log_DSV",
  y_col = "log_Ecoli",
  comparison_name = "Desulfovibrio vs Escherichia coli"
)

stats_all <- bind_rows(stats_curli, stats_ecoli) %>%
  select(Comparison, Group, N, Rho, P_value, P_fmt, Significance, Label)

write_xlsx(
  stats_all,
  file.path(dir_out, "figure2ab_spearman_statistics.xlsx")
)

print(stats_all)

# ----------------------------------------------------------
# 5. Plot settings
# ----------------------------------------------------------

pal <- c(
  "Control" = "#4575B4",
  "PD" = "#D73027"
)

base_theme <- theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.4),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

make_annotation_df <- function(data, stats_tbl, x_col, y_col) {

  range_x <- range(data[[x_col]], na.rm = TRUE)
  range_y <- range(data[[y_col]], na.rm = TRUE)

  stats_tbl %>%
    mutate(
      x = range_x[1] + 0.04 * diff(range_x),
      y = case_when(
        Group == "Control" ~ range_y[2] - 0.06 * diff(range_y),
        Group == "PD" ~ range_y[2] - 0.16 * diff(range_y),
        TRUE ~ range_y[2] - 0.26 * diff(range_y)
      )
    )
}

ann_curli <- make_annotation_df(df_curli, stats_curli, "log_DSV", "log_Curli")
ann_ecoli <- make_annotation_df(df_ecoli, stats_ecoli, "log_DSV", "log_Ecoli")

# ----------------------------------------------------------
# 6. Generate scatter plots
# ----------------------------------------------------------

p_curli <- ggplot(df_curli, aes(x = log_DSV, y = log_Curli, color = Group)) +
  geom_point(alpha = 0.60, size = 2) +
  geom_smooth(
    method = "lm",
    aes(fill = Group),
    alpha = 0.15,
    linewidth = 1,
    se = TRUE
  ) +
  geom_text(
    data = ann_curli,
    aes(x = x, y = y, label = Label, color = Group),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.7,
    show.legend = FALSE
  ) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(
    title = "Desulfovibrio vs. curli genes",
    x = expression(log[10]~"(Desulfovibrio abundance)"),
    y = expression(log[10]~"(curli-gene abundance)")
  ) +
  base_theme

p_ecoli <- ggplot(df_ecoli, aes(x = log_DSV, y = log_Ecoli, color = Group)) +
  geom_point(alpha = 0.60, size = 2) +
  geom_smooth(
    method = "lm",
    aes(fill = Group),
    alpha = 0.15,
    linewidth = 1,
    se = TRUE
  ) +
  geom_text(
    data = ann_ecoli,
    aes(x = x, y = y, label = Label, color = Group),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.7,
    show.legend = FALSE
  ) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(
    title = "Desulfovibrio vs. Escherichia coli",
    x = expression(log[10]~"(Desulfovibrio abundance)"),
    y = expression(log[10]~"(Escherichia coli abundance)")
  ) +
  base_theme

combined_plot <- p_curli + p_ecoli +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 18)
    )
  ) &
  theme(legend.position = "bottom")

# ----------------------------------------------------------
# 7. Save outputs
# ----------------------------------------------------------

ggsave(
  filename = file.path(dir_out, "figure2ab_specific_coupling.pdf"),
  plot = combined_plot,
  width = 12,
  height = 6,
  device = "pdf"
)

ggsave(
  filename = file.path(dir_out, "figure2ab_specific_coupling.tiff"),
  plot = combined_plot,
  width = 12,
  height = 6,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

saveRDS(
  list(
    df_curli = df_curli,
    df_ecoli = df_ecoli,
    stats_all = stats_all,
    plot = combined_plot
  ),
  file.path(dir_out, "figure2ab_specific_coupling_objects.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(dir_out, "sessionInfo.txt")
)

cat("Figure 2A-B specific coupling analysis completed successfully.\n")
cat("Outputs written to:", dir_out, "\n")
