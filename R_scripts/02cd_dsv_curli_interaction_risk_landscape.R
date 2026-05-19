############################################################
# Figure 2C-D: Desulfovibrio-curli interaction and PD risk
#
# Aim:
# Evaluate whether the interaction between Desulfovibrio and
# curli-gene abundance is associated with Parkinson's disease
# status after accounting for cohort-level heterogeneity.
#
# Main analyses:
# 1. Continuous interaction model:
#    PD status ~ Desulfovibrio * curli genes + (1 | Cohort)
#
# 2. Comparator interaction model:
#    PD status ~ Desulfovibrio * Escherichia coli + (1 | Cohort)
#
# 3. Predicted PD-risk landscape based on the Desulfovibrio-curli
#    interaction model.
#
# Additional output:
# A quadrant-based combined exposure analysis is generated as a
# sensitivity analysis.
############################################################

suppressWarnings(suppressMessages({
  library(tidyverse)
  library(lme4)
  library(broom.mixed)
  library(writexl)
  library(patchwork)
  library(viridis)
}))

set.seed(2025)

# ----------------------------------------------------------
# 1. Path configuration
# ----------------------------------------------------------

dir_base <- "C:/Users/chifang/OneDrive - University of Helsinki/Manuscripts/Tommi_Vatanen/3Figures_share"

dir_in <- file.path(dir_base, "Figure2_Outputs", "Step0_Preprocessed")
dir_out <- file.path(dir_base, "Figure2_Outputs", "Fig2CD_interaction_risk_landscape")

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

file_input <- file.path(dir_in, "df0_preprocessed.rds")

stopifnot(file.exists(file_input))

cat("Input file:", file_input, "\n")
cat("Output directory:", dir_out, "\n")

# ----------------------------------------------------------
# 2. Load preprocessed data
# ----------------------------------------------------------

df <- readRDS(file_input)

required_cols <- c(
  "SampleID",
  "Group",
  "Cohort",
  "DSV_raw",
  "Curli_raw",
  "Ecoli_raw"
)

missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    "Input data missing required columns: ",
    paste(missing_cols, collapse = ", "),
    "\nPlease check the Figure 2 preprocessing script."
  )
}

df <- df %>%
  mutate(
    SampleID = as.character(SampleID),
    Cohort = factor(Cohort),
    Group = factor(Group, levels = c("Control", "PD")),
    PD_bin = if_else(Group == "PD", 1L, 0L),
    DSV_raw = as.numeric(DSV_raw),
    Curli_raw = as.numeric(Curli_raw),
    Ecoli_raw = as.numeric(Ecoli_raw)
  )

if (anyNA(df$Group)) {
  stop("Group contains NA values. Please check group labels.")
}

# ----------------------------------------------------------
# 3. Log transformation and standardization
# ----------------------------------------------------------

get_half_min_pseudocount <- function(x) {
  nz <- x[x > 0 & is.finite(x)]
  if (length(nz) == 0) return(1e-12)
  min(nz) / 2
}

pc_dsv <- get_half_min_pseudocount(df$DSV_raw)
pc_curli <- get_half_min_pseudocount(df$Curli_raw)
pc_ecoli <- get_half_min_pseudocount(df$Ecoli_raw)

cat("Pseudocount DSV  :", pc_dsv, "\n")
cat("Pseudocount Curli:", pc_curli, "\n")
cat("Pseudocount Ecoli:", pc_ecoli, "\n")

df <- df %>%
  mutate(
    DSV_log10 = log10(DSV_raw + pc_dsv),
    Curli_log10 = log10(Curli_raw + pc_curli),
    Ecoli_log10 = log10(Ecoli_raw + pc_ecoli),

    DSV_z = as.numeric(scale(DSV_log10)),
    Curli_z = as.numeric(scale(Curli_log10)),
    Ecoli_z = as.numeric(scale(Ecoli_log10))
  )

cat("Group counts:\n")
print(table(df$Group))

cat("Cohort counts:\n")
print(table(df$Cohort))

# ----------------------------------------------------------
# 4. Helper functions
# ----------------------------------------------------------

format_p_value <- function(p) {
  case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "< 0.001",
    p < 0.01 ~ sprintf("%.3f", p),
    TRUE ~ sprintf("%.2f", p)
  )
}

extract_fixed_effects <- function(model, model_label) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    mutate(
      Model = model_label,
      OR = exp(estimate),
      CI_low = exp(estimate - 1.96 * std.error),
      CI_high = exp(estimate + 1.96 * std.error),
      P_formatted = format_p_value(p.value)
    ) %>%
    select(
      Model,
      term,
      estimate,
      std.error,
      statistic,
      p.value,
      P_formatted,
      OR,
      CI_low,
      CI_high
    )
}

extract_interaction <- function(model, model_label) {
  extract_fixed_effects(model, model_label) %>%
    filter(str_detect(term, ":"))
}

extract_model_fit <- function(model, model_label) {
  tibble(
    Model = model_label,
    AIC = AIC(model),
    BIC = BIC(model),
    logLik = as.numeric(logLik(model)),
    nobs = nobs(model),
    is_singular = lme4::isSingular(model)
  )
}

# ----------------------------------------------------------
# 5. Continuous interaction models
# ----------------------------------------------------------

m_target <- glmer(
  PD_bin ~ DSV_z * Curli_z + (1 | Cohort),
  data = df,
  family = binomial,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

m_comparator <- glmer(
  PD_bin ~ DSV_z * Ecoli_z + (1 | Cohort),
  data = df,
  family = binomial,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

tbl_interaction <- bind_rows(
  extract_interaction(m_target, "DSV x Curli"),
  extract_interaction(m_comparator, "DSV x E. coli")
)

tbl_fixed_effects <- bind_rows(
  extract_fixed_effects(m_target, "DSV x Curli"),
  extract_fixed_effects(m_comparator, "DSV x E. coli")
)

tbl_model_fit <- bind_rows(
  extract_model_fit(m_target, "DSV x Curli"),
  extract_model_fit(m_comparator, "DSV x E. coli")
)

# ----------------------------------------------------------
# 6. Figure 2C: interaction odds-ratio plot
# ----------------------------------------------------------

p_interaction <- ggplot(
  tbl_interaction,
  aes(
    x = Model,
    y = OR,
    ymin = CI_low,
    ymax = CI_high,
    color = Model
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_pointrange(
    size = 1.1,
    linewidth = 1.1
  ) +
  geom_text(
    aes(
      label = paste0("P = ", P_formatted),
      y = CI_high * 1.15
    ),
    size = 4,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_y_log10(
    limits = c(0.4, max(tbl_interaction$CI_high, na.rm = TRUE) * 1.35)
  ) +
  coord_flip() +
  scale_color_manual(
    values = c(
      "DSV x Curli" = "#D73027",
      "DSV x E. coli" = "#4575B4"
    )
  ) +
  labs(
    title = "Continuous interaction effects on PD risk",
    subtitle = "Mixed-effects logistic regression",
    x = NULL,
    y = "Odds ratio (95% CI, log scale)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

# ----------------------------------------------------------
# 7. Figure 2D: predicted PD-risk landscape
# ----------------------------------------------------------

grid_df <- expand.grid(
  DSV_z = seq(-2.5, 2.5, length.out = 100),
  Curli_z = seq(-2.5, 2.5, length.out = 100)
)

grid_df$Cohort <- levels(df$Cohort)[1]

grid_df$PD_prob <- predict(
  m_target,
  newdata = grid_df,
  type = "response",
  re.form = NA,
  allow.new.levels = TRUE
)

set.seed(20251224)

df_points <- df %>%
  slice_sample(n = min(300, nrow(df)))

p_landscape <- ggplot(grid_df, aes(x = DSV_z, y = Curli_z)) +
  geom_tile(aes(fill = PD_prob)) +
  geom_contour(
    aes(z = PD_prob),
    breaks = c(0.1, 0.2, 0.3, 0.4, 0.5),
    color = "white",
    alpha = 0.5
  ) +
  geom_point(
    data = df_points,
    aes(x = DSV_z, y = Curli_z, color = Group),
    alpha = 0.25,
    size = 1,
    inherit.aes = FALSE
  ) +
  scale_fill_viridis(
    option = "plasma",
    limits = c(0, 1),
    name = "Pr(PD)"
  ) +
  scale_color_manual(
    values = c(
      "Control" = "#74ADD1",
      "PD" = "#D73027"
    ),
    guide = "none"
  ) +
  labs(
    title = "DSV x Curli",
    subtitle = "Predicted PD probability",
    x = "DSV abundance (z-score)",
    y = "Curli abundance (z-score)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# ----------------------------------------------------------
# 8. Sensitivity analysis: quadrant-based exposure groups
# ----------------------------------------------------------

df_quadrant <- df %>%
  group_by(Cohort) %>%
  mutate(
    DSV_high = DSV_raw > median(DSV_raw, na.rm = TRUE),
    Curli_high = Curli_raw > median(Curli_raw, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    Exposure_code = case_when(
      !DSV_high & !Curli_high ~ "Both_Low",
      DSV_high & !Curli_high ~ "DSV_Only",
      !DSV_high & Curli_high ~ "Curli_Only",
      DSV_high & Curli_high ~ "Both_High"
    ),
    Exposure_code = factor(
      Exposure_code,
      levels = c(
        "Both_Low",
        "DSV_Only",
        "Curli_Only",
        "Both_High"
      )
    ),
    Exposure_label = recode(
      Exposure_code,
      Both_Low = "Both Low",
      DSV_Only = "DSV Only",
      Curli_Only = "Curli Only",
      Both_High = "Both High"
    )
  )

exposure_counts <- df_quadrant %>%
  count(Exposure_code, Exposure_label) %>%
  mutate(
    Label = paste0(Exposure_label, "\n(n = ", n, ")")
  )

m_quadrant <- glmer(
  PD_bin ~ Exposure_code + (1 | Cohort),
  data = df_quadrant,
  family = binomial,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

tbl_quadrant <- broom.mixed::tidy(m_quadrant, effects = "fixed") %>%
  filter(str_detect(term, "Exposure_code")) %>%
  mutate(
    Exposure_code = str_remove(term, "Exposure_code"),
    OR = exp(estimate),
    CI_low = exp(estimate - 1.96 * std.error),
    CI_high = exp(estimate + 1.96 * std.error),
    P_formatted = format_p_value(p.value)
  ) %>%
  left_join(exposure_counts, by = "Exposure_code")

p_quadrant <- ggplot(
  tbl_quadrant,
  aes(
    y = fct_reorder(Label, OR),
    x = OR,
    xmin = CI_low,
    xmax = CI_high,
    color = Exposure_label
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_errorbar(
    orientation = "y",
    height = 0.2,
    linewidth = 1.1
  ) +
  geom_point(size = 4) +
  geom_text(
    aes(label = sprintf("OR = %.2f", OR)),
    vjust = -1.2,
    size = 4,
    fontface = "bold"
  ) +
  geom_text(
    aes(
      label = paste0("P = ", P_formatted),
      x = CI_high * 1.12
    ),
    size = 4,
    fontface = "bold"
  ) +
  scale_x_log10(
    limits = c(0.4, max(tbl_quadrant$CI_high, na.rm = TRUE) * 1.35)
  ) +
  scale_color_manual(
    values = c(
      "Both Low" = "#4575B4",
      "DSV Only" = "#74ADD1",
      "Curli Only" = "#FDAE61",
      "Both High" = "#D73027"
    )
  ) +
  labs(
    title = "Categorical risk stratification",
    subtitle = "Within-cohort median splits; reference: Both Low",
    x = "Odds ratio (95% CI, log scale)",
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

# ----------------------------------------------------------
# 9. Save statistics
# ----------------------------------------------------------

write_xlsx(
  list(
    interaction_terms = tbl_interaction,
    fixed_effects = tbl_fixed_effects,
    model_fit = tbl_model_fit,
    quadrant_terms = tbl_quadrant,
    quadrant_counts = exposure_counts
  ),
  file.path(dir_out, "figure2cd_model_statistics.xlsx")
)

# ----------------------------------------------------------
# 10. Save figures
# ----------------------------------------------------------

ggsave(
  filename = file.path(dir_out, "figure2c_interaction_odds_ratio.pdf"),
  plot = p_interaction,
  width = 6.8,
  height = 4.8,
  device = "pdf"
)

ggsave(
  filename = file.path(dir_out, "figure2c_interaction_odds_ratio.tiff"),
  plot = p_interaction,
  width = 6.8,
  height = 4.8,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(dir_out, "figure2d_dsv_curli_pd_risk_landscape.pdf"),
  plot = p_landscape,
  width = 6.8,
  height = 5.5,
  device = "pdf"
)

ggsave(
  filename = file.path(dir_out, "figure2d_dsv_curli_pd_risk_landscape.tiff"),
  plot = p_landscape,
  width = 6.8,
  height = 5.5,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

combined_main <- p_interaction + p_landscape +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    tag_levels = "c",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 18)
    )
  ) &
  theme(legend.position = "right")

ggsave(
  filename = file.path(dir_out, "figure2cd_combined_interaction_landscape.pdf"),
  plot = combined_main,
  width = 12,
  height = 6,
  device = "pdf"
)

ggsave(
  filename = file.path(dir_out, "figure2cd_combined_interaction_landscape.tiff"),
  plot = combined_main,
  width = 12,
  height = 6,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(dir_out, "supplementary_quadrant_combined_exposure_or.pdf"),
  plot = p_quadrant,
  width = 6.8,
  height = 4.8,
  device = "pdf"
)

ggsave(
  filename = file.path(dir_out, "supplementary_quadrant_combined_exposure_or.tiff"),
  plot = p_quadrant,
  width = 6.8,
  height = 4.8,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

saveRDS(
  list(
    data = df,
    m_target = m_target,
    m_comparator = m_comparator,
    m_quadrant = m_quadrant,
    interaction_terms = tbl_interaction,
    fixed_effects = tbl_fixed_effects,
    model_fit = tbl_model_fit,
    quadrant_terms = tbl_quadrant,
    p_interaction = p_interaction,
    p_landscape = p_landscape,
    p_quadrant = p_quadrant
  ),
  file.path(dir_out, "figure2cd_analysis_objects.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(dir_out, "sessionInfo.txt")
)

cat("Figure 2C-D interaction risk analysis completed successfully.\n")
cat("Outputs written to:", dir_out, "\n")
