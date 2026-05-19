############################################################
# Supplementary Figure 2A: DSV-E. coli comparator landscape
#
# Aim:
# Visualize the predicted Parkinson's disease probability surface
# for the comparator interaction model between Desulfovibrio and
# total Escherichia coli abundance.
#
# Models:
# 1. Target model:
#    PD status ~ Desulfovibrio * curli genes + (1 | Cohort)
#
# 2. Comparator model:
#    PD status ~ Desulfovibrio * Escherichia coli + (1 | Cohort)
#
# Output:
# Predicted-risk heatmaps for the target and comparator models.
# The Desulfovibrio x Escherichia coli comparator landscape is
# used as Supplementary Figure 2A.
############################################################

suppressWarnings(suppressMessages({
  library(tidyverse)
  library(lme4)
  library(broom.mixed)
  library(writexl)
}))

set.seed(2025)

# ----------------------------------------------------------
# 1. Path configuration
# ----------------------------------------------------------

dir_base <- "C:/Users/chifang/OneDrive - University of Helsinki/Manuscripts/Tommi_Vatanen/3Figures_share"

dir_in <- file.path(dir_base, "Figure2_Outputs", "Step0_Preprocessed")
dir_out <- file.path(dir_base, "Figure2_Outputs", "SuppFig2A_dsv_ecoli_comparator_landscape")

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
    Group = factor(Group, levels = c("Control", "PD")),
    Cohort = factor(Cohort),
    PD_bin = if_else(Group == "PD", 1L, 0L),
    DSV_raw = as.numeric(DSV_raw),
    Curli_raw = as.numeric(Curli_raw),
    Ecoli_raw = as.numeric(Ecoli_raw)
  )

if (anyNA(df$Group)) {
  stop("Group contains NA values. Please check group labels.")
}

# ----------------------------------------------------------
# 3. Log transformation and independent standardization
# ----------------------------------------------------------

get_half_min_pseudocount <- function(x) {
  nz <- x[x > 0 & is.finite(x)]

  if (length(nz) == 0) {
    return(1e-12)
  }

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
    log_DSV = log10(DSV_raw + pc_dsv),
    log_Curli = log10(Curli_raw + pc_curli),
    log_Ecoli = log10(Ecoli_raw + pc_ecoli),

    DSV_z = as.numeric(scale(log_DSV)),
    Curli_z = as.numeric(scale(log_Curli)),
    Ecoli_z = as.numeric(scale(log_Ecoli))
  )

cat("Group counts:\n")
print(table(df$Group))

cat("Cohort counts:\n")
print(table(df$Cohort))

# ----------------------------------------------------------
# 4. Fit mixed-effects logistic regression models
# ----------------------------------------------------------

m_target <- glmer(
  PD_bin ~ DSV_z * Curli_z + (1 | Cohort),
  data = df,
  family = binomial,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

m_comparator <- glmer(
  PD_bin ~ DSV_z * Ecoli_z + (1 | Cohort),
  data = df,
  family = binomial,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

# ----------------------------------------------------------
# 5. Export model statistics
# ----------------------------------------------------------

format_p_value <- function(p) {
  case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "< 0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

extract_model_statistics <- function(model, model_label) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    mutate(
      Model = model_label,
      OR = exp(estimate),
      CI_low = exp(estimate - 1.96 * std.error),
      CI_high = exp(estimate + 1.96 * std.error),
      P_label = format_p_value(p.value)
    ) %>%
    select(
      Model,
      term,
      estimate,
      std.error,
      statistic,
      p.value,
      P_label,
      OR,
      CI_low,
      CI_high
    )
}

stats_target <- extract_model_statistics(
  m_target,
  "DSV x Curli"
)

stats_comparator <- extract_model_statistics(
  m_comparator,
  "DSV x E. coli"
)

stats_all <- bind_rows(
  stats_target,
  stats_comparator
)

model_fit <- tibble(
  Model = c("DSV x Curli", "DSV x E. coli"),
  AIC = c(AIC(m_target), AIC(m_comparator)),
  BIC = c(BIC(m_target), BIC(m_comparator)),
  logLik = c(as.numeric(logLik(m_target)), as.numeric(logLik(m_comparator))),
  nobs = c(nobs(m_target), nobs(m_comparator)),
  is_singular = c(
    lme4::isSingular(m_target),
    lme4::isSingular(m_comparator)
  )
)

write_xlsx(
  list(
    fixed_effects = stats_all,
    model_fit = model_fit
  ),
  file.path(dir_out, "suppfig2a_model_statistics.xlsx")
)

get_interaction_p_label <- function(stats_table) {
  p_label <- stats_table %>%
    filter(str_detect(term, ":")) %>%
    pull(P_label)

  if (length(p_label) == 0) {
    return("NA")
  }

  p_label[1]
}

pval_target <- get_interaction_p_label(stats_target)
pval_comparator <- get_interaction_p_label(stats_comparator)

# ----------------------------------------------------------
# 6. Generate predicted-risk landscapes
# ----------------------------------------------------------

create_prediction_grid <- function(x_name, y_name) {
  grid <- expand.grid(
    x = seq(-2.5, 2.5, length.out = 100),
    y = seq(-2.5, 2.5, length.out = 100)
  )

  colnames(grid) <- c(x_name, y_name)
  grid$Cohort <- levels(df$Cohort)[1]

  grid
}

predict_surface <- function(model, grid) {
  grid$predicted_PD_probability <- predict(
    model,
    newdata = grid,
    type = "response",
    re.form = NA,
    allow.new.levels = TRUE
  )

  grid
}

grid_target <- create_prediction_grid("DSV_z", "Curli_z") %>%
  predict_surface(m_target)

grid_comparator <- create_prediction_grid("DSV_z", "Ecoli_z") %>%
  predict_surface(m_comparator)

plot_surface <- function(grid, x_name, y_name, title, subtitle, x_lab, y_lab) {
  ggplot(grid, aes(x = .data[[x_name]], y = .data[[y_name]])) +
    geom_raster(aes(fill = predicted_PD_probability), interpolate = TRUE) +
    geom_contour(
      aes(z = predicted_PD_probability),
      color = "black",
      alpha = 0.20,
      linewidth = 0.40
    ) +
    scale_fill_gradientn(
      colors = c("#3288BD", "white", "#D53E4F"),
      limits = c(0, 1),
      name = "Pr(PD)"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_lab,
      y = y_lab
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right",
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      axis.line = element_blank(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1
      )
    )
}

p_target <- plot_surface(
  grid = grid_target,
  x_name = "DSV_z",
  y_name = "Curli_z",
  title = "DSV x Curli",
  subtitle = paste0("Interaction P = ", pval_target),
  x_lab = "DSV abundance (z-score)",
  y_lab = "Curli abundance (z-score)"
)

p_comparator <- plot_surface(
  grid = grid_comparator,
  x_name = "DSV_z",
  y_name = "Ecoli_z",
  title = "DSV x E. coli",
  subtitle = paste0("Interaction P = ", pval_comparator),
  x_lab = "DSV abundance (z-score)",
  y_lab = "Total E. coli abundance (z-score)"
)

# ----------------------------------------------------------
# 7. Save outputs
# ----------------------------------------------------------

ggsave(
  filename = file.path(dir_out, "target_dsv_curli_reference_landscape.pdf"),
  plot = p_target,
  width = 6.5,
  height = 5.5,
  device = "pdf"
)

ggsave(
  filename = file.path(dir_out, "target_dsv_curli_reference_landscape.tiff"),
  plot = p_target,
  width = 6.5,
  height = 5.5,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(dir_out, "suppfig2a_dsv_ecoli_comparator_landscape.pdf"),
  plot = p_comparator,
  width = 6.5,
  height = 5.5,
  device = "pdf"
)

ggsave(
  filename = file.path(dir_out, "suppfig2a_dsv_ecoli_comparator_landscape.tiff"),
  plot = p_comparator,
  width = 6.5,
  height = 5.5,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

write_xlsx(
  list(
    target_prediction_grid = grid_target,
    comparator_prediction_grid = grid_comparator
  ),
  file.path(dir_out, "suppfig2a_prediction_grids.xlsx")
)

saveRDS(
  list(
    data = df,
    m_target = m_target,
    m_comparator = m_comparator,
    fixed_effects = stats_all,
    model_fit = model_fit,
    grid_target = grid_target,
    grid_comparator = grid_comparator,
    p_target = p_target,
    p_comparator = p_comparator
  ),
  file.path(dir_out, "suppfig2a_comparator_landscape_objects.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(dir_out, "sessionInfo.txt")
)

cat("Supplementary Figure 2A comparator landscape analysis completed successfully.\n")
cat("Outputs written to:", dir_out, "\n")
