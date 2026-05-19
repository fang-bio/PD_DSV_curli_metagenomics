############################################################
# Taxonomic feature preprocessing from MetaPhlAn profiles
#
# Aim:
# Prepare targeted taxonomic features for downstream analyses and
# visualization by extracting genus-level Desulfovibrio abundance
# and species-level Escherichia coli abundance from merged MetaPhlAn
# profiles.
#
# Main steps:
# 1. Load the merged MetaPhlAn abundance table.
# 2. Extract and collapse genus-level abundance profiles.
# 3. Extract Desulfovibrio abundance from the genus-level table.
# 4. Extract and collapse species-level abundance profiles.
# 5. Extract Escherichia coli abundance from the species-level table.
# 6. Convert the combined Desulfovibrio and Escherichia coli table
#    into MaAsLin3-compatible metadata and feature tables.
#
# Input:
# 1. 1715_merged_metaphlan_data.csv
# 2. Step3_DSV_Ecoli_Long_With_Group.csv
#
# Output:
# Genus-level abundance table, species-level abundance table,
# Desulfovibrio abundance table, Escherichia coli abundance table,
# and MaAsLin3-formatted input files.
############################################################

suppressWarnings(suppressMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
}))

set.seed(2025)

# ==========================================================
# 1. Path configuration
# ==========================================================

dir_metaphlan <- "C:/Users/chifang/OneDrive - University of Helsinki/Manuscripts/Tommi Vatanen/metaphlan/Desulfovibrio"

dir_target <- "E:/OneDrive - University of Helsinki/Manuscripts/Tommi Vatanen/metaphlan/Desulfovibrio"

dir_maaslin <- "E:/OneDrive - University of Helsinki/Manuscripts/Tommi Vatanen/metaphlan/Desulfovibrio and Escherichia_coli"

file_metaphlan <- file.path(
  dir_metaphlan,
  "1715_merged_metaphlan_data.csv"
)

file_dsv_ecoli_long <- file.path(
  dir_maaslin,
  "Step3_DSV_Ecoli_Long_With_Group.csv"
)

stopifnot(file.exists(file_metaphlan))

dir.create(dir_metaphlan, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_target, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_maaslin, recursive = TRUE, showWarnings = FALSE)

cat("Start taxonomic feature preprocessing\n")
cat("MetaPhlAn merged file:", file_metaphlan, "\n")

# ==========================================================
# 2. Load merged MetaPhlAn abundance table
# ==========================================================

abundance_data <- read.csv(
  file_metaphlan,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"Taxa" %in% colnames(abundance_data)) {
  colnames(abundance_data)[1] <- "Taxa"
}

cat("Loaded MetaPhlAn table:", nrow(abundance_data), "taxa/features\n")
cat("Sample columns:", ncol(abundance_data) - 1, "\n")

# Convert abundance columns to numeric
abundance_data <- abundance_data %>%
  mutate(
    across(
      .cols = -Taxa,
      .fns = ~ suppressWarnings(as.numeric(.x))
    )
  )

if (anyNA(abundance_data[, -1])) {
  warning("Some abundance values were converted to NA. Please check the input table.")
}

# ==========================================================
# 3. Extract genus-level abundance profiles
# ==========================================================
# MetaPhlAn taxonomy strings can include multiple taxonomic levels.
# Genus-level rows are identified as entries containing g__ but not s__.

genus_data <- abundance_data %>%
  filter(grepl("g__", Taxa), !grepl("s__", Taxa)) %>%
  mutate(
    Genus = sub(".*g__", "", Taxa),
    Genus = sub("\\|.*", "", Genus)
  ) %>%
  select(Genus, everything(), -Taxa) %>%
  group_by(Genus) %>%
  summarise(
    across(
      .cols = everything(),
      .fns = ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

cat("Genus-level features:", nrow(genus_data), "\n")

output_genus <- file.path(
  dir_metaphlan,
  "Step2_Genus_Level_Data_Deduplicated.csv"
)

write.csv(
  genus_data,
  output_genus,
  row.names = FALSE
)

cat("Genus-level abundance table saved to:", output_genus, "\n")

# ==========================================================
# 4. Extract Desulfovibrio abundance
# ==========================================================

desulfovibrio_data <- genus_data %>%
  filter(Genus == "Desulfovibrio")

if (nrow(desulfovibrio_data) == 0) {
  warning("Desulfovibrio was not found in the genus-level abundance table.")
}

output_dsv <- file.path(
  dir_target,
  "Step3_Only_Desulfovibrio_level.csv"
)

write.csv(
  desulfovibrio_data,
  output_dsv,
  row.names = FALSE
)

cat("Desulfovibrio abundance table saved to:", output_dsv, "\n")

# ==========================================================
# 5. Extract species-level abundance profiles
# ==========================================================
# Species-level rows are identified as entries containing s__.
# Strain-level suffixes are removed before collapsing duplicated species.

species_data <- abundance_data %>%
  filter(grepl("s__", Taxa)) %>%
  mutate(
    Species = gsub("\\|t__.*", "", Taxa),
    Species = sub(".*s__", "", Species)
  ) %>%
  select(Species, everything(), -Taxa) %>%
  group_by(Species) %>%
  summarise(
    across(
      .cols = everything(),
      .fns = ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

cat("Species-level features:", nrow(species_data), "\n")

output_species <- file.path(
  dir_target,
  "Step2_Species_Level_Data_Deduplicated.csv"
)

write.csv(
  species_data,
  output_species,
  row.names = FALSE
)

cat("Species-level abundance table saved to:", output_species, "\n")

# ==========================================================
# 6. Extract Escherichia coli abundance
# ==========================================================

ecoli_data <- species_data %>%
  filter(Species == "Escherichia_coli")

if (nrow(ecoli_data) == 0) {
  warning("Escherichia_coli was not found in the species-level abundance table.")
}

output_ecoli <- file.path(
  dir_target,
  "Step3_Only_Escherichia_coli.csv"
)

write.csv(
  ecoli_data,
  output_ecoli,
  row.names = FALSE
)

cat("Escherichia coli abundance table saved to:", output_ecoli, "\n")

# ==========================================================
# 7. Convert combined DSV and E. coli table to MaAsLin3 format
# ==========================================================

if (file.exists(file_dsv_ecoli_long)) {

  cat("Preparing MaAsLin3 input files from:", file_dsv_ecoli_long, "\n")

  dsv_ecoli_data <- read_csv(
    file_dsv_ecoli_long,
    show_col_types = FALSE
  )

  required_cols <- c(
    "SampleID",
    "Group",
    "Country",
    "StudyID",
    "Escherichia_coli",
    "Desulfovibrio"
  )

  missing_cols <- setdiff(required_cols, colnames(dsv_ecoli_data))

  if (length(missing_cols) > 0) {
    stop(
      "Combined DSV/E. coli table missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nDetected columns: ",
      paste(colnames(dsv_ecoli_data), collapse = ", ")
    )
  }

  maaslin_metadata <- dsv_ecoli_data %>%
    mutate(
      SampleID = as.character(SampleID),
      Group = factor(Group, levels = c("Control", "PD")),
      Country = as.factor(Country),
      StudyID = as.factor(StudyID)
    ) %>%
    select(SampleID, Group, Country, StudyID) %>%
    distinct(SampleID, .keep_all = TRUE)

  maaslin_features <- dsv_ecoli_data %>%
    mutate(
      SampleID = as.character(SampleID),
      Escherichia_coli = as.numeric(Escherichia_coli),
      Desulfovibrio = as.numeric(Desulfovibrio)
    ) %>%
    select(SampleID, Escherichia_coli, Desulfovibrio) %>%
    distinct(SampleID, .keep_all = TRUE) %>%
    pivot_longer(
      cols = c(Escherichia_coli, Desulfovibrio),
      names_to = "Feature",
      values_to = "Abundance"
    ) %>%
    pivot_wider(
      names_from = SampleID,
      values_from = Abundance,
      values_fill = list(Abundance = 0)
    )

  output_maaslin_metadata <- file.path(
    dir_maaslin,
    "MaAsLin3_Metadata.csv"
  )

  output_maaslin_features <- file.path(
    dir_maaslin,
    "MaAsLin3_Features_DSV_Ecoli.csv"
  )

  write_csv(
    maaslin_metadata,
    output_maaslin_metadata
  )

  write_csv(
    maaslin_features,
    output_maaslin_features
  )

  cat("MaAsLin3 metadata saved to:", output_maaslin_metadata, "\n")
  cat("MaAsLin3 feature table saved to:", output_maaslin_features, "\n")

} else {

  warning(
    "Combined DSV/E. coli long-format table was not found. ",
    "Skipping MaAsLin3 input conversion: ",
    file_dsv_ecoli_long
  )
}

# ==========================================================
# 8. Export preprocessing summary
# ==========================================================

summary_table <- tibble(
  Metric = c(
    "Merged MetaPhlAn taxa/features",
    "Sample columns in MetaPhlAn table",
    "Genus-level features",
    "Desulfovibrio rows",
    "Species-level features",
    "Escherichia_coli rows"
  ),
  Value = c(
    nrow(abundance_data),
    ncol(abundance_data) - 1,
    nrow(genus_data),
    nrow(desulfovibrio_data),
    nrow(species_data),
    nrow(ecoli_data)
  )
)

output_summary <- file.path(
  dir_target,
  "Taxonomic_Feature_Preprocessing_Summary.csv"
)

write_csv(
  summary_table,
  output_summary
)

print(summary_table)

writeLines(
  capture.output(sessionInfo()),
  file.path(dir_target, "sessionInfo_taxonomic_preprocessing.txt")
)

cat("Taxonomic feature preprocessing completed successfully.\n")
cat("Summary saved to:", output_summary, "\n")
