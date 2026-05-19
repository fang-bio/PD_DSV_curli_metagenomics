# Preprocessing, analysis, and figure-generation scripts

This folder contains scripts used for upstream feature preprocessing, statistical analysis, and figure generation for the Parkinson’s disease metagenomic meta-analysis.

The workflow focuses on three major microbial features: **Desulfovibrio**, **Escherichia coli**, and **E. coli curli-gene abundance**. Taxonomic features were derived from MetaPhlAn profiles, whereas curli-gene abundance was derived from HUMAnN genefamilies output using curated UniRef90 identifiers corresponding to E. coli csgA–csgG genes.

---

## Upstream feature preprocessing

- `00_preprocess_metaphlan_target_taxa.R`  
  Preprocesses targeted taxonomic features from merged MetaPhlAn profiles. The script extracts and collapses genus-level profiles for **Desulfovibrio** and species-level profiles for **Escherichia coli**, and prepares target taxonomic abundance tables for downstream statistical analyses and visualization.

- `00a_convert_humann_genefamilies_rpk_to_relative_abundance.slurm`  
  Converts HUMAnN genefamilies RPK profiles to sample-wise relative abundance values before targeted curli-gene extraction. The script excludes unmapped and unintegrated features when calculating total mapped RPK values and generates relative-abundance gene-family profiles.

- `00b_extract_ecoli_curli_genes_from_relative_abundance.slurm`  
  Extracts E. coli curli-related UniRef90 gene families from HUMAnN gene-family relative-abundance profiles. The script uses curated UniRef90 identifiers corresponding to E. coli csgA–csgG genes, performs exact ID matching across samples, and exports sample-gene-level curli abundance for downstream targeted MaAsLin3 and Desulfovibrio–curli interaction analyses.

- `00_extract_ecoli_curli_genes_from_humann.slurm`  
  Provides the earlier HUMAnN-based curli-gene extraction workflow used during feature development. This script documents the upstream extraction of E. coli curli-related UniRef90 gene families from HUMAnN genefamilies output. The final two-step workflow is represented by `00a_convert_humann_genefamilies_rpk_to_relative_abundance.slurm` and `00b_extract_ecoli_curli_genes_from_relative_abundance.slurm`.

---

## Figure 1 scripts

- `01b_global_microbiome_structure_pcoa_permanova.R`  
  Generates Figure 1b. The script evaluates global species-level microbiome structure using CLR-transformed relative abundance profiles and Aitchison distance. It performs PCoA ordination and tests the disease-status effect using cohort-stratified PERMANOVA with 9,999 permutations. The script exports the Figure 1b PCoA plot, PCoA coordinates, Aitchison distance matrix, PERMANOVA result tables, and source data.

- `01c_targeted_maaslin3_desulfovibrio_ecoli_curli.R`  
  Generates Figure 1c. The script performs targeted MaAsLin3 mixed-effects modeling for Desulfovibrio, Escherichia coli, and bacterial curli-gene abundance, testing differences between Parkinson’s disease patients and controls while accounting for cohort-level heterogeneity.

- `01d_classical_abundance_comparison.R`  
  Generates Figure 1d. The script performs classical abundance comparisons for the targeted microbial features and provides a complementary non-model-based visualization of differences between Parkinson’s disease and control groups.

---

## Figure 2 scripts

- `02_step0_preprocess_figure2_features.R`  
  Preprocesses metadata and targeted microbial features for Figure 2 and Supplementary Figure 2 analyses. The script merges sample-level metadata with Desulfovibrio, Escherichia coli, and curli-gene abundance profiles; generates detection indicators; applies feature-specific pseudocounts; and produces log-transformed standardized variables for downstream correlation and interaction analyses.

- `02ab_dsv_curli_ecoli_correlation.R`  
  Generates Figure 2a–b. The script evaluates disease-stratified Spearman correlations between Desulfovibrio abundance and curli-gene abundance among carrier samples, and compares this relationship with the association between Desulfovibrio and total Escherichia coli abundance.

- `02cd_dsv_curli_interaction_risk_landscape.R`  
  Generates Figure 2c–d. The script fits mixed-effects logistic regression models to estimate the interaction effect of Desulfovibrio and curli-gene abundance on Parkinson’s disease status while accounting for cohort-level heterogeneity. It compares the target Desulfovibrio × curli-gene interaction with a Desulfovibrio × total Escherichia coli comparator model and visualizes the predicted PD-risk landscape for the Desulfovibrio × curli model.

- `02_suppfig2a_dsv_ecoli_comparator_landscape.R`  
  Generates Supplementary Figure 2a, which supports the specificity analysis of Figure 2c–d. The script fits mixed-effects logistic regression models for the target Desulfovibrio × curli-gene interaction and the comparator Desulfovibrio × total Escherichia coli interaction, with cohort identity included as a random effect. Predicted Parkinson’s disease probabilities are visualized across standardized abundance gradients, and the Desulfovibrio × total Escherichia coli landscape is used as the comparator output for Supplementary Figure 2a.

---

## Figure 3 scripts

- `03_step0_preprocess_functional_pathways.R`  
  Preprocesses pathway-level functional profiles for Figure 3 analyses. The script integrates sample metadata with Desulfovibrio and curli-gene abundance, defines cohort-specific DSV–curli combined exposure groups, retains unstratified HUMAnN pathway profiles, removes unmapped and unintegrated features, aligns pathway profiles with metadata, and applies a 10% pathway prevalence filter before downstream functional profiling.

- `03a_global_functional_reorganization_pcoa_permanova.R`  
  Generates Figure 3a. The script applies CLR transformation to pathway-level HUMAnN profiles, computes Aitchison distances, performs PCoA, and uses PERMANOVA with 9,999 permutations to evaluate cohort, Desulfovibrio, curli-gene, and interaction effects on global functional pathway profiles. A cohort-colored ordination and a variance-explained plot are also generated as supplementary or diagnostic outputs. The permutation-based PERMANOVA step was run on the Puhti high-performance computing environment because of its computational cost.

- `03b_dsv_curli_functional_shift_lollipop.R`  
  Generates Figure 3b. The script visualizes the top pathway-level functional shifts associated with the Desulfovibrio–curli interaction. It reads the differential pathway association results, retains pathways passing the q-value threshold, selects the top positively and negatively associated pathways based on absolute effect size, generates a diverging lollipop plot, and exports Figure 3b source data and summary statistics.

- `03cd_prepare_synergy_pathway_heatmap_matrices.R`  
  Prepares the source matrices for Figure 3c–d heatmaps. The script reads the pathway-level Desulfovibrio–curli interaction results, selects the top positively and negatively associated pathways passing the q-value threshold, aligns these pathways with the HUMAnN pathway abundance matrix and sample metadata, applies log transformation and row-wise z-score scaling, and exports group-averaged heatmap matrices across DSV–curli combined exposure groups.

- `03cd_plot_synergy_pathway_heatmaps.R`  
  Generates Figure 3c–d heatmaps from the group-averaged pathway matrices prepared in the previous step. The script visualizes positively and negatively associated pathways across DSV–curli combined exposure groups, exports the final heatmaps, and saves the corresponding source data and plotting statistics.

---

## Workflow summary

The recommended script order is:

```text
00_preprocess_metaphlan_target_taxa.R
00a_convert_humann_genefamilies_rpk_to_relative_abundance.slurm
00b_extract_ecoli_curli_genes_from_relative_abundance.slurm

01b_global_microbiome_structure_pcoa_permanova.R
01c_targeted_maaslin3_desulfovibrio_ecoli_curli.R
01d_classical_abundance_comparison.R

02_step0_preprocess_figure2_features.R
02ab_dsv_curli_ecoli_correlation.R
02cd_dsv_curli_interaction_risk_landscape.R
02_suppfig2a_dsv_ecoli_comparator_landscape.R

03_step0_preprocess_functional_pathways.R
03a_global_functional_reorganization_pcoa_permanova.R
03b_dsv_curli_functional_shift_lollipop.R
03cd_prepare_synergy_pathway_heatmap_matrices.R
03cd_plot_synergy_pathway_heatmaps.R
