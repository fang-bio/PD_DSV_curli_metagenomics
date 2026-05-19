# PD_DSV_curli_meta-analysis_metagenomics

This repository contains code and processed results for the manuscript:

**Multi-cohort metagenomic analysis identifies prevalence-driven *Desulfovibrio*–curli co-enrichment and functional reorganization in Parkinson’s disease**

This study performs a multi-cohort fecal metagenomic meta-analysis of Parkinson’s disease, focusing on prevalence-driven microbial signals, ecological interactions, and functional reorganization associated with *Desulfovibrio*, *Escherichia coli*, and *E. coli* curli-gene abundance.

## Repository contents

- `R_scripts/`  
  Scripts used for upstream feature preprocessing, statistical analysis, and figure generation. This folder includes MetaPhlAn-based taxonomic preprocessing, HUMAnN-based *E. coli* curli-gene feature construction, MaAsLin3 mixed-effects analyses, PERMANOVA, logistic regression models, and figure-generation scripts.

- `upstream_processing/`  
  Representative shell and SLURM scripts used for upstream metagenomic processing, including raw-data download, quality control, host-read removal, taxonomic profiling, functional profiling, and optional assembly-related steps.

- `results/`  
  Processed results associated with the manuscript, including main figures, supplementary figures, supplementary tables, and source-data outputs.

## Notes

Raw sequencing data are not included in this repository because of their large file size and public availability through public sequencing repositories. The repository provides processed result files and scripts required to reproduce the statistical analyses and figures reported in the manuscript.
