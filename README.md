# PD_DSV_curli_metagenomics

This repository contains analysis code and processed results associated with the manuscript:

**Multi-cohort metagenomic analysis identifies joint enrichment of Desulfovibrio and curli-associated genes in Parkinson’s disease**

This study presents a multicohort reanalysis of fecal shotgun metagenomic data from 1,609 samples across 10 Parkinson’s disease cohorts. It evaluates global microbiome structure and targeted associations of *Desulfovibrio*, *Escherichia coli*, and curli-associated genes with Parkinson’s disease case status.

Detection-prevalence and abundance associations were evaluated separately using cohort-aware mixed-effects models. The study also examined ecological correlations, continuous and categorical joint associations of *Desulfovibrio* and curli-associated genes, and variation in inferred metagenomic functional profiles.

The curli-associated gene measure was derived from community-level unstratified HUMAnN gene-family profiles using selected UniRef90 annotations corresponding to the *Escherichia coli csg* curli system. It therefore represents community-level curli-associated gene content rather than exclusively *E. coli*-derived curli gene abundance.

## Repository contents

- `R_scripts/`  
  R scripts used for feature preprocessing, statistical analysis, and figure generation. These include preparation of MetaPhlAn taxonomic profiles and HUMAnN functional profiles, construction of curli-associated gene features, MaAsLin3 mixed-effects analyses, Spearman correlation analyses, mixed-effects logistic regression, PERMANOVA, pathway-level interaction models, and generation of the main and supplementary figures.

- `upstream_processing/`  
  Representative shell and SLURM scripts used for raw-data download, quality control, host-read removal, taxonomic profiling with MetaPhlAn, and functional profiling with HUMAnN.

- `results/`  
  Processed results associated with the manuscript, including main figures, supplementary figures, supplementary data, and statistical output files.

- `results/supplementary tables/`  
  The Excel workbook containing Supplementary Tables 1–15, including dataset and sample-level information, quality-control summaries, targeted microbial-feature analyses, ecological correlations, mixed-effects regression outputs, predicted probabilities, and pathway-level functional results.

## Data availability

Raw shotgun metagenomic sequencing data are not included in this repository because of their large file size and availability through public sequencing repositories. The project accession numbers for all included datasets are provided in the manuscript and Supplementary Table 1.

Harmonized metadata, processed feature tables, statistical outputs, and analysis scripts supporting the reported results are available within this repository.

## Reproducibility

The scripts and processed files are organized to support reproduction of the statistical analyses and figures reported in the manuscript. File paths and computational settings may require adjustment for use in a different computing environment.
