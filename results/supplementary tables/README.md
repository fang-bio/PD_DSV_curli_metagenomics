# Supplementary data

This folder contains the supplementary data associated with the multicohort analysis of fecal metagenomic profiles in Parkinson’s disease.

The Excel workbook `supplementary_tables.xlsx` contains dataset and sample-level information, quality-control summaries, targeted microbial-feature analyses, ecological correlation analyses, mixed-effects regression outputs, predicted probabilities, and inferred pathway-level functional profiles supporting the main and supplementary results.

## File

### `supplementary_tables.xlsx`

Excel workbook containing Supplementary Tables 1–15.

## Table descriptions

| Table                  | Description                                                                                                                                                                                                       |
| ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Supplementary Table 1  | Summary of publicly available fecal metagenomic datasets included in this study                                                                                                                                   |
| Supplementary Table 2  | Sample-level manifest of all downloaded metagenomic runs, including accession identifiers, cohort and disease metadata, sequencing and quality-control metrics, final inclusion status, and reasons for exclusion |
| Supplementary Table 3  | Sample-level metadata for metagenomic samples retained after quality control and host-read removal and included in the final analyses                                                                             |
| Supplementary Table 4  | Cohort-level sample accounting, quality-control reconciliation, sequencing characteristics, and read-depth summary for the downloaded and finally included metagenomic datasets                                   |
| Supplementary Table 5  | PERMANOVA assessment of the associations of disease status, cohort, geographic region, and sequencing platform with gut microbial community structure                                                             |
| Supplementary Table 6  | Targeted MaAsLin3 analysis of *Desulfovibrio*, *Escherichia coli*, and curli-associated genes in Parkinson’s disease                                                                                              |
| Supplementary Table 7  | Disease-stratified Spearman correlation analysis of ecological relationships among *Desulfovibrio*, *Escherichia coli*, and curli-associated genes                                                                |
| Supplementary Table 8  | Mixed-effects logistic regression models assessing continuous interaction terms between *Desulfovibrio*, curli-associated genes, and total *Escherichia coli* abundance                                           |
| Supplementary Table 9  | Predicted Parkinson’s disease probabilities from mixed-effects logistic regression models across standardized *Desulfovibrio* and curli-associated gene abundance gradients                                       |
| Supplementary Table 10 | Odds ratios for Parkinson’s disease case status across categorical *Desulfovibrio*–curli exposure groups based on within-cohort median stratification                                                             |
| Supplementary Table 11 | PERMANOVA assessment of inferred metagenomic functional-profile variation associated with *Desulfovibrio*–curli co-elevation and cohort effects                                                                   |
| Supplementary Table 12 | Group-mean standardized abundances of pathway-level features across *Desulfovibrio*–curli exposure groups                                                                                                         |
| Supplementary Table 13 | Significant pathway-level features associated with the *Desulfovibrio* × curli-associated gene interaction term                                                                                                   |
| Supplementary Table 14 | Significant pathway-level features with positive *Desulfovibrio* × curli-associated gene interaction coefficients                                                                                                 |
| Supplementary Table 15 | Significant pathway-level features with negative *Desulfovibrio* × curli-associated gene interaction coefficients                                                                                                 |

## Notes

The supplementary tables were generated or assembled from the analysis scripts and processed outputs provided in this repository. They support reproducibility of the dataset accounting, quality-control procedures, targeted microbial-feature analyses, ecological correlation analyses, mixed-effects regression models, categorical joint-association analyses, and inferred pathway-level functional analyses reported in the manuscript.

The curli-associated gene measure was derived from community-level unstratified HUMAnN gene-family profiles and therefore represents the summed abundance of gene families corresponding to the selected curli-associated UniRef90 annotations rather than exclusively *Escherichia coli*-derived curli gene content.

Raw metagenomic sequencing data are available from the original public repositories under the project accession numbers listed in Supplementary Table 1.
