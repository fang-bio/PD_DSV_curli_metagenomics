# Upstream processing

This folder contains representative shell and SLURM scripts used for upstream metagenomic processing in the Parkinson's disease gut metagenome meta-analysis.

The full dataset was processed on the Puhti high-performance computing environment. For scheduling and resource-management purposes, raw-data download, quality control, taxonomic profiling, functional profiling, and optional assembly were performed in multiple cohort-specific or batch-specific jobs. The scripts provided here are representative examples of the actual workflow. Across batches, the same software environments, database paths, and processing parameters were used; only cohort identifiers or sample ID ranges differed.

Raw shotgun metagenomic reads were obtained from publicly available repositories listed in Supplementary Table 1 of the manuscript. Raw sequencing files are not included in this repository because of their large file size and public availability.

## Scripts

- `01_raw_fastq_download_representative_cohort.slurm`: Representative cohort-specific SLURM script for downloading raw FASTQ files from public sequencing repositories.
- `02_kneaddata_qc_dehost_representative_batch.slurm`: Representative KneadData batch script for read-level quality control and host-read removal.
- `03_metaphlan_taxonomic_profiling_representative_batch.slurm`: Representative MetaPhlAn batch script for taxonomic profiling of KneadData-filtered shotgun metagenomic reads.
- `04_humann_functional_profiling_representative_batch.slurm`: Representative HUMAnN batch script for functional profiling of KneadData-filtered shotgun metagenomic reads.
- `05_megahit_assembly_representative_batch.slurm`: Representative MEGAHIT batch script for per-sample metagenomic assembly. This assembly step documents an exploratory upstream workflow and was not used in the final statistical analyses reported in the manuscript.

File paths and SLURM settings reflect the Puhti computing environment used in this study and may need to be adapted before reuse in other computing environments.
