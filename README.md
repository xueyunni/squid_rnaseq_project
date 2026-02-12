# Squid RNA-seq Project
# Doryteuthis pealeii – Statocyst vs Gill

## Overview

This repository contains code and documentation for a bulk RNA-seq analysis comparing statocyst and gill tissues from the squid Doryteuthis pealeii.
The goal of this project is to identify candidate mechanotransduction genes in squid statocyst hair cells.

This repository is designed to be:
- reproducible
- script-driven
- safe for version control (no large raw data tracked)

## Project structure 

```
squid_rnaseq_project/
|-- data/ 
|-- docs/
|-- README.md
|-- renv/
|-- renv.lock
|-- results/
|-- sample_sheet.r            
|--scripts/
 |--00_setup_packages.R
 |--01_deseq2_analysis.R
 |--02_enrichment_GO
 |--03_GSEA_GO.R
 |--preprocessing
|--squid_bulk.Rproj
               
```

## Data and Results

Sequencing data (FASTQ), QC reports, trimming results, alignment files (BAM), and raw read counts files, tables and plots are not tracked in this repository.
These files are stored separately on local disks or HPC clusters and are excluded via .gitignore.

Only scripts and sample metadata are tracked here.

## Analysis workflow

1. Quality control (fastqc) and trimming (fastp)
2. Aligment (STAR)
- Mapping reads to the squid genome
3. Differential expression analysis
- Statocyst vs gill comparison
- Idenfication of significantly enriched genes
4. Pathway enrichment analysis
- GO enrichment analysis
- GSEA using ranked gene lists
 
## Reproducibility

This project uses renv to ensure reproducible R environment

To store the exact package version used in the analysis:
renv::restoreO()

R version and package versions are recorded in renv.lock

## Requirements

- R (≥ 4.x)
- renv
- Common R packages
 - DESeq2
 - clusterProfiler
 - ggplot2
 - tidyverse
 
 ## Author
 
 Max Ni 
 Biology Master's student 
 Case Western Reserve University




	




