# ZER1 Interactome and CRL2 Complex Analysis

Analysis code for manuscript investigating ZER1 protein interactions and CRL2 complex expression in cardiac hypertrophy and cardiomyopathy.

## Repository Contents

| File | Description |
|------|-------------|
| `ZER1_interactome_tier_analysis.R` | IP-MS data integration, pathway annotation, and tier classification |
| `CRL2_expression_analysis.R` | CRL2 complex expression analysis (human DCM & mouse TAC) |
| `nterm_cleavage_analysis.py` | N-terminal methionine cleavage prediction |

## Requirements

### R (≥ 4.2.1)
```r
# CRAN packages
install.packages(c("readxl", "dplyr", "pheatmap", "ggplot2", "ggrepel", 
                   "reshape2", "openxlsx", "readr", "tidyr", "ggpubr", 
                   "stringr", "tibble"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db", "GO.db"))
