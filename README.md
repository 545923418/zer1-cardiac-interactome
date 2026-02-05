# ZER1 Interactome and CRL2 Complex Analysis

Analysis code for manuscript investigating ZER1 protein interactions and CRL2 complex expression in cardiac hypertrophy and cardiomyopathy.

## Repository Contents

| File | Description | Output |
| :--- | :--- | :--- |
| `CRL2_expression_analysis.R` | CRL2 complex expression analysis (human DCM & mouse TAC) | **Fig 1** (A-B), **Supp Fig I** (A-C) |
| `ZER1_interactome_tier_analysis.R` | IP-MS data integration, pathway annotation, and tier classification | **Supp Fig III**, **Supp Table 1-2** |
| `nterm_cleavage_analysis.py` | N-terminal methionine cleavage prediction | Sequence analysis results |

---

## Analysis Workflow

### 1. Bulk Module Score Analysis
Evaluation of the CRL2 complex expression signatures across different cardiac disease states.
* **Input Data:** * `GSE116250`: Human Dilated Cardiomyopathy (DCM).
    * `GSE203083`: Mouse Transverse Aortic Constriction (TAC) model.
* **Outputs:** Figure 1 panel A-B; Supplementary Figure I panel A-C.

### 2. IP–MS Prioritization
Filtering and ranking of ZER1-interacting proteins to identify high-confidence candidates.
* **Input Data:** IP–MS result table.
* **Outputs:** Supplementary Figure III; Supplementary table 1-2.

---

## Requirements

### R (≥ 4.2.1)
The following packages are required for data processing and visualization:

```R
# CRAN packages
install.packages(c("readxl", "dplyr", "pheatmap", "ggplot2", "ggrepel", 
                   "reshape2", "openxlsx", "readr", "tidyr", "ggpubr", 
                   "stringr", "tibble"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db", "GO.db"))
