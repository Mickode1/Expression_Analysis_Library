# Title

**Differential Gene Expression Analysis: Normal vs. UUO Mouse Models**

## Background

**Unilateral Ureteral Obstruction (UUO)** is a commonly used experimental model for studying renal injury. While it effectively models the kidney's response to urinary tract blockage, its applications extend beyond this. UUO is an excellent model for investigating irreversible acute kidney disease and the mechanisms underlying chronic kidney diseases [1](https://link.springer.com/article/10.1007/s11255-013-0520-1#citeas)

Therefore, this project explores and uncovers the genes significantly expressed in the UUO mouse model and the associated biological processes. The goal is to understand the molecular mechanisms involved in these disease states and to identify potential therapeutic interventions.

This repository is divided into these sections:

- Project Staging: Details project dependencies, data loading and preparation, quality control, filtering, and normalization steps.
- Analysis: The analysis includes differential expression analysis and functional enrichment analysis.
- Results and Conclusion: The biological implications of the analyses were discussed.

## Data Sources and Analyses
The datasets used in this project were obtained from an experimental study titled [Silencing SMOC2 ameliorates kidney fibrosis by inhibiting fibroblast to myofibroblast transformation](https://insight.jci.org/articles/view/90299). The datasets are publicly available under the accession number [GSE85209](https://www.ncbi.nlm.nih.gov/gds/?term=GSE85209). The sequencing data for the cell lines can be accessed through the following accession numbers: GSM2260466 GSM2260467 GSM2260468 GSM2260469 GSM2260470 GSM2260471 GSM2260472

In this project, I will explore the genes differentially expressed in the Unilateral Ureteral Obstruction (UUO) mouse model. This study aims to understand the molecular functions of the differentially expressed genes and their associated biological processes.

The following analyses will be performed:

- **Differential Expression Analysis (DEA)**: This analysis will focus on identifying the genes that are significantly expressed in the UUO mouse model and the control model. Four replicates of the UUO mice and three replicates of the healthy mouse models will be analyzed.

- **Functional Enrichment Analysis (FEA)**: This analysis will determine the biological significance of the differentially expressed genes. By categorizing these genes into known biological functions and pathways, FEA will help uncover the molecular mechanisms and pathways involved in the disease state. This can provide insights into the processes involved and potentially lead to the identification of therapeutic targets for intervention.

  # Project Staging

  This section is broken into three sub-sections:

- Project Dependencies: Import the necessary libraries for downstream analyses.
- Import and Prepare Data: Import raw datasets and prepare them for downstream processing.
- Quality Control, Filtering, and Normalization: Applying quality control measures and data transformations to generate high-quality, reliable data for downstream analysis.

## Project dependencies

To start, we will need to import the R packages needed for the project

```r
library(dplyr)
library(tidyverse)
library(pheatmap)
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(pathview)
library(gage)
library(annotables)

## Import and Prepare Data






