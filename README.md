# Glioblastoma Multiforme (GBM) RNA-Seq Analysis

## Overview

This repository contains R scripts for analyzing Glioblastoma Multiforme (GBM) RNA-Seq data from The Cancer Genome Atlas (TCGA). The project aims to identify differentially expressed genes between GBM samples and perform pathway enrichment analysis to understand the underlying biological processes.

## Key Features

- Data acquisition from TCGA
- Data preprocessing and normalization
- Differential expression analysis
- Pathway enrichment analysis (GO and KEGG)
- Visualization of results
- Export of analysis results

## Prerequisites

Ensure you have R installed along with the following packages:

- BiocManager
- TCGAbiolinks
- SummarizedExperiment
- edgeR
- limma
- ggplot2
- pheatmap
- org.Hs.eg.db
- clusterProfiler

Install the required packages using:


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "edgeR", "limma",
                       "ggplot2", "pheatmap", "org.Hs.eg.db", "clusterProfiler"))


## Usage
Run the main script:


gbm_analysis.R: https://github.com/Derleen123/Glioblastoma-Multiforme-GBM-RNA-Seq-data-analysis/blob/main/gbm_analysis.R


## Workflow

1. **Data Acquisition**: Query and download GBM RNA-Seq data from TCGA.
2. **Data Preprocessing**: Normalize data and filter low-expression genes.
3. **Differential Expression Analysis**: Identify significantly differentially expressed genes (DEGs).
4. **Pathway Enrichment Analysis**: Perform GO and KEGG enrichment analysis on DEGs.
5. **Visualization**: Generate heatmaps and dot plots for enriched pathways.
6. **Results Export**: Save analysis results and visualizations.

## Output

- `heatmap_top_genes.png`: Heatmap of top differentially expressed genes
- `GO_dotplot.pdf`: Dot plot of top GO terms
- `KEGG_dotplot.pdf`: Dot plot of top KEGG pathways
- `differential_expression_results.csv`: Full list of differentially expressed genes
- `go_enrichment_results.csv`: GO enrichment analysis results
- `kegg_pathway_enrichment_results.csv`: KEGG pathway enrichment analysis results

## Acknowledgments

This project is part of the Advanced Genomics Course: Bioinformatics for Cancer Biology from HackBio. For more information, visit [HackBio Course](https://course.thehackbio.com/).

Data provided by [The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/ccg/research/genome-sequencing/tcga).

