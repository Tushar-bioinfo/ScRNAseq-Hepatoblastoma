# Single-Nucleus-RNA-Seq-Analysis-of-Hepatoblastoma-Tumors
This repository contains the code, data processing pipeline, and visualizations for a single-nucleus RNA-seq (snRNA-seq) study of hepatoblastoma.The dataset includes nuclei from patient tumor tissues, background liver, and patient-derived xenografts (PDX). The project replicates and extends analyses from the hepatoblastoma study published in Nature Cancer (PMID: 34497364).

Key Steps:
Preprocessing and quality control of snRNA-seq data using Seurat.

Batch correction and clustering with Harmony, yielding 10 distinct clusters (C0–C9).

Manual annotation of clusters based on canonical marker genes using GeneCards and literature references.

Differential expression analysis between tumor and PDX cells in specific clusters (e.g., fetal-like hepatoblastoma subpopulations).

Visualization of gene expression (e.g., AFP, CFH) and cell-type distributions.

Quantification of immune infiltration (e.g., T cell proportions) across sample types.

Tools & Packages:
Seurat, Harmony, ggplot2, dplyr

Feature plots, violin plots, volcano plots, and cell composition bar plots

Goals:
Characterize intratumoral heterogeneity and fetal-like expression patterns

Compare tumor vs PDX gene expression profiles


