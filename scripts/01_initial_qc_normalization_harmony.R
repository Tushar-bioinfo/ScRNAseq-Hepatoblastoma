# Script: 01_initial_qc_normalization_harmony.R
# Purpose: Load and merge 10X data, perform QC filtering, normalization,
#          PCA, and batch correction using Harmony.

# Load required libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(harmony)
library(dplyr)

# Base directory containing all 10X sample folders
basedir <- "/Users/tusharsingh/Downloads/single_RNAseq/workflow/hb"
dirs <- list.dirs(basedir, full.names = FALSE, recursive = FALSE)

# Initialize list to store Seurat objects
seurat_list <- list()

# Step 1: Create Seurat objects for each sample
for (i in dirs) {
  read <- Read10X(file.path(basedir, i))
  seurat_object <- CreateSeuratObject(counts = read, min.cells = 3, min.features = 200)
  name <- gsub("_filtered_feature_bc_matrix", "", i)
  seurat_list[[name]] <- seurat_object
}
rm(read)

# Step 2: Merge all samples into one object
merged_seurat_obj <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "HB"
)
rm(seurat_object)
rm(seurat_list)

# Step 3: Annotate metadata
merged_seurat_obj$sample <- rownames(merged_seurat_obj@meta.data)
merged_seurat_obj$sample_id <- sapply(merged_seurat_obj$sample, function(x) strsplit(x, "_")[[1]][1])
merged_seurat_obj$sample_type <- sapply(merged_seurat_obj$sample, function(x) strsplit(x, "_")[[1]][2])
merged_seurat_obj$sample <- NULL

# Step 4: Calculate % mitochondrial genes
merged_seurat_obj$mt_percent <- PercentageFeatureSet(merged_seurat_obj, pattern = "^MT-")

# Step 5: QC Plots
VlnPlot(merged_seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "mt_percent"), ncol = 3)
FeatureScatter(merged_seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")

# Step 6: Filter low-quality cells
filtered_mg_seurat_obj <- subset(
  merged_seurat_obj,
  subset = nFeature_RNA > 500 & nCount_RNA > 1000 & mt_percent < 10
)
rm(merged_seurat_obj)

# Step 7: Normalize and identify variable features
filtered_mg_seurat_obj <- NormalizeData(filtered_mg_seurat_obj)
filtered_mg_seurat_obj <- FindVariableFeatures(filtered_mg_seurat_obj, selection.method = "vst", nfeatures = 2000)

# Step 8: Scale and PCA
filtered_mg_seurat_obj <- ScaleData(filtered_mg_seurat_obj)
filtered_mg_seurat_obj <- RunPCA(filtered_mg_seurat_obj)

# Step 9: Elbow Plot
ElbowPlot(filtered_mg_seurat_obj)

# Step 10: Visualize UMAP before integration
before1 <- DimPlot(filtered_mg_seurat_obj, reduction = "umap", group.by = "sample_id", label = FALSE,
                   cols = c("red", "green", "blue"))
before2 <- DimPlot(filtered_mg_seurat_obj, reduction = "umap", group.by = "sample_type", label = FALSE,
                   cols = c("red", "green", "blue"))
before1 + before2

# Step 11: Harmony Integration
harmony_obj <- filtered_mg_seurat_obj %>%
  RunHarmony(group.by.vars = "sample_id", plot_convergence = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.1)

# Step 12: Set clustering identity
Idents(harmony_obj) <- "RNA_snn_res.0.1"

# Step 13: Visualize UMAP post-Harmony
after <- DimPlot(harmony_obj, reduction = "umap", group.by = "sample_id")
after2 <- DimPlot(harmony_obj, reduction = "umap", group.by = "sample_type")
after + after2
after2

# Step 14: Optional cluster visualization
clusters <- DimPlot(harmony_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
clusters

# Step 15: Check sample composition by cluster
table(Idents(harmony_obj), harmony_obj$sample_type)

# Step 16: Subset for downstream DE analysis (e.g., fetal-like HB)
subset_obj <- subset(harmony_obj, subset = sample_type %in% c("PDX", "tumor"))
subset_obj <- JoinLayers(subset_obj)
