# Script: 02_annotation_and_differential_expression.R
# Purpose: Perform manual cell type annotation and differential expression (DE)
#          analysis between tumor and PDX cells within the fetal-like HB cluster.

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Step 1: Set identity to clusters
Idents(harmony_obj) <- "seurat_clusters"

# Step 2: Join layers to enable DE testing
harmony_obj <- JoinLayers(harmony_obj)

# Step 3: Find markers for all clusters
markers <- FindAllMarkers(
  object = harmony_obj,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

# Step 4: Select top 10 markers per cluster
top10 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
head(top10)

# Step 5: Plot heatmap for top markers
DoHeatmap(harmony_obj, features = top10$gene) + NoLegend()

# Step 6: Manual annotation of clusters based on known marker genes
harmony_obj$celltype <- "unknown"
harmony_obj$celltype[which(Idents(harmony_obj) == "0")] <- "Tumor-like (stem)"
harmony_obj$celltype[which(Idents(harmony_obj) == "1")] <- "Fetal-like hepatoblastoma"
harmony_obj$celltype[which(Idents(harmony_obj) == "2")] <- "Embryonic tumor-like"
harmony_obj$celltype[which(Idents(harmony_obj) == "3")] <- "Bulk hepatoblastoma"
harmony_obj$celltype[which(Idents(harmony_obj) == "4")] <- "Normal hepatocytes"
harmony_obj$celltype[which(Idents(harmony_obj) == "5")] <- "T cells"
harmony_obj$celltype[which(Idents(harmony_obj) == "6")] <- "Liver endothelial cells"
harmony_obj$celltype[which(Idents(harmony_obj) == "7")] <- "Proliferating tumor cells"
harmony_obj$celltype[which(Idents(harmony_obj) == "8")] <- "Lymphatic cells"
harmony_obj$celltype[which(Idents(harmony_obj) == "9")] <- "Stromal / smooth muscle"

# Step 7: Visualize annotated clusters
clusters <- DimPlot(harmony_obj, reduction = "umap", group.by = "celltype", label = TRUE)
clusters

# Step 8: Feature expression plots
FeaturePlot(harmony_obj, features = c("RTN4RL1", "NPSR1", "MLIP-AS1", "RERGL"), reduction = "umap")
VlnPlot(harmony_obj, features = c("RTN4RL1", "NPSR1", "MLIP-AS1", "RERGL"))

FeaturePlot(harmony_obj, features = c("LYZ", "CD68", "CD163", "HBB", "CD302"), reduction = "umap")

# Step 9: Differential expression in fetal-like cluster (cluster 1)
fetal_like <- subset(harmony_obj, subset = seurat_clusters == 1)
table(fetal_like$sample_type)

deg_fetal <- FindMarkers(
  object = fetal_like,
  ident.1 = "tumor",
  ident.2 = "PDX",
  min.pct = 0.1,
  logfc.threshold = 0.25
)

head(deg_fetal, 20)

# Step 10: Plot selected DE genes
FeaturePlot(fetal_like, features = c("CFH", "C5", "GREM2", "AFP"), split.by = "sample_type")
VlnPlot(fetal_like, features = c("CYP3A5", "KCNT2", "CFH", "FAM13A", "C5"), group.by = "sample_type")

# Optional: Subset to tumor and PDX
fetal_like_tumor_pdx <- subset(fetal_like, subset = sample_type %in% c("tumor", "PDX"))
VlnPlot(fetal_like_tumor_pdx, features = c("CFH", "C5", "CYP3A5", "AFP"), group.by = "sample_type")
FeaturePlot(fetal_like_tumor_pdx, features = c("CFH", "CYP3A5", "AFP"), split.by = "sample_type", reduction = "umap")

# Optional housekeeping genes
VlnPlot(fetal_like_tumor_pdx, features = c("GAPDH", "RPLP0", "ACTB"), split.by = "sample_type", pt.size = 0.1)

# Step 11: Heatmap of top DE genes
top_genes <- rownames(head(deg_fetal[order(deg_fetal$p_val_adj), ], 20))
DoHeatmap(fetal_like_tumor_pdx, features = top_genes, group.by = "sample_type") + NoLegend()

# Step 12: Optional downsampling for balanced DE
pdxcells <- WhichCells(fetal_like, expression = sample_type == "PDX")
tumorcells <- WhichCells(fetal_like, expression = sample_type == "tumor")

set.seed(123)
tumor_sampled <- sample(tumorcells, length(pdxcells))
balanced_cells <- c(pdxcells, tumor_sampled)

balanced_obj <- subset(fetal_like, cells = balanced_cells)
deg_balanced <- FindMarkers(
  object = balanced_obj,
  group.by = "sample_type",
  ident.1 = "tumor",
  ident.2 = "PDX"
)
head(deg_balanced, 20)

# Optional: EnhancedVolcano plot
# if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
#   BiocManager::install("EnhancedVolcano")
# }
library(EnhancedVolcano)
EnhancedVolcano(
  deg_balanced,
  lab = rownames(deg_balanced),
  x = "avg_log2FC",
  y = "p_val_adj",
  title = "Balanced DE: Tumor vs PDX (Cluster 1)",
  pCutoff = 1e-5,
  FCcutoff = 1.5
)
