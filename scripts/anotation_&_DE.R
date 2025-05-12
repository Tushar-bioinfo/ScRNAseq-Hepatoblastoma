###############################################################
# this script will perform manual anotation and then Differential 
# expression test


###############################################################

# Setting cluster identities
Idents(harmony_obj) <- "seurat_clusters"
Idents(harmony_obj)

### joinlayers ######## had to do this as findAllMarkers would't run without it

harmony_obj <- JoinLayers(harmony_obj)


# Find markers for each cluster


markers <- FindAllMarkers(harmony_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.1, 
                          logfc.threshold = 0.25)

# Top markers per cluster
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
head(top10)


##### heatmap #######

DoHeatmap(harmony_obj, features = top10$gene) + NoLegend()

table(Idents(harmony_obj), harmony_obj$sample_type)


##### Manually anotating them #####

Idents(harmony_obj)
harmony_obj$celltype <- "unknown"
harmony_obj$celltype[which(Idents(harmony_obj) == "0")] <- "Tumor-like (stem)"
harmony_obj$celltype[which(Idents(harmony_obj) == "1")] <-  "Fetal-like hepatoblastoma"
harmony_obj$celltype[which(Idents(harmony_obj) == "2")] <- "Embryonic tumor-like"
harmony_obj$celltype[which(Idents(harmony_obj) == "3")] <- "Bulk hepatoblastoma"
harmony_obj$celltype[which(Idents(harmony_obj) == "4")] <- "Normal hepatocytes"
harmony_obj$celltype[which(Idents(harmony_obj) == "5")] <- "T cells"
harmony_obj$celltype[which(Idents(harmony_obj) == "6")] <- "Liver endothelial cells"
harmony_obj$celltype[which(Idents(harmony_obj) == "7")] <- "Proliferating tumor cells"
harmony_obj$celltype[which(Idents(harmony_obj) == "8")] <- "Lymphatic cells"
harmony_obj$celltype[which(Idents(harmony_obj) == "9")] <- "Stromal / smooth muscle"




harmony_obj$celltype <- levels(harmony_obj)

##### testing plots #######

clusters <- DimPlot(harmony_obj, reduction = 'umap',group.by = "celltype", label = T)
clusters

FeaturePlot(harmony_obj, features = c("RTN4RL1", "NPSR1", "MLIP-AS1", "RERGL"), reduction = "umap")
VlnPlot(harmony_obj, features = c("RTN4RL1", "NPSR1", "MLIP-AS1", "RERGL"))

clusters <- DimPlot(harmony_obj, reduction = 'umap',group.by = "celltype", label = T)

FeaturePlot(harmony_obj, features = c("LYZ", "CD68", "CD163", "HBB", "CD302"), reduction = "umap")


###### Differntial expression #######



# Set sample type as the identity
Idents(harmony_obj) <- "sample_type"
Idents(harmony_obj) <- "celltype"
Idents(harmony_obj)



# Subsetting only cluster 1 (fetal-like hepatoblastoma)
fetal_like <- subset(harmony_obj, subset = seurat_clusters == 1)

#  counts
table(fetal_like$sample_type)

# Running DEG: tumor vs PDX within cluster 1
deg_fetal <- FindMarkers(fetal_like,
                         ident.1 = "tumor",
                         ident.2 = "PDX",
                         min.pct = 0.1,
                         logfc.threshold = 0.25)  # fast and conservative

head(deg_fetal,20)


FeaturePlot(fetal_like, features = c("CFH", "C5", "GREM2", "AFP"), split.by = "sample_type")

VlnPlot(fetal_like, features = c("CYP3A5", "KCNT2", "CFH", "FAM13A","C5"), group.by = "sample_type")




BiocManager::install("EnhancedVolcano")

library(EnhancedVolcano)
EnhancedVolcano(deg_balanced,
                lab = rownames(deg_balanced),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Balanced DE: Tumor vs PDX (Cluster 1)',
                pCutoff = 1e-5,
                FCcutoff = 1.5)


# Subset to tumor and PDX only
fetal_like_tumor_pdx <- subset(fetal_like, subset = sample_type %in% c("tumor", "PDX"))


VlnPlot(fetal_like_tumor_pdx, features = c("CFH", "C5", "CYP3A5", "AFP"), group.by = "sample_type")

VlnPlot(balanced_obj, features = c("AFP", "CFH", "CYP3A5", "DPYS", "C5"), 
        group.by = "sample_type", pt.size = 0.1)

FeaturePlot(fetal_like_tumor_pdx, features = c("CFH", "CYP3A5", "AFP"), split.by = "sample_type", reduction = "umap")
VlnPlot(fetal_like_tumor_pdx, features = c("GAPDH", "RPLP0","ACTB"), split.by = "sample_type", group.by = "sample_type", pt.size = 0.1)


top_genes <- rownames(head(deg_fetal[order(deg_fetal$p_val_adj), ], 20))
DoHeatmap(fetal_like_tumor_pdx, features = top_genes, group.by = "sample_type") + NoLegend()



# Get indices
pdxcells <- WhichCells(fetal_like, expression = sample_type == "PDX")
tumorcells <- WhichCells(fetal_like, expression = sample_type == "tumor")



####### Testing : donwsampling our Tumor cells counts #####


# Downsample tumor
set.seed(123)
tumor_sampled <- sample(tumorcells, length(pdxcells))

# Combine
balanced_cells <- c(pdxcells, tumor_sampled)
balanced_obj <- subset(fetal_like, cells = balanced_cells)

#  FindMarkers
deg_balanced <- FindMarkers(balanced_obj, group.by = "sample_type", ident.1 = "tumor", ident.2 = "PDX")



head(deg_balanced,20)
