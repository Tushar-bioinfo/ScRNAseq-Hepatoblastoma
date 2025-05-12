library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(harmony)
library(dplyr)

basedir <- "/Users/tusharsingh/Downloads/single_RNAseq/workflow/hb"

dirs <- list.dirs(basedir,full.names=F,recursive = F)

seurat_list <- list()


# creating seurat objects 
for (i in dirs){
  
  read <-  Read10X(file.path(basedir,i) )
  
  #GE <- read$`Gene Expression`
  
  seurat_object <- CreateSeuratObject(counts = read , min.cells = 3 , min.features = 200 )
  
  name <- gsub("_filtered_feature_bc_matrix","",i)
  
  seurat_list[[name]] <- seurat_object  
}

rm(read)
### merger

merged_seurat_obj <- merge( x= seurat_list[[1]], 
                            y = seurat_list[-1],
                            add.cell.ids = names(seurat_list),
                            project = "HB")


rm(seurat_object)
rm(seurat_list)
merged_seurat_obj

#### creating sample column 


merged_seurat_obj$sample <- rownames(merged_seurat_obj@meta.data)

merged_seurat_obj$sample_id <- sapply(merged_seurat_obj$sample,function (x) {strsplit(x,"_")[[1]][1]} )

merged_seurat_obj$sample_type <- sapply(merged_seurat_obj$sample,function (x) {strsplit(x,"_")[[1]][2]} )

merged_seurat_obj$sample <- NULL


#unique(merged_seurat_obj@meta.data$sample_id)


merged_seurat_obj$mt_percent <- PercentageFeatureSet(merged_seurat_obj, pattern  = "^MT-")


######### Quality check ####

VlnPlot(filtered_mg_seurat_obj , features = c("nFeature_RNA","nCount_RNA","mt_percent"),ncol=3)

FeatureScatter(filtered_mg_seurat_obj , feature1 ="nCount_RNA",feature2 = "nFeature_RNA" )+ geom_smooth(method = "lm")



########### Filtering ######

filtered_mg_seurat_obj <- subset(merged_seurat_obj, subset = nFeature_RNA >500 & nCount_RNA > 1000 & mt_percent < 10 )

rm(merged_seurat_obj)
filtered_mg_seurat_obj

###### normalize ####
filtered_mg_seurat_obj <- NormalizeData(filtered_mg_seurat_obj)

##### Variable features ######
filtered_mg_seurat_obj <- FindVariableFeatures(filtered_mg_seurat_obj, selection.method = "vst", nfeatures = 2000)



#####plot variables 

#plot <- VariableFeaturePlot(filtered_mg_seurat_obj)
#LabelPoints(plot = plot , points = top, repel = T)


#### Scale ######

#all.genes <- rownames(filtered_mg_seurat_obj)     ## this took too long so we instead ised a smaller genes set
#filtered_mg_seurat_obj <-ScaleData(filtered_mg_seurat_obj, features = all.genes)

filtered_mg_seurat_obj <-ScaleData(filtered_mg_seurat_obj)

#### PCA ####

#filtered_mg_seurat_obj <-RunPCA(filtered_mg_seurat_obj , features = VariableFeatures(filtered_mg_seurat_obj))
filtered_mg_seurat_obj <-RunPCA(filtered_mg_seurat_obj)


####### Elbow plot & heat #########

#DimHeatmap(filtered_mg_seurat_obj ,cells = 1000 , balanced = T , dims =40:50  )

ElbowPlot(filtered_mg_seurat_obj)

##################################
#we did all these steps below like clustering,PCA,etc and Found batch effects. So had to comment these all out 
# and run harmony ( integration ) to remove them
################################


##### cultering ######

#filtered_mg_seurat_obj <- FindNeighbors(filtered_mg_seurat_obj , dims = 1:20)

#filtered_mg_seurat_obj <- FindClusters(filtered_mg_seurat_obj , resolution = c(0.1))


##### checking best resolution 

#DimPlot(filtered_mg_seurat_obj, group.by = "RNA_snn_res.0.1", label = T)


#Idents(filtered_mg_seurat_obj)
#Idents(filtered_mg_seurat_obj) <- "RNA_snn_res.0.1"

#### UMAP ########

filtered_mg_seurat_obj <- RunUMAP(filtered_mg_seurat_obj , dims = 1:20, reduction = 'pca')


#### Plots #########

#DimPlot(filtered_mg_seurat_obj , reduction = "umap", label=T)

before1<- DimPlot(filtered_mg_seurat_obj , reduction = "umap", group.by = "sample_id",label=F ,
                  cols = c("red","green","blue"))

before2 <- DimPlot(filtered_mg_seurat_obj , reduction = "umap", group.by = "sample_type",label=F ,
              cols = c("red","green","blue"))

before1 + before2 

rm(filtered_mg_seurat_obj)

######## Harmony ####
harmony_obj <- filtered_mg_seurat_obj %>%
  RunHarmony(group.by.vars = 'sample_id', plot_convergence = FALSE)

harmony_obj@reductions


harmony_obj

Idents(harmony_obj)
# Harmony embeddings
harmony_obj <- harmony_obj %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = c(0.1))

Idents(harmony_obj)
Idents(harmony_obj) <- "RNA_snn_res.0.1"

# visualize 

## using updated umap values which we got form the harmony 
after <- DimPlot(harmony_obj, reduction = 'umap', group.by = 'sample_id')
after2 <- DimPlot(harmony_obj, reduction = 'umap', group.by = 'sample_type')

after+after2
after2



##### just trying out plotting diff plots
clusters <- DimPlot(harmony_obj, reduction = 'umap',group.by= "seurat_clusters", label = T)
clusters <- DimPlot(harmony_obj, reduction = 'umap', label = T)
clusters


table(Idents(harmony_obj), harmony_obj$sample_type)

########

# Subset to only PDX and tumor
subset_obj <- subset(harmony_obj, subset = sample_type %in% c("PDX", "tumor"))
subset_obj <- JoinLayers(subset_obj)



