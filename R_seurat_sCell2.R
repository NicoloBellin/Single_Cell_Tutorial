library(Seurat)
library(ggplot2)

#Load the data
load("Data/df_subset.RData")

#Normalization
df_sub <- NormalizeData(df_sub3, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)

#Feature Selection 
df_sub <- FindVariableFeatures(df_sub, 
                               selection.method = "vst", 
                               nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(df_sub), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(df_sub)
plot1 

#Scale Data
df_sub <- ScaleData(df_sub)

#Run the PCA
df_sub <- RunPCA(df_sub, 
                 npcs = 100,
                 features = VariableFeatures(object = df_sub))

#Check the Loadings
VizDimLoadings(df_sub, dims = 1:2, reduction = "pca")

#Visualize the Pca
DimPlot(df_sub, reduction = "pca") + NoLegend()

#Visualize the HeatMap 
#DimHeatmap(df_sub, dims = 1:3, cells = 100, balanced = TRUE)

#Elbow plot
ElbowPlot(df_sub)

#Run the Algorithm of Clustering
df_sub <- FindNeighbors(df_sub, reduction = "pca", dims = 1:50)
df_sub <- FindClusters(df_sub, resolution = 0.8)

#Run UMAP
df_sub <- RunUMAP(df_sub, dims = 1:50)

#Plot 
DimPlot(df_sub, reduction = "umap") +xlab("UMAP 1") + ylab("UMAP 2")

#Find all marker genes
df.markers <- FindAllMarkers(df_sub, only.pos = TRUE)


#Dplyr check and manage the genes
top10 <- df.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() 

#Make an heatmap 
DoHeatmap(df_sub, features = top10$gene) + NoLegend()

#Wnta gene:Hcha.Hcha1001G7804 and cortex gene:Hcha.Hcha1500G10526 
wnta <- FeaturePlot(df_sub, features = "gene:Hcha.Hcha1001G7804") + xlab('UMAP1') + ylab("UMAP2")
cortex <- FeaturePlot(df_sub, features = "gene:Hcha.Hcha1500G10526") + xlab('UMAP1') + ylab("UMAP2")


#Check the co-espressione 
wnta_m <- subset(df_sub, features = 'gene:Hcha.Hcha1001G7804')
cortex_m <- subset(df_sub, features = 'gene:Hcha.Hcha1500G10526')

#Make a data frame
df_coexp <- data.frame(WNTA=wnta_m@assays$RNA@layers$data, CORTEX=cortex_m@assays$RNA@layers$data)

#Subset the cluster
df_coexp <- df_coexp[df_sub@meta.data$RNA_snn_res.0.5 == 0,]


#Plot the relation
plot(df_coexp$WNTA, df_coexp$CORTEX, col='black', pch = 16)
model <- lm(df_coexp$CORTEX ~ df_coexp$WNTA)
summary(model)
abline(model, col= 'red')
