library(Seurat)
library(ggplot2)

#Load the Dataset
load("Data/df_5iFW.RData")

#Take the Raw Expression matrix 
df_5iFW@assays$RNA@layers$counts
LayerData(df_5iFW, assay = "RNA", layer = "counts")

#Visualize the dimension
colnames(df_5iFW) #Visualize the cells ID
rownames(df_5iFW) #Visualize the Genes ID

#Another method to visualize the Genes ID
Features(df_5iFW[["RNA"]])

#Get the dimension 
nrow(df_5iFW) #Genes
ncol(df_5iFW) #Cells

#Accessing metadata
metadata_5iFW <- df_5iFW@meta.data

#nFeature = UMI (total mRNA) in a cell
#nCounts = Genes detected

#Visualize the Variables
VlnPlot(df_5iFW, 
        features = c("nFeature_RNA", "nCount_RNA"), 
        ncol = 2)

#Scatterplot between nFeature_RNA and nCount_RNA
FeatureScatter(df_5iFW, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA")

#Using ggplot2
#install.packages("ggplot2")

ggplot(metadata_5iFW, aes(x=nCount_RNA, y=nFeature_RNA))+
  geom_point() +
  theme_classic()

#Using log 10 transformation
ggplot(metadata_5iFW, aes(x= log10(nCount_RNA), y= log10(nFeature_RNA)))+
  geom_point() +
  geom_smooth(method='lm')+
  xlab("UMI") +
  ylab("nGene") +
  theme_classic()


#Compute the Novelty Score
df_5iFW$Novelty_Score <- log10(df_5iFW$nFeature_RNA) / log10(df_5iFW$nCount_RNA)

#Accessing metadata
metadata_5iFW <- df_5iFW@meta.data


################### Quality Control ############
ggplot(metadata_5iFW)+
  geom_point(aes(x= log10(nCount_RNA), 
                 y= log10(nFeature_RNA), 
                 color=Novelty_Score)) +
  xlab("UMI") +
  ylab("nGene") +
  scale_color_continuous(type = "viridis") +
  theme_classic()

#Remove cells :
#-UMI == Use the IQR
#-nGene == Use the IQR
#-Novelty Score 0.80

#### interquartile range (IQR)

############ UMI #############
# Calculate quartiles
q_UMI <- quantile(df_5iFW$nCount_RNA, probs = c(0.05, 0.95))
q_Genes <- quantile(df_5iFW$nFeature_RNA, probs = c(0.05, 0.95))
out_novelty <- ifelse(df_5iFW$Novelty_Score > 0.80 , "Good", "Bad")
  
  
#Plot MAD
ggplot(metadata_5iFW)+
  geom_point(aes(x= log10(nCount_RNA), 
                 y= log10(nFeature_RNA), 
                 color = out_novelty)) +
  geom_vline(xintercept = log10(q_UMI[1]) , color = "red", linetype = "dashed") +
  geom_vline(xintercept = log10(q_UMI[2]), color = "red", linetype = "dashed") +
  geom_hline(yintercept = log10(q_Genes[1]) , color = "blue", linetype = "dashed") +
  geom_hline(yintercept = log10(q_Genes[2]), color = "blue", linetype = "dashed") +
  xlab("UMI") +
  ylab("nGene") +
  theme_classic()

#Remove the Outliers
ncells_start <- ncol(df_5iFW)

df_sub <- subset(df_5iFW, subset = (nCount_RNA < 13105 & nFeature_RNA < 3339))
df_sub <- subset(df_sub, subset = ((nCount_RNA > 530 & nFeature_RNA > 296)))
df_sub <- subset(df_sub, subset = (((Novelty_Score > 0.8))))

#How many cells
n_cells_qc <- ncol(df_sub)

#Save the file
save(df_sub, file="C:/Users/nicol/OneDrive/Desktop/SCell_Tutorial/Data/df_sub.RData")





#Save the file
save(df_sub, file="C:/Users/nicol/OneDrive/Desktop/SCell_Tutorial/Data/df_sub.RData")
