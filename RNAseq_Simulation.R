library(ggplot2)

################### RNA-seq analysis simulation ################

#Import the dataset
df <- read.table("Data/rnaseq_bulk.txt")


#Check some descriptive statistics
media_per_gruppo <- aggregate(. ~ Status, data = df, FUN = mean)
ds_per_gruppo <- aggregate(. ~ Status, data = df, FUN = sd)

#Make an Histogram
hist(df$PLK1, breaks = 30, 
     col = "skyblue", 
     main = "PLK1 distribution", 
     xlab = "Expression level", ylab = "Frequency")

hist(df$ERBB2, breaks = 30, 
     col = "skyblue", 
     main = "ERBB2 distribution", 
     xlab = "Expression level", ylab = "Frequency")

#Normalize the data
max_normalize <- function(x){
  x / max(x)
}

# Apply Max normalization to each column in the data frame
df_norm<- as.data.frame(lapply(df[,1:15], max_normalize))

#Visualize a gene with a boxplot
boxplot(df_norm$WNT1, xlab="CDK4", ylab="Normalized Expression Level")

#Scale the data
df_scaling <- scale(df_norm,center=TRUE, scale=TRUE) 

#Make the PCA
# Esegui la PCA sui dati
pca <- prcomp(df_scaling, scale. = FALSE)

var_exp <- summary(pca)

#Visualize the proportion of standard deviation explained by each axis
pca_scores <- as.data.frame(pca$x)

#plot
plot(pca_scores$PC1, pca_scores$PC2, 
     col=ifelse(df$Status == "Health", 'blue', 'red'), 
     pch=ifelse(df$Status == "Health", 16, 17), 
     xlab= paste('PC1', round(var_exp$importance[2,1] * 100), '%'),
     ylab= paste('PC2', round(var_exp$importance[2,2] * 100), '%'),
     main='Principal Component RNA-seq')

#Check the most important genes for the first two component
sort(abs(pca$rotation[,1]), decreasing = TRUE)
sort(abs(pca$rotation[,2]), decreasing = TRUE)

#Make the Marker test (Wilcoxon)
# Funzione per eseguire il test di Wilcoxon per una singola variabile
wilcoxon_test <- function(variabile){
  wilcox.test(variabile ~ df$Status)
}

#Esegui il test di Wilcoxon per ciascuna colonna (variabile) del dataframe
risultati <- lapply(df_norm, wilcoxon_test)

# Estrai p-value e statistiche dai risultati
p_values <- sapply(risultati, function(x) x$p.value)
statistiche <- sapply(risultati, function(x) x$statistic)

# Crea un dataframe con i risultati
risultati_df <- data.frame(
  Variabile = names(p_values),
  P_Value = p_values,
  Statistica = statistiche
)

#Visualizza il dataframe dei risultati
print(risultati_df)

#Add the pvalue corrected
risultati_df$P_values_adj <- p.adjust(risultati_df$P_Value,
                                     method = 'bonferroni')

#Check the genes with p values < 0.01
risultati_df[risultati_df$P_values_adj < 0.01, ]

#Make a density with ggplot2
notch1_density <- ggplot(df, aes(x=NOTCH1, fill=Status)) +
                  geom_density(alpha=0.4)+
                  theme_classic()

plk1_density <- ggplot(df, aes(x=PLK1, fill=Status)) +
                geom_density(alpha=0.4)+
                theme_classic()

#Make a violin plot with ggplot2
notch1_violiny <- ggplot(df, aes(y=NOTCH1, x=Status,fill=Status)) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic()

plk1_violiny <- ggplot(df, aes(y=PLK1,  x=Status, fill=Status)) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic()

#Evaluate the co expression
cor.test(df$NOTCH1, df$PLK1, method='spearman')

#Plot the relationship
df_sub <- df[df$Status == 'Illness',]

plot(df_sub$NOTCH1, df_sub$PLK1, type = 'p', pch= 16,col='black',
     xlab='NOTCH1', ylab='PLK1')

#make a linear model 
regressione <- lm(PLK1 ~ NOTCH1, data = df_sub)
summary(regressione)
abline(regressione, col='red')


