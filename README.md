# Single_Cell_Tutorial
This repository store the information to set up the working environment

# Introduction
Single-cell genomics is an emerging technique used to study cells at high resolution in the framework of Next Generation Sequencing and Omics. We considered an individual of H. charithonia (Lepidoptera) at a specific developmental stage, the 5th instar. Our analysis will focus on the creation of the forewing disc cell atlas using scRNA-seq data (transcriptome cells' layer).

# Set up the programming environment
## 1. Download R and R studio
R is a free software environment for statistical computing and graphics and runs on a wide variety of UNIX platforms, Windows and MacOS. To downloads R follows the instructions in https://www.r-project.org/.
RStudio is an integrated development environment (IDE) for R and Python. It includes a console, syntax-highlighting editor that supports direct code execution, and tools for plotting, history, debugging, and workspace management. To download R studio follow the instructions in https://posit.co/products/open-source/rstudio/. R studio should be downloaded after the R installation.

## 2. R packages
To configure the computing environment we use Seurat package stored in CRAN repository (https://cran.r-project.org/).

```{r}
install.packages(Seurat) #Install Package Seurat from CRAN
library(Seurat)
```
### 3. Download data and create a working directory
Data generated by 10X cellranger pipeline are stored in google drive following the link:

https://

The dataset is in RData format.

