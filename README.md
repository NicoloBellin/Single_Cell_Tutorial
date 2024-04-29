# Single Cell Tutorial
This repository store the information to set up the working environment in R.

# Introduction
Single-cell genomics is an emerging technique used to study cells at high resolution in the framework of Next Generation Sequencing and Omics. We considered an individual of H. charithonia (Lepidoptera) at a specific developmental stage, the 5th instar. Our analysis will focus on the creation of the forewing disc cell atlas using scRNA-seq data (transcriptome cells' layer).

# Set up the programming environment
## 1. Download R and R studio
*nR is a free software environment for statistical computing and graphics and runs on a wide variety of UNIX platforms, Windows and MacOS. To downloads R follows the instructions in https://www.r-project.org/.
* RStudio is an integrated development environment (IDE) for R and Python. It includes a console, syntax-highlighting editor that supports direct code execution, and tools for plotting, history, debugging, and workspace management. To download R studio follow the instructions in https://posit.co/products/open-source/rstudio/. R studio should be downloaded after the R installation.

## 2. R packages
To configure the computing environment we use Seurat package stored in CRAN repository (https://cran.r-project.org/).

```{r}
install.packages(Seurat) #Install Package Seurat from CRAN
library(Seurat) #Load the package during the session
```
### 3. Download data and create a working directory
Create a working directory (fold) named SCell_Tutorial in your Desktop.
Data generated by 10X cellranger pipeline are stored in google drive following the link:

https://drive.google.com/file/d/1SWxqWCNb-UGmeGz_mPJyO-0naqtWy6F0/view?usp=drive_link

The dataset is in RData format and it must be stored in SCell_Tutorial fold. 


