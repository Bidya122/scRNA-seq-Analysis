# scRNA-seq-Analysis
Downstream single-cell RNA-seq analysis starting from an RDS Seurat object, performing dimensionality reduction, clustering, and cell-type marker discovery.
# Dataset: GSE183276, 49 samples   [Click to view Dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183276)
Single Cell RNA-Seq performed on kidney tissues from 15 Chronic Kidney Disease (CKD), 12 Acute Kidney Disease (AKI) and 20 healthy reference (Ref) individuals.  
Tools Used: RStudio, Python
*Introduction of Dataset:* ##A subset of**49 scRNA-seq samples** was analyzed from a larger kidney atlas study.
- Data type: Single-cell RNA sequencing    
- Samples: Healthy donors and kidney disease patients    
- Scale: High-dimensional gene expression data across thousands of cells    

This dataset enables:    
- Identification of cell types via clustering    
- Detection of differentially expressed genes (DEGs)    
- Characterization of disease-associated cellular states    
<img width="920" height="891" alt="image" src="https://github.com/user-attachments/assets/7c03d620-e9d5-4af6-ac08-f326efbcb99c" />


I begin by loading the .rds file, which stores the single-cell dataset (typically a Seurat object) including expression matrices and prior preprocessing steps.
Alongside this, I import the metadata file containing cell-wise annotations (e.g., sample origin, condition labels, batch information), which is essential for downstream grouping, visualization, and biological interpretation.
The code to do so on Linux is given here but could be also done on R in the command platform. 
```bash
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183276/suppl/GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Counts_03282022.RDS.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183276/suppl/GSE183276%5FKidney%5FHealthy%2DInjury%5FCell%5FAtlas%5FscCv3%5FMetadata%5F03282022.txt.gz
```
Seurat is a tool that helps us organize and understand data from thousands of individual cells. It takes raw gene activity data and turns it into something we can explore — like grouping similar cells, finding different cell types, and visualizing patterns. The .rds file is used to store a Seurat object, which encapsulates the entire single-cell RNA-seq dataset in a structured format. It typically contains raw or normalized gene expression data, cell metadata, feature annotations, and results from prior analyses such as PCA, clustering, and UMAP. This allows the dataset to be easily reloaded without repeating computationally intensive preprocessing steps.    
# For further steps, I have provided the R code here in chunks. The chunks are self explanatory with the use of comments.
# Installation of Libraries & Packages
```bash
## Check R and Bioconductor versions
getRversion()
BiocManager::version()

## Install core CRAN packages for scRNA-seq analysis
install.packages(c("Seurat","SeuratDisk","dplyr","R.utils","ggplot2","ggExtra","RColorBrewer","openxlsx","scales","HGNChelper","dittoSeq","harmony"))

## Install remotes to allow installation of packages from GitHub
## (required here to install SeuratDisk)
## NOTE (Windows): If this step fails, make sure RTools is installed
## and restart R before re-running the command.
install.packages("remotes")

## Install SeuratDisk from GitHub
remotes::install_github("mojaveazure/seurat-disk")

## Install rlang (needed by ggplot2 and tidyverse)
install.packages("rlang")

## Check whether make is available (from RTools on Windows)
Sys.which("make")

## Install jsonlite from source (requires RTools)
install.packages("jsonlite", type = "source")

## Reinstall ggplot2 to match updated dependencies
install.packages("ggplot2")

## Reinstall dittoSeq after dependency updates
BiocManager::install("dittoSeq")

## Reinstall SeuratDisk from GitHub
remotes::install_github("mojaveazure/seurat-disk")

## Load key packages to confirm installation
library(rlang)
library(ggplot2)
library(dittoSeq)
library(SeuratDisk)

## ---------------------------------------------------------
## Install Bioconductor single-cell infrastructure packages
## This produced permission warnings on Windows
## ---------------------------------------------------------

BiocManager::install(c("SingleCellExperiment","batchelor","zellkonverter"))

## If above doesn't work, Install the same packages into the user library
## (alternative to running R as administrator)
BiocManager::install(c("SingleCellExperiment", "batchelor", "zellkonverter"),
  lib = Sys.getenv("R_LIBS_USER"))

## Install HDF5 support for AnnData / zellkonverter
install.packages("hdf5r")

## Install SummarizedExperiment (core Bioconductor class)
BiocManager::install("SummarizedExperiment")

## Load remaining packages
library(SummarizedExperiment)
library(hdf5r)
```
# Loading the packages
```bash
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(R.utils)  
library(ggplot2)
library(ggExtra)
library(RColorBrewer)
library(openxlsx)
library(dplyr)
library(scales)
library(HGNChelper)
library(dittoSeq)
library(harmony)
library(batchelor)
library(zellkonverter)
library(SingleCellExperiment)
```
The following libraries are loaded to support different aspects of single-cell RNA-seq analysis, including data handling, visualization, batch correction, and interoperability.  

🔬 Core scRNA-seq Analysis  
Seurat – Main framework for single-cell data processing, clustering, and visualization  
SingleCellExperiment – Standard data structure for single-cell data in Bioconductor workflows  

🔄 Data Conversion & Interoperability  
SeuratDisk – Converts between Seurat and other formats (e.g., h5ad)  
zellkonverter – Enables conversion between R and Python (Scanpy/AnnData) formats  

⚖️ Batch Correction & Integration  
harmony – Integrates datasets by correcting batch effects  
batchelor – Alternative methods for batch correction and data integration  

📊 Data Manipulation & Utilities  
dplyr – Data wrangling and filtering  
rlang – Supports tidy evaluation and programming features  
R.utils – General utility functions for file and data handling  

📈 Visualization
ggplot2 – Core plotting library  
ggExtra – Adds marginal plots to ggplot visualizations  
RColorBrewer – Provides color palettes for plots  
scales – Improves axis formatting and scaling  

🧬 Annotation & Reporting  
HGNChelper – Standardizes and corrects gene symbols  
dittoSeq – Simplifies visualization of single-cell data    
openxlsx – Exports results to Excel files  

# Defining the Directories
```bash
inputDir <- "D:/Bidya Work/single/GSE183276/input"
outputDir <- "D:/Bidya Work/single/GSE183276/output"
plotDir <- "D:/Bidya Work/single/GSE183276/plots"
```
In this step, I define directory paths for input data, output files, and generated plots. These directories help organize the workflow by separating raw data, processed results, and visualizations. If you want to create directories, you can use dir.create().    

# Loading the Dataset 
```bash
##readRDS reads the RDS file. Basically we are telling it to read it the way its written and hand it to us in R. This step loads raw scRNA-seq counts, removes obvious low-quality genes and cells, and organizes individual kidney cells into a structured object for downstream quality control and biological analysis.

data_GSE183276 <- readRDS("GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Counts_03282022.RDS") ##reads file and saves it as an object
class(data_GSE183276) ##tells the class of the object
dim(data_GSE183276) # rows(genes)37080 columns(cells)109741

seurat_all <- CreateSeuratObject(counts = data_GSE183276, project = "GSE183276", min.cells = 3, min.features = 200) ##wraps raw counts into a seurat object counts=raw UMI count matrix, project=sets the project name, min.cells=remove genes detected in fewer than 3 cells, min.features=remove cells with fewer than 200 detected genes

colnames(seurat_all@meta.data) ##Lists cell-level metadata columns
seurat_oi <- SplitObject(seurat_all, split.by = "orig.ident") ##Splits one big Seurat object into a list of Seurat objects,Each element corresponds to one value of orig.ident
print(names(seurat_oi)) ##prints the names of the list elements.
data_GSE183276
```
Load raw single-cell RNA-seq count data and convert it into a structured Seurat object using CreateSeuratObject() in which cells were removed with fewer than 200 detected genes and remove the genes present in fewer than 3 cells. This was done to make the data more cleaner and less noisy and saved it as seurat_all with column names as "orig.ident"   "nCount_RNA"   "nFeature_RNA".      
<img width="350" height="105" alt="image" src="https://github.com/user-attachments/assets/a0b00c21-0c4b-4cbf-89b3-48d3c6b288ad" />     

The combined dataset was split into multiple Seurat objects based on orig.ident, which represents the original sample identity assigned during the experiment. This enables sample-wise handling of data, facilitating quality control, normalization, and downstream comparison across different biological conditions. The list of seurat object looks like below:    
<img width="1345" height="170" alt="image" src="https://github.com/user-attachments/assets/65624589-ee86-45ff-804a-6d2bd3c5c7a6" />    

<img width="1486" height="746" alt="image" src="https://github.com/user-attachments/assets/aa570207-0784-43df-b300-165f9e869bb3" />    
The displayed matrix is a sparse representation of gene expression data, where rows correspond to genes and columns to individual cells. The dots (.) indicate zero or negligible expression values, meaning that most genes are not expressed in most cells. This sparsity is a defining feature of single-cell RNA-seq data, reflecting cell-type specificity and transcriptional heterogeneity.  Dots (.) represent zero expression values, highlighting the sparse nature of single-cell gene expression data.     
One cell ≠ expresses all genes,     
Most genes = OFF in a given cell,       
👉 So matrix = mostly zeros → shown as(.)        








