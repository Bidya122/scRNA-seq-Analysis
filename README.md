# scRNA-seq-Analysis
Downstream single-cell RNA-seq analysis starting from an RDS Seurat object, performing dimensionality reduction, clustering, and cell-type marker discovery.
# Dataset: GSE183276, 49 samples   [Click to view Dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183276)
Single Cell RNA-Seq performed on kidney tissues from 15 Chronic Kidney Disease (CKD), 12 Acute Kidney Disease (AKI) and 20 healthy reference (Ref) individuals.  
Tools Used: RStudio, Python    
**Introduction of Dataset:** A subset of 49 scRNA-seq samples was analyzed from a larger kidney atlas study.
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

# Mapping the IDs
```bash
id_map <- c(
  "AKI3010018" = "GSM5554468",
  "AKI3010034" = "GSM5554469",
  "AKI3010123" = "GSM5554470",
  "AKI3010125" = "GSM5554471",
  "AKI3210003" = "GSM5554472",
  "AKI3210034" = "GSM5554473",
  "AKI3210074" = "GSM5554474",
  "AKI3310005" = "GSM5554475",
  "AKI3310006" = "GSM5554476",
  "AKI3410050" = "GSM5554477",
  "AKI3410184" = "GSM5554478",
  "AKI3410187" = "GSM5554479",
  "DKD2710039" = "GSM5554480",
  "DKD2810051" = "GSM5554481",
  "DKD2910006a" = "GSM5554482",
  "DKD2910006b" = "GSM5554483",
  "DKD2910006c" = "GSM5554484",
  "DKD2910010" = "GSM5554485",
  "DKD2910011" = "GSM5554486",
  "DKD2910012" = "GSM5554487",
  "DKD2910013" = "GSM5554488",
  "DKD2910016" = "GSM5554489",
  "DKD3110001" = "GSM5554490",
  "DKD3110035" = "GSM5554491",
  "DKD3110040" = "GSM5554492",
  "DKD3110042" = "GSM5554493",
  "HCKD2910008" = "GSM5554494",
  "HCKD3110000" = "GSM5554495",
  "HCKD3110013" = "GSM5554496",
  "LDPRE0181" = "GSM5554497",
  "LDPRE0194" = "GSM5554498",
  "LDPRE027" = "GSM5554499",
  "LDPRE038" = "GSM5554500",
  "LDPRE0551" = "GSM5554501",
  "LDPRE0621" = "GSM5554502",
  "LDPRE19025" = "GSM5554503",
  "LDPRE1905" = "GSM5554504",
  "LDPRE98sc" = "GSM5554505",
  "LDSample1153EO1" = "GSM5554506",
  "LDSample1153EO2" = "GSM5554507",
  "LDSample1153EO3" = "GSM5554508",
  "LDSample1157EO1" = "GSM5554509",
  "LDSample1157EO2" = "GSM5554510",
  "LDSample1157EO3" = "GSM5554511",
  "LDSample1158EO1" = "GSM5554512",
  "LDSample1158EO2" = "GSM5554513",
  "LDSample1158EO3" = "GSM5554514",
  "LDSample1162EO1" = "GSM5554515",
  "LDSample1162EO2" = "GSM5554516"
)
```
This chunk creates a mapping between the internal sample IDs (e.g., AKI3010018) and their corresponding GEO accession numbers (e.g., GSM5554468). I did this to ensure consistency when working with multiple datasets or metadata sources, allowing to programmatically link experimental data with publicly available information. This mapping can later be used to rename Seurat objects, merge datasets, or subset samples accurately. These GSM IDs could be found in the metadata file. 

# GSM ID Mapping and Per Sample Saving
```bash
# Robust GSM mapping replace internal sample IDs with GEO GSM IDs
mapped_names <- sapply(names(seurat_oi), function(x) {  # For each sample name in the Seurat object list:
  if(x %in% names(id_map)) return(id_map[[x]])  # If sample exists in the ID map, replace with GSM ID
  else return(x)  # keep original if no mapping found
})
# Rename the list elements to the mapped GSM IDs
names(seurat_oi) <- mapped_names 
# This ensures that each Seurat object is now labeled with its GEO accession ID,necessary for merging with downloaded metadata or reporting.

seurat_list <- seurat_oi #Assign the renamed Seurat list to a new variable

length(seurat_list) 
names(seurat_list)

for (s in names(seurat_list)) {  #Loop through each sample in the Seurat list to save it individually
  saveRDS(
    seurat_list[[s]],    # The Seurat object for the current sample
    file = file.path(inputDir, paste0(s, "_seurat.rds"))   # File path and name 
  ) 
  # Saves the object as an .rds file named by the GSM ID
  # This allows loading each sample individually later without processing the full dataset
}
print(names(seurat_list))  
```
In this step, we replace the internal sample IDs with their GEO GSM IDs so that each Seurat object has a consistent, public identifier. After renaming, we save each sample individually as a .seurat.rds file. This way, any changes or updates can be made to a single sample without affecting the rest of the dataset, making the workflow easier to manage.  
<img width="1118" height="444" alt="image" src="https://github.com/user-attachments/assets/5b5e73d9-2739-47f0-9647-ee1020525eea" />   
This proves that we have all the data (49 samples) with their correct GSM IDs.

# Metadata load
```bash
meta <- read.delim("GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Metadata_03282022.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
meta
##header=the first row contains column names, stringasfactors=keeps text as character strings, not factors, check.names= keeps column names exactly as they are in the file
```
<img width="1105" height="397" alt="image" src="https://github.com/user-attachments/assets/e30a6bed-e335-4b4e-8234-3eca33c619be" />    

Then I loaded the sample metadata using read.delim, keeping column names and text exactly as in the file. This metadata provides key information about each sample, like patient, tissue, or condition, for annotating the Seurat objects. This is the metadata file that we downloaded in the start.

# Read all the individual seurat objects
```bash
 ## This chunk looks inside a folder (input) and finds all files that contain single-cell data saved as Seurat objects (.rds files). Each of these files corresponds to one sample. It then loads all these Seurat objects into R and stores them together in a list, so that multiple samples can be handled at once. After loading, each Seurat object is given a meaningful name (usually based on sample) so that it is easy to identify which object belongs to which sample. Finally, it prints the names to confirm that all samples were loaded correctly.
rds_files <- list.files(inputDir, pattern = "_seurat\\.rds$", full.names = TRUE) ##Get full paths to all Seurat RDS files in inputDir and Only files ending with "_seurat.rds" are matched
seurat_list <- lapply(rds_files, readRDS) ##read each RDS file into R.This creates a list of Seurat objects (one per file)
names(seurat_list) <- gsub("_seurat\\.rds$", "", basename(rds_files))  # pattern to remove  # replace with empty string  # extract filename from full path
print(names(seurat_list)) # Print the names of the Seurat objects
```
<img width="1114" height="227" alt="image" src="https://github.com/user-attachments/assets/eef93fae-91c6-482f-b06b-41d557e4f1af" />    

This chunk reads back all the individual .seurat.rds files we saved earlier, loading each sample separately into a list of Seurat objects. Each object is given a meaningful name based on the sample, and printing the names confirms that all samples were loaded correctly.

# combining the seurat objects
```bash
sapply(seurat_list, function(obj) {    # For each Seurat object in the list, count how many duplicated gene names exist
  sum(duplicated(rownames(GetAssayData(obj, layer = "counts"))))
})

# Got 0s in all samples (meaning no duplicate genes are present in any of the samples) and thus no need to run and save clean_seurat_list. Proceed with seurat_list. But because in downstream, we have 'clean_seurat_list' everywhere, we'll do

clean_seurat_list <- seurat_list
print(names(clean_seurat_list))
# clean_seurat_list <- lapply(seurat_list, function(obj) {
#   # Get unique genes from RNA assay
#   counts_mat <- GetAssayData(obj, layer = "counts")
#   unique_genes <- !duplicated(rownames(counts_mat))
#   counts_mat_unique <- counts_mat[unique_genes, ]
#   obj[["RNA"]] <- CreateAssayObject(counts_mat_unique)
#   return(obj)
# })
```
<img width="1116" height="595" alt="image" src="https://github.com/user-attachments/assets/a1511da7-02b5-480e-9385-0ed15153b32d" />    
This step checks each Seurat object for duplicated gene names in the counts matrix. Since no duplicates were found, I used seurat_list as is. To keep downstream code consistent, I assigned it to clean_seurat_list, which will be used in all subsequent analyses.

# Merge all Seurat objects into one because we are not doing it for per sample
```bash
# Update orig.ident to GSM IDs inside each Seurat object
for (s in names(clean_seurat_list)) {
  clean_seurat_list[[s]]$orig.ident <- s  # s is GSM ID
}

seurat_combined <- merge(clean_seurat_list[[1]], y = clean_seurat_list[-1], add.cell.ids = names(clean_seurat_list))
unique(seurat_combined$orig.ident) # Merge all Seurat objects into a single object and prefix cell IDs with sample names

setdiff(id_map, names(clean_seurat_list)) # should return 0 # Check that expected sample IDs match the loaded Seurat objects

saveRDS(seurat_combined, file =  paste0(outputDir,"01_GSE183276_seurat_combined.rds")) # Save the merged Seurat object to disk
#seurat_combined <- readRDS(file = paste0(outputDir, "01_GSE183276_seurat_combined.rds"))

# Explore Combined Seurat Object}
class(seurat_combined)            # 'Seurat'
Assays(seurat_combined)           # (RNA) assays available in the rds
DefaultAssay(seurat_combined)     # default/active (RNA) assay in the rds
dim(seurat_combined)              # 29447 genes/features x 109741 cells
```
<img width="1371" height="328" alt="image" src="https://github.com/user-attachments/assets/93fe8732-7a4f-469c-b690-edc7958f3fc9" />    

Individual Seurat objects (each representing a GSM/sample) were first annotated with their respective sample IDs by updating the orig.ident metadata. This ensures that the origin of each cell is preserved after merging. All objects were then combined into a single Seurat object using merge(), with add.cell.ids used to prefix cell barcodes with sample names to avoid duplication.    
This step brings all cells into a unified dataset, enabling cross-sample comparison while retaining sample-level identity — which is essential for downstream analysis like batch correction, clustering, and differential expression. The final merged object contains 29,447 genes across 109,741 cells, representing the combined cellular landscape of all samples.    

# Add biological condition metadata to merged Seurat object
```bash
# Adding 'Condition' col to the seurat object to do qc plot much better
table(seurat_combined$orig.ident)  # Show the number of cells per sample (orig.ident)
gsm_ids <- unique(seurat_combined$orig.ident)  # Extract unique GSM sample IDs from the merged Seurat object
# Define the biological condition corresponding to each GSM sample
condition_map <- c(
 rep("AKI", 12),
 rep("DKD", 14),
 rep("HCKD", 3),
 rep("Healthy", 20)
)
# Assign GSM IDs as names to the condition mapping vector
names(condition_map) <- gsm_ids

# Create cell-level condition vector (drop names!)
seurat_combined$Condition <- unname(condition_map[seurat_combined$orig.ident])

table(seurat_combined$Condition)    # Display the number of cells per biological condition
sum(is.na(seurat_combined$Condition))    # Check for cells with missing (NA) condition annotations
```
<img width="1366" height="316" alt="image" src="https://github.com/user-attachments/assets/afd43189-4a75-4cf9-afaa-4220a6d26abe" />    
To enable biologically meaningful comparisons, a new metadata column Condition was added to the merged Seurat object. Each sample (orig.ident, i.e., GSM ID) was mapped to its corresponding biological condition (AKI, DKD, HCKD, Healthy). This mapping was then expanded to the cell level, so that every individual cell inherits the condition of the sample it originated from.
Basic checks were performed:    
table(orig.ident) → verifies cell distribution across samples    
table(Condition) → confirms correct assignment across conditions    
sum(is.na(Condition)) → ensures no cells are missing condition labels    

# Calculate QC Metrics for Seurat Object
```bash
head(rownames(seurat_combined), 20)  # View the first 20 row names (typically gene names) in the Seurat object
seurat_combined[["percent.mt"]] <- PercentageFeatureSet(   # Calculate the percentage of mitochondrial genes per cell
 seurat_combined,
 pattern = "^MT-"
)  ## "^MT-" means genes starting with "MT-" The result is stored in the metadata slot as "percent.mt"

seurat_combined[["percent.rb"]] <- PercentageFeatureSet(  # Calculate the percentage of ribosomal genes per cell
 seurat_combined,
 pattern = "^RPL|^RPS"
)  # "^RPL|^RPS" matches genes starting with "RPL" or "RPS". Stored in metadata as "percent.rb"

colnames(seurat_combined@meta.data)  # Show all column names in the metadata to confirm our new columns are added

summary(seurat_combined$percent.mt) 
summary(seurat_combined$percent.rb)
```
<img width="1373" height="164" alt="image" src="https://github.com/user-attachments/assets/786ea519-1e80-4c10-9fed-9e26a56b3ddf" />    

To assess cell quality, additional QC metrics were computed and added to the Seurat object metadata.  
percent.mt: percentage of mitochondrial gene expression per cell  
percent.rb: percentage of ribosomal gene expression per cell  
These were calculated using PercentageFeatureSet() by identifying genes based on naming patterns:  
"^MT-" → mitochondrial genes  
"^RPL|^RPS" → ribosomal protein genes   
Mitochondrial percentage (percent.mt): High mitochondrial gene expression is often a sign of stressed or dying cells, as damaged cells tend to leak cytoplasmic RNA and retain mitochondrial transcripts.      
Ribosomal percentage (percent.rb): Reflects the protein synthesis activity of a cell. Extremely high values may indicate technical bias or overrepresentation of housekeeping functions.    
The relatively wide range of mitochondrial content suggests the presence of both healthy cells (low mt%) and stressed/dying cells (high mt%). Ribosomal content also shows high variability, indicating differences in protein synthesis activity or potential technical bias across cells. These distributions highlight the need for careful threshold selection during QC filtering, rather than applying arbitrary cutoffs.    

# PreQC Visualization










