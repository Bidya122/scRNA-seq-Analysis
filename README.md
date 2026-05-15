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

 [Click to view .RMD file](https://github.com/Bidya122/scRNA-seq-Analysis/blob/main/scRNAseq%20Bidya.Rmd)

 
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

All individual Seurat objects corresponding to different samples were automatically loaded into R. Each object represents gene expression data from a single biological sample. The objects were stored in a named list structure to preserve sample identity, which is essential for downstream comparative analysis, integration, and batch effect correction. Sample names were extracted from file names to ensure traceability across the workflow.

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
```bash
# These parameters define quality control rules that help distinguish real, healthy single cells from empty droplets, dying cells, or technical artifacts.
study_id <- "GSE183276"

seurat_combined <- AddMetaData(seurat_combined, metadata = seurat_combined$orig.ident, col.name = "Sample") # Add a new metadata column named "Sample" to the Seurat object
# This duplicates the 'orig.ident' column (original sample identity)

save_violin_plots_separate <- function(      ###Generate violin + boxplots for specified QC metrics
 seurat_obj,    #   - seurat_obj: Seurat object containing data and metadata
 plotDir,       #   - plotDir: directory where plots will be saved
 study_id,      #   - study_id: ID to include in file names
 features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")    #   - features: vector of metadata features to plot 
) {   

 if (!dir.exists(plotDir)) {                
 dir.create(plotDir, recursive = TRUE)
 }

 meta <- seurat_obj@meta.data                   # Extract metadata from Seurat object
 meta$sample <- as.factor(meta$orig.ident)             # Convert sample column to factor for plotting

 for (feat in features) {   # Loop over each feature and generate plots
 if (!feat %in% colnames(meta)) {
 warning(paste("Skipping", feat, "- not found in metadata"))     # Skip the feature if it does not exist in metadata
 next
 }
  # Create a violin plot with overlaid boxplot
 p <- ggplot(meta, aes(x = sample, y = .data[[feat]])) +
 geom_violin(trim = TRUE, fill = "steelblue", alpha = 0.7) +
 geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6) +
 labs(title = feat,x = "Sample",y = feat) + 
 theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
 plot.title = element_text(hjust = 0.5))

 ggsave(filename = file.path(plotDir,paste0(study_id, "_preQC_", feat, "_violin.png")),
 plot = p,width = 16,height = 9,dpi = 400,bg = "white")
 }
}

save_violin_plots_separate(seurat_combined, plotDir, study_id)  # Run the QC plotting function: generates and saves violin + boxplots for key metrics (UMI counts, gene counts, mitochondrial % and ribosomal %) for each sample before filtering

# Scatter plot of nCount vs nFeature with marginal histograms
density_scatter_plot <- function(seurat_obj, filename){     ##again a general function like above; each row = one cell
 df <- data.frame(
 log1p_nCount_RNA = log1p(seurat_obj$nCount_RNA),      #total RNA molecules
 log1p_nFeature_RNA = log1p(seurat_obj$nFeature_RNA),    #no. of genes detected
 Condition = seurat_obj$Condition                        #control vs disease (or groups)
 )
 
 p <- ggplot(df, aes(x = log1p_nCount_RNA, y = log1p_nFeature_RNA, colour = Condition)) +
 geom_point(alpha = 0.3, size = 0.5) +
 theme_minimal() +
 theme(legend.position.inside = c(0.05, 0.95), legend.justification = c("left", "top"), legend.key.size = unit(0.5, 'cm')) +
 guides(color = guide_legend(override.aes = list(size = 5))) + # 'size' here controls the symbol size
 labs(x = "log1p(nCount_RNA)", y = "log1p(nFeature_RNA)")
 
 p <- ggMarginal(p, type = "histogram", fill = "skyblue", bins = 40)
 ggsave(filename,plot = p,width = 8,height = 10,dpi = 400,bg = "white")}

density_scatter_plot(seurat_combined, file.path(plotDir, paste0(study_id, "_preQC_density-scatter.png")))
```
<img width="1336" height="798" alt="image" src="https://github.com/user-attachments/assets/b5012f60-10ac-4e07-bc8c-32bb7c6f3411" />    

<img width="795" height="966" alt="image" src="https://github.com/user-attachments/assets/76562d77-b695-444f-999c-b4deb0ffd1ea" />


For pre-quality-control (pre-QC) assessment, I generated four main visualizations, **Violin + Boxplots** for:  
   - `nCount_RNA` (UMI counts)  
   - `nFeature_RNA` (genes detected per cell)  
   - `percent.mt` (mitochondrial content)   
   - `percent.rb` (ribosomal content)    
The violin plots show a broad distribution of nFeature_RNA and nCount_RNA, with a small subset of cells exhibiting extremely high values, suggesting potential doublets. The majority of cells had moderate gene and count values, indicating good sequencing depth.    
The mitochondrial percentage (percent.mt) remained relatively low across most cells, suggesting minimal cell stress or apoptosis. However, a few cells with elevated mitochondrial content were observed and considered for removal.    
Hence, Pre-QC quality assessment across all samples revealed broadly consistent distributions for nCount_RNA, nFeature_RNA, mitochondrial percentage, and ribosomal gene content. Most cells exhibited moderate gene detection and sequencing depth, indicating good overall capture efficiency. A subset of cells showed elevated RNA counts and gene numbers, suggesting potential doublets, while a small fraction of cells exhibited higher mitochondrial content, indicating possible low-quality or stressed cells. Ribosomal gene proportions remained relatively stable across samples, suggesting minimal technical variation between libraries.

# Quality Control Filtering Functions for Seurat Object
```bash
# Remove cells with abnormal total RNA counts; automatically detects outliers
## Defines functions to remove low-quality or outlier cells, based on RNA counts, feature counts, mitochondrial and ribosomal RNA percentages. 

# Steps: 1. Extract nCount_RNA from metadata
#   2. Compute Mahalanobis distance for each cell
#   3. Identify outlier cells exceeding chi-squared threshold
#   4. Return Seurat object with outliers removed
fil_nCounts <- function(seurat_obj, threshold_mahalanobis){    ##Remove cells with abnormal total RNA counts (nCount_RNA) using Mahalanobis distance 
 ncounts_data <- as.data.frame(seurat_obj@meta.data$nCount_RNA)    # Convert nCount_RNA to a data frame
 colnames(ncounts_data) <- "nCount_RNA"
 
# Compute Mahalanobis distance (distance from multivariate mean considering covariance
##The Mahalanobis distance is a way to measure how far a cell is from a group of cells, taking correlations between genes into account.
 mahalanobis_dist <- mahalanobis(x = ncounts_data, center = colMeans(ncounts_data), cov = cov(ncounts_data))
 
 seurat_obj$mahal_dist_nCount <- mahalanobis_dist    # Add Mahalanobis distances to Seurat metadata
 mahal.fil <- qchisq(threshold_mahalanobis, df = ncol(ncounts_data))   # Compute chi-squared cutoff for outlier detection
 message(paste0("Mahalanobis threshold: ", round(mahal.fil, 3)))
 
 cells_to_remove <- rownames(seurat_obj@meta.data)[seurat_obj$mahal_dist_nCount > mahal.fil]     # Identify cells exceeding threshold
 message(paste0("Removing ", length(cells_to_remove), " cells (outliers by nCount_RNA)"))
 
 return(subset(seurat_obj, cells = setdiff(colnames(seurat_obj), cells_to_remove)))    # Return Seurat object with outlier cells removed
}

# Filter low (empty droplets or poor capture) & high (likely doublets) gene cells
## Remove cells with too few or too many detected genes
#   - seurat_obj: Seurat object
#   - min_feat: minimum number of features required (filters low-quality/empty droplets)
#   - max_feat: maximum number of features allowed (filters likely doublets)
filter_nFeatures <- function(seurat_obj, min_feat, max_feat, ...) {
 return(subset(seurat_obj, subset = nFeature_RNA > min_feat & nFeature_RNA < max_feat))
}

# Filter high mitochondrial RNA- Remove cells with high mitochondrial RNA content
# Rationale: High percent.mt often indicates stressed or dying cells
filter_mt <- function(seurat_obj, max_mt) {
 return(subset(seurat_obj, subset = percent.mt < max_mt))
}

# Filter high ribosomal RNA- Remove cells with high ribosomal RNA content
# Rationale: Extremely high ribosomal RNA may indicate technical artifacts
filter_rb <- function(seurat_obj, max_rb) {
 return(subset(seurat_obj, subset = percent.rb < max_rb))
}
```
A multi-layered quality control strategy was implemented combining multivariate outlier detection and biologically informed thresholds. Mahalanobis distance was used to identify global outliers, while gene count, mitochondrial percentage, and ribosomal content filters were applied to remove low-quality, stressed, or technically biased cells. This approach ensures robust retention of biologically meaningful single-cell profiles while minimizing technical noise.      
To identify globally aberrant cells, a multivariate outlier detection approach using Mahalanobis distance was applied on nCount_RNA. This method accounts for covariance structure in the dataset and identifies cells that deviate significantly from the population distribution. Cells exceeding a chi-square–based threshold were considered potential technical artifacts or doublets and were removed. The nCount_RNA distribution plot in PreQC revealed the presence of cells with extreme RNA content values, suggesting potential doublets or technical artifacts. To robustly capture such deviations, a Mahalanobis distance-based outlier detection was applied instead of a simple univariate threshold. This allows identification of cells that deviate from the global distribution of RNA counts.   
Cells were filtered based on detected gene counts (nFeature_RNA) to remove empty droplets (low gene counts) and potential doublets (abnormally high gene counts), ensuring retention of high-quality single-cell transcriptomes. The nFeature_RNA distribution plot showed a clear central population with low-gene and high-gene outliers. Low gene counts were interpreted as empty droplets or poor-quality cells, while high gene counts likely represented doublets. Based on this distribution, minimum and maximum feature thresholds were defined to retain high-confidence single-cell profiles.    
Cells with elevated mitochondrial RNA content were excluded, as high mitochondrial proportion is indicative of stressed, dying, or low-quality cells resulting from compromised membrane integrity. The mitochondrial gene fraction plot showed that most cells had low mitochondrial content, indicating good cell viability. A subset of cells exhibited elevated mitochondrial percentages, consistent with cellular stress or apoptosis. These cells were removed using a mitochondrial threshold to ensure exclusion of low-quality cells.  
Cells with unusually high ribosomal RNA content were filtered to reduce potential technical artifacts and transcriptional bias associated with abnormal ribosomal enrichment. Ribosomal gene content remained relatively consistent across samples, suggesting minimal technical variation. However, cells with abnormally high ribosomal RNA levels were excluded to reduce potential transcriptional bias and technical artifacts.   

# QC thresholds application + reporting
```bash
##Applies filtering thresholds to remove low-quality or abnormal cells and generates a summary table showing the effect of QC.
study_id <- "GSE183276"
min_features <- 200 # A cell must express at least 200 genes to be considered real. <200 --> empty droplets/cause most noise
max_features <- 7000 # Cells with more than 7000 genes are removed. Extremely high gene counts often indicate doublets
max_percent_mt <- 15 # Remove cells where >20% of RNA comes from mitochondria
max_percent_rb <- 10 # Remove cells dominated by ribosomal expression. High rRNA --> low info content/technical bias
threshold_mahalanobis <- 0.95 # Remove the worst 5% of cells that look abnormal overall (catches subtle outliers)

seurat_filtered <- filter_nFeatures(seurat_combined, min_features, max_features) # Filter cells by number of features

seurat_filtered <- filter_mt(seurat_filtered, max_percent_mt)  # Filter cells by mitochondrial RNA content

seurat_filtered <- fil_nCounts(seurat_filtered, threshold_mahalanobis)  # Filter cells with abnormal total RNA counts (Mahalanobis distance)

seurat_filtered <- filter_rb(seurat_filtered,max_percent_rb)  # Filter cells by ribosomal RNA content

# Cell counts before QC
df_pre_qc <- as.data.frame(table(seurat_combined$Sample))
colnames(df_pre_qc) <- c("Sample_ID", "Cells_Before_QC")

# Cell counts after QC
df_post_qc <- as.data.frame(table(seurat_filtered$Sample))
colnames(df_post_qc) <- c("Sample_ID", "Cells_After_QC")

# Merge into final QC table (left join behavior)
final_qc_table <- merge(df_pre_qc, df_post_qc, by = "Sample_ID", all.x = TRUE)
final_qc_table
total_cells_before_qc <- sum(final_qc_table$Cells_Before_QC, na.rm = TRUE)
total_cells_after_qc <- sum(final_qc_table$Cells_After_QC, na.rm = TRUE)

total_cells_before_qc
total_cells_after_qc
```
<img width="949" height="339" alt="image" src="https://github.com/user-attachments/assets/938747a1-0227-4b67-bc83-0e09f5fcce92" />  

The QC process ensures that downstream analyses (clustering, DEGs, trajectory inference) are performed only on high-quality, biologically meaningful cells, reducing technical noise and improving interpretability. 
For nCounts_RNA, min_features = 200 → removes empty/low-quality droplets and Mahalanobis-based filtering → removes global outliers in RNA counts.   
For the hight nFearture_count Cells with extremely high gene counts which represent doublets/multiplets. This justifies min_features = 200 and max_features = 7000.  
High mitochondrial content indicates stressed or apoptotic cells. Cells above this threshold max_percent_mt = 15% are removed to improve biological signal clarity.  
Excessively high values of rbRNA% suggest low informational content or technical bias hence max_percent_rb = 10%. 
The long right tails in nCount_RNA and nFeature_RNA motivate removal of extreme outliers (doublets). The low-count region on the left supports removal of empty droplets. The spread in percent.mt and percent.rb highlights stressed or low-quality cells. Together, these distributions define the thresholds used in filtering.    
# Save the summary file QC metrics
```bash
qc_metrics_file <- file.path(outputDir, paste0(study_id, "_QC_metrics.csv"))
write.csv(final_qc_table, qc_metrics_file, row.names = FALSE)
cat("Saved QC summary table to:", qc_metrics_file, "\n")
```
<img width="161" height="433" alt="image" src="https://github.com/user-attachments/assets/bd542ab4-1498-4aff-881b-c87a20fd5041" /> 
This QC file makes it clear about the dataset after the QC filters. Quality control filtering significantly reduced the number of cells across samples, retaining approximately 5–10% of high-quality cells. This indicates the presence of substantial low-quality or stressed cells in the raw dataset, which were removed based on feature counts, RNA counts, and mitochondrial content. A few samples showed extremely low retention, suggesting either poor initial quality or overly stringent filtering thresholds.    

# Post QC Visualization
```bash
##Generates violin plots and scatter plots for QC metrics after filtering, allowing visual inspection of cell quality.
study_id <- "GSE183276"

save_violin_plots_separate <- function(seurat_obj,plotDir,study_id,features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")) 

{
 # Ensure output directory exists
 if (!dir.exists(plotDir)) {
 dir.create(plotDir, recursive = TRUE)
 }

 meta <- seurat_obj@meta.data    # Extract metadata for plotting
 meta$sample <- as.factor(meta$orig.ident)      # Create a factor for sample identity (used for x-axis)

 for (feat in features) {         # Loop through features and generate plots
 if (!feat %in% colnames(meta)) {     
 warning(paste("Skipping", feat, "- not found in metadata"))                  # Skip features not found in metadata
 next 
 }
  # Create violin plot with overlaid boxplot
 p <- ggplot(meta, aes(x = sample, y = .data[[feat]])) +
 geom_violin(trim = TRUE, fill = "steelblue", alpha = 0.7) +
 geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6) +
 labs(title = feat, x = "Sample", y = feat) +
 theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
 plot.title = element_text(hjust = 0.5))

 ggsave(filename = file.path(plotDir,paste0(study_id, "_postQC_", feat, "_violin.png")),
 plot = p, width = 16, height = 9, dpi = 400, bg = "white" ) }
}

save_violin_plots_separate(seurat_filtered, plotDir, study_id)

# Scatter plot of nCount vs nFeature with marginal histograms
density_scatter_plot <- function(seurat_obj, filename){
 df <- data.frame(
 log1p_nCount_RNA = log1p(seurat_obj$nCount_RNA),
 log1p_nFeature_RNA = log1p(seurat_obj$nFeature_RNA),
 Condition = seurat_obj$Condition
 )
 
 p <- ggplot(df, aes(x = log1p_nCount_RNA, y = log1p_nFeature_RNA, colour = Condition)) +
 geom_point(alpha = 0.3, size = 0.5) +
 theme_minimal() +
 theme(legend.position.inside = c(0.05, 0.95), legend.justification = c("left", "top"), legend.key.size = unit(0.5, 'cm')) +
 guides(color = guide_legend(override.aes = list(size = 5))) + # 'size' here controls the symbol size
 labs(x = "log1p(nCount_RNA)", y = "log1p(nFeature_RNA)")
 
 p <- ggMarginal(p, type = "histogram", fill = "skyblue", bins = 40)
 ggsave(filename,plot = p,width = 8,height = 10, dpi = 400, bg = "white")
}

density_scatter_plot(seurat_filtered, file.path(plotDir, paste0(study_id, "_postQC_density-scatter.png")))
```
<img width="883" height="519" alt="image" src="https://github.com/user-attachments/assets/1ff82556-3f5f-45cf-9adc-088d8020111a" />    

<img width="797" height="962" alt="image" src="https://github.com/user-attachments/assets/6d16d497-969a-4aac-90ec-4b772b5e5eb2" />    

Post-quality control filtering resulted in a marked improvement in data quality. Cells with low gene counts were effectively removed, leading to a more compact and biologically meaningful distribution of detected features (nFeature_RNA). Similarly, cells with extremely low transcript counts (nCount_RNA) were excluded, eliminating long-tailed distributions associated with low-quality or empty droplets. Additionally, cells exhibiting high mitochondrial gene expression (percent.mt), indicative of cellular stress or apoptosis, were largely removed. Overall, the filtering process enriched for high-quality, transcriptionally active cells suitable for downstream analysis. This indicates that the retained cells represent biologically meaningful transcriptional profiles rather than technical artifacts. 

# Keeping only protein coding gene
```bash
## Restrict Seurat Object to Protein-Coding Genes
## Keeps only protein-coding genes prior to normalization and inspects the filtered Seurat object.
protein_coding_genes <- readLines(file.path(inputDir, "Protein_coding_genes.txt"))

## subsetting a Seurat object by features (genes), not by cells before normalization
seurat_filtered <- subset(seurat_filtered, features = protein_coding_genes)

####  Explore Filtered Seurat Object after QC Filter  ####
class(seurat_filtered)            # Seurat class
Assays(seurat_filtered)           # assays available in the rds
DefaultAssay(seurat_filtered)     # default/active assay in the rds
dim(seurat_filtered)# 16442 genes/features x 5405 cells
```
<img width="528" height="215" alt="image" src="https://github.com/user-attachments/assets/053fadf4-af5d-46b1-89af-e31760f21b46" />

The dataset was restricted to protein-coding genes to reduce transcriptional noise and improve biological interpretability. A curated list of protein-coding genes was used to subset the Seurat object at the feature level, while retaining all cells. This step ensured that downstream analyses focused on functionally relevant genes. I got this list from ensemble.

# Download the GEO metadata.csv file
```bash
library(GEOquery)
gse <- getGEO("GSE183276", GSEMatrix = TRUE)
metadata <- pData(gse[[1]])
colnames(metadata)
setwd("D://Bidya Work/single/GSE183276")
write.csv(metadata, "GSE183276_metadata.csv", row.names = TRUE)
```
# Merge this GEO metadata file into seurat object
```bash
metadata <- read.csv("D://Bidya Work/single/GSE183276/GSE183276_metadata.csv", header = TRUE, row.names = 1)
metadata
# Add GEO metadata to seurat metadata slot}
# Extract current Seurat metadata and keep barcodes
seurat_meta <- seurat_filtered@meta.data
seurat_meta$barcode <- rownames(seurat_meta)  # store barcodes as a column. why? coz merge() will shuffle rows and remove rownames so we temporarily save sample IDs
colnames(seurat_meta)
colnames(metadata)
head(metadata)

# Merge by sample <-> geo_accession
merged_meta <- merge(seurat_meta, metadata, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)   #Join based on sample ID

# Restore barcodes as rownames
rownames(merged_meta) <- merged_meta$barcode   
merged_meta$barcode <- NULL 

# Assign back to Seurat object
seurat_filtered@meta.data <- merged_meta
```
<img width="1453" height="235" alt="image" src="https://github.com/user-attachments/assets/92da2365-18bd-4071-ad36-3a9d7a2d3843" />    

<img width="1462" height="214" alt="image" src="https://github.com/user-attachments/assets/d14253d1-9c90-40f0-ae15-e622ff7c6bfb" />    
So, both in both the files, seurat_meta and the GEO metadata file there is a key difference. The seurat meta file is cell level and the GEO meta file has sample level data. Trick is to map the cells as samples. so in seurat_meta$Sample is same as GSE183276_metadata$geo_accession. So both these need to be linked. So first we stored the barcodes and then merged it, restored the rownames as they are cell barcodes. 

# Save the Seurat Object
```bash
saveRDS(seurat_filtered,file = file.path(outputDir, "02_GSE183276_seurat_qc_filtered.rds"))
```
 # Normalization and Scaling, Standard Workflow
 ```bash
# Normalize gene expression values for each cell by total expression, multiply by a scale factor (default 10,000), and log-transform. This makes expression levels comparable across cells.
cat("Normalizing data using LogNormalize method (scale factor = 10,000)...\n")
seurat_processed <- NormalizeData(seurat_filtered,
                                 normalization.method = "LogNormalize")  # What percentage of total RNA does each gene contribute in this cell?

# Identify highly variable genes across cells using the VST method. These genes capture the most biological signal and are used for downstream analyses.
seurat_processed <- FindVariableFeatures(seurat_processed, selection.method = "vst", nfeatures = 2500)
seurat_processed <- ScaleData(seurat_processed)  # Scale and center the data so that each gene has mean = 0 and variance = 1. This ensures that highly expressed genes do not dominate PCA.

seurat_processed <- RunPCA(seurat_processed, dims=1:100, npcs = 100) # Perform Principal Component Analysis (PCA) for dimensionality reduction. PCA is run on the scaled expression values of variable genes. Here, up to 100 PCs are computed to capture major sources of variation.
```
<img width="965" height="429" alt="image" src="https://github.com/user-attachments/assets/2b24c898-a645-4778-987c-80ac63677906" />    
<img width="845" height="426" alt="image" src="https://github.com/user-attachments/assets/54d94511-e2b8-4124-92ef-691881d9a68c" />
<img width="919" height="399" alt="image" src="https://github.com/user-attachments/assets/2b7075a7-531e-48ec-a0c0-10ca016af197" />    

*Normalization* - Gene expression data were normalized using the LogNormalize method from the Seurat package. During this step, raw UMI counts for each cell were scaled by the total expression, multiplied by a factor of 10,000, and log-transformed. This ensures comparability of gene expression across cells with differing sequencing depths. Normalization was performed across multiple data layers (e.g., counts.1, counts.2, counts.3), with successful completion indicated by progress bars reaching 100% for each layer.

*Identification of Variable Genes* - Highly variable genes (HVGs) were identified using the variance stabilizing transformation (VST) method.   
For each data layer:  
Gene-wise mean expression and variance were computed  
Variance was standardized and clipped to reduce the effect of outliers  
A subset of the most variable genes was selected for downstream analysis  
The progress output (0–100%) confirms successful computation of gene variances and standardized variance across all layers. 

*Data Scaling* - Following identification of highly variable genes, the dataset was scaled using the ScaleData function from Seurat. During this step:  
Gene expression values were centered (mean = 0)  
Expression values were scaled (variance = 1)  
This standardization ensures that genes with high absolute expression levels do not dominate downstream analyses, allowing biologically meaningful variation to drive dimensionality reduction.  

*Dimensionality Reduction using PCA* - Principal Component Analysis (PCA) was performed on the scaled expression matrix of the selected highly variable genes.  
A total of 100 principal components (PCs) were computed, PCA captures the major sources of variation in the dataset. Each PC represents a combination of genes contributing to a specific pattern of variation across cells.   
This step reduces the complexity of high-dimensional gene expression data while preserving biologically relevant signals.    

| Step                 | Parameter    | Meaning                               |
| -------------------- | ------------ | ------------------------------------- |
| NormalizeData        | LogNormalize | Converts counts → relative expression |
| FindVariableFeatures | vst          | Finds biologically variable genes     |
| FindVariableFeatures | 2500         | Number of genes used downstream       |
| ScaleData            | default      | Standardizes gene expression          |
| RunPCA               | npcs = 100   | Computes 100 variation axes           |
| RunPCA               | dims = 1:100 | (Mostly redundant here)               |

# Add custom colors for Conditions in metadata
```bash
custom_colors <- c("AKI"="purple3", "DKD"="darkorange3", "HCKD"="green4", "Healthy"="blue3")
```
# Principal Component Selection and Variance Analysis_Elbow Plot
```bash
# Get the standard deviation for each PC then get the variance explained by each PC, then calculate the cumulative variance
stdev <- seurat_processed[["pca"]]@stdev # Extract the standard deviation of each principal component (PC) from the PCA reduction object.
var_explained <- stdev^2 / sum(stdev^2) # Calculate the proportion of variance explained by each PC. Variance explained is the squared standard deviation of each PC divided by the total variance across all PCs.
cum_var <- cumsum(var_explained)  # Compute cumulative variance explained across PCs.
pca_var_df <- data.frame(PC = 1:length(stdev),  Variance = var_explained,  CumulativeVariance = cum_var) # Store variance metrics in a data frame for inspection and reporting.
pca_var_df
num_PCs <- min(which(cum_var >= 0.95)) # Identify the minimum number of PCs required to explain at least 95% of the total variance in the dataset.
num_PCs

elbow_plot <- ElbowPlot(seurat_processed, ndims = 100, reduction = 'pca')+
  labs(title = "PCA Elbow Plot for Dimensionality Selection")

ggsave(file.path(plotDir, paste0(study_id, "_ElbowPlot.png")),
       elbow_plot, width = 8, height = 6, bg = 'white')
```
<img width="981" height="445" alt="image" src="https://github.com/user-attachments/assets/68bcee07-2217-44bd-80d0-e29c8edd663d" />    
num_PCs = <img width="98" height="72" alt="image" src="https://github.com/user-attachments/assets/3fcb0f59-36ff-4613-bb70-1fa98b4e20ca" />    

<img width="1205" height="905" alt="image" src="https://github.com/user-attachments/assets/b97fe759-c710-44dd-ab17-b495d6ee02eb" />
The Elbow plot indicated a sharp decline in variance explained within the first ~15–20 principal components, followed by a plateau. However, 35 principal components were retained for downstream analysis to ensure sufficient capture of biological variability, as supported by cumulative variance analysis.

# Total variance explained across the top 35 PCs
```bash
pca_stdev <- seurat_processed[["pca"]]@stdev  # Extract the standard deviation for each principal component (PC) from the PCA reduction object.
pca_var_explained <- (pca_stdev^2) / sum(pca_stdev^2) * 100  # Calculate the percentage of variance explained by each PC. Variance explained is computed as the squared standard deviation and divided by the total variance, multiplied by 100.
total_var<- sum(pca_var_explained[1:35])  # Compute the total variance explained by the top 35 PCs. This helps quantify how much biological variation is retained when using 35 dimensions for downstream analyses.

cat("Total variance explained by top PCs:", round(total_var, 2), "%\n")   # Print the total variance explained by the top PCs.
##60–85% variance explained is very typical
##100% is neither possible nor desirable (that would mean you kept all the noise)
##80% is a healthy balance between signal retention and noise reduction
```
<img width="530" height="69" alt="image" src="https://github.com/user-attachments/assets/22230aee-d2b9-4848-8071-c62089d42181" />

The cumulative variance explained by the top 35 principal components (PCs) was calculated using Seurat to assess the proportion of total variability captured after dimensionality reduction. This analysis confirmed that the selected PCs capture a substantial proportion of the total variance. In single-cell RNA-seq analysis, retaining approximately 60–85% of total variance is considered optimal, as it balances biological signal retention with noise reduction. Therefore, 35 principal components were retained for downstream analysis, ensuring that sufficient biological variation is preserved for clustering and visualization.    

# Perform UMAP & Clustering on unintegrated data and UMAP plot save. 
```bash
library(glue)
seurat_processed <- RunUMAP(seurat_processed, dims = 1:35)  # UMAP provides a low-dimensional (2D) representation of the data that preserves local neighborhood structure for visualization.
seurat_processed <- FindNeighbors(seurat_processed, dims = 1:35, reduction = "pca", graph.name = "pca_nn") # Construct a shared nearest-neighbor (SNN) graph using PCA embeddings. This graph represents cell–cell similarity and is used for clustering.

seurat_processed <- FindClusters(seurat_processed, resolution = 0.8,graph.name = "pca_nn", cluster.name = "pca_clusters")  # Perform graph-based clustering on the PCA-derived neighbor graph. The resolution parameter controls cluster granularity (higher values yield more clusters). 
n_clusters <- length(unique(seurat_processed$pca_clusters))  # Count and report the number of identified clusters.
cat(glue("Clustering complete. Number of clusters: {n_clusters}\n"))  # 19 clusters found
table(seurat_processed$pca_clusters)
```
<img width="896" height="637" alt="image" src="https://github.com/user-attachments/assets/0f59acec-7081-439b-9aaf-d5b60ecd3923" />    

table(seurat_processed$pca_clusters)
<img width="633" height="39" alt="image" src="https://github.com/user-attachments/assets/d7a2905b-96b7-4a6b-a15a-37250e5ded88" />    

```bash
# UMAP colored by PCA-based clusters, shows the cell types, so cells that are similar stay close together in 2D. Each cluster (represented by a color) groups cells with similar overall gene expression patterns.

umap_clusters <- DimPlot(
  seurat_processed,
  reduction = "umap",
  group.by = "pca_clusters",
  label = TRUE
) +
  labs(title = "UMAP: PCA-Based Clusters")

ggsave(filename = file.path(plotDir, paste0(study_id, "_UMAP_pca_clusters.png")), plot = umap_clusters,  width = 8,  height = 6,  bg = "white")
```
<img width="1205" height="908" alt="image" src="https://github.com/user-attachments/assets/b0464fb1-16b1-44a2-b3f9-dc3f75613933" />    

UMAP visualization was performed using the top 35 principal components to project high-dimensional gene expression data into two-dimensional space using Seurat. The resulting plot revealed 19 distinct clusters, corresponding to transcriptionally distinct cell populations. Clusters were observed to be well-separated and compact, indicating effective grouping of cells with similar gene expression profiles. Minimal overlap between clusters further supports the robustness of the clustering approach. These observations validate the selection of 35 principal components and a clustering resolution of 0.8, confirming that the chosen parameters effectively capture biologically meaningful variation while minimizing noise. 

```bash
head(seurat_processed@meta.data, 3)
```
<img width="1506" height="177" alt="image" src="https://github.com/user-attachments/assets/2755691d-d081-4670-9e56-71dd5f236bc4" />    

# Clustering by condition
```bash
# UMAP colored by condition (if available in metadata) 
umap_condition <- DimPlot( seurat_processed, reduction = "umap", group.by = "Condition" ) + labs(title = "UMAP: Condition") 
ggsave( filename = file.path(plotDir, paste0(study_id, "_UMAP_condition.png")), plot = umap_condition, width = 8, height = 6, bg = "white" )
```
<img width="1206" height="905" alt="image" src="https://github.com/user-attachments/assets/0c307d50-e681-49ef-9235-5b43ddea7de9" />    

To assess the influence of experimental conditions on clustering, UMAP visualization was colored by sample condition (AKI, DKD, HCKD, and Healthy) using Seurat. Cells from different conditions were observed to be well distributed across clusters, with no clear condition-specific segregation. Each cluster contained a mixture of cells from multiple conditions. This indicates that the clustering is primarily driven by underlying transcriptional identity (cell types) rather than condition-specific effects. These results suggest minimal batch or condition-driven bias, confirming the robustness of the clustering approach. The absence of condition-specific clustering suggests that downstream differential expression analysis can be performed within clusters to identify condition-associated transcriptional changes.    

# Save the unintegrated Seurat Object
```bash
saveRDS(seurat_processed, file = paste0(outputDir,"03_GSE183276_seurat_pca_umap.rds"))
#seurat_processed <- readRDS(file = paste0(outputDir, "03_GSE183276_seurat_pca_umap.rds"))
#seurat_processed contains Gene expression matrix (genes × cells), Cell clusters (tumor subtypes, immune cells, etc.), Metadata (sample info, condition), Dimensional reductions (PCA, UMAP)

merged = JoinLayers(seurat_processed) # Join layers in the Seurat object (e.g., counts, normalized data, scaled data). Combines multiple assay layers into one unified structure
sce <- as.SingleCellExperiment(merged, assay = "RNA")  # Convert the Seurat object to a SingleCellExperiment (SCE) object
# This allows usage of Bioconductor tools and compatibility with other R workflows

dim(sce) #16442  5405
assayNames(sce)
reducedDimNames(sce)
head(colData(sce),2)

outputDir <-  "D:/Bidya Work/single/GSE183276/output"
h5seurat_name <- "03_GSE183276_seurat_pca_umap.h5seurat"
h5ad_name <- "03_GSE183276_seurat_pca_umap.h5ad"

SaveH5Seurat(object = seurat_processed, filename = file.path(outputDir, h5seurat_name), overwrite = TRUE, version = "3")  # Save Seurat object in H5Seurat format (efficient storage for large datasets)

writeH5AD(sce, file = paste0(outputDir, "03_GSE183276_seurat_pca_umap.h5ad")) # Save SCE object in H5AD format (AnnData, used in Python / Scanpy)
```
<img width="1442" height="392" alt="image" src="https://github.com/user-attachments/assets/2df6d45f-cd71-4f4a-bd9a-d60980b726d9" />

The processed Seurat object was converted to a SingleCellExperiment after joining assay layers to ensure all gene expression data (raw and normalized) were aligned correctly. This preserves biologically relevant information (~16K genes across ~5.4K cells), including cell metadata and PCA/UMAP embeddings for studying tumor heterogeneity in OSCC. The data was then exported to H5Seurat and H5AD formats to enable flexible downstream analysis across R (Seurat/Bioconductor) and Python (Scanpy).    

# Run Harmony & Perform Clustering on Harmony integrated data
```bash
# Run Harmony to integrate data across batches (here, "Sample" is the batch variable)
# Harmony adjusts PCA embeddings to remove batch effects while preserving biological variation
harmony_processed <- RunHarmony(seurat_processed, c("Sample"), plot_convergence = TRUE)

# Compute UMAP based on Harmony-corrected embeddings (low-dimensional visualization)
# Using the first 50 Harmony dimensions
harmony_processed <- RunUMAP(harmony_processed, reduction = "harmony", dims = 1:50)

# Construct a nearest-neighbor graph from Harmony embeddings for clustering
harmony_processed <- FindNeighbors(harmony_processed, reduction = "harmony", dims = 1:50, graph.name = "harmony_nn")

# Perform graph-based clustering on the Harmony-corrected neighbor graph
# The resolution parameter controls the number of clusters (higher = more clusters)
harmony_processed <- FindClusters(harmony_processed, graph.name = "harmony_nn", resolution = 0.8, group.name = "Harmony_clusters")
length(unique(harmony_processed$seurat_clusters))

#found 20 clusters

saveRDS(harmony_processed, file = paste0(outputDir,"04_GSE183276_harmony_corrected.rds"))
harmony_processed <- readRDS(file = paste0(outputDir, "04_GSE183276_harmony_corrected.rds"))

merged1 = JoinLayers(harmony_processed) 
sce_harmony <- as.SingleCellExperiment(merged1, assay = "RNA")

dim(sce_harmony)
assayNames(sce_harmony)
reducedDimNames(sce_harmony)
head(colData(sce_harmony),2)

outputDir <- "D:/Bidya Work/single/GSE183276/output"
h5seurat_name1 <- "04_GSE183276_harmony_corrected.h5seurat"
h5ad_name1 <- "04_GSE183276_harmony_corrected.h5ad"

SaveH5Seurat(
  object = harmony_processed,
  filename = file.path(outputDir, h5seurat_name1),
  overwrite = TRUE,
  version = "3"
)

writeH5AD(sce_harmony, file = paste0(outputDir,  "04_GSE183276_harmony_corrected.h5ad"))
```
<img width="379" height="471" alt="image" src="https://github.com/user-attachments/assets/0b3b406e-7c9d-4148-9646-83908481880e" />
<img width="1322" height="775" alt="image" src="https://github.com/user-attachments/assets/99b5aa5c-91f5-42ee-a377-2469cdfc2fdc" />
Harmony convergence plot shows a steady decrease in the objective function across iterations, indicating effective removal of batch effects. The plateau at later steps confirms convergence, suggesting that further correction does not significantly improve integration. Each point represents the value of the Harmony objective function at a given iteration step during batch correction. The progressive decrease in these values indicates reduction of batch effects across samples. The plateau at later steps shows convergence, meaning further iterations do not significantly improve integration.    

Although clustering was not strongly driven by experimental condition, technical batch effects arising from sample-specific variation may still influence the data. To address this, Harmony integration was applied using “Sample” as the batch variable.

Harmony was applied to the processed Seurat object to correct for batch effects across samples (using the “Sample” metadata), ensuring that downstream analysis reflects biological variation rather than technical differences. The corrected embeddings were then used to compute a UMAP projection (using the first 50 dimensions) for low-dimensional visualization of cell populations. A k-nearest neighbor graph was constructed in the Harmony space to capture cell–cell similarity, followed by graph-based clustering (resolution = 0.8), which identified 20 distinct transcriptional clusters. The processed object was saved for reproducibility and later reloaded, after which it was converted into a SingleCellExperiment format to enable compatibility with Bioconductor tools. Basic structure and metadata were inspected to verify integrity, and finally, the dataset was exported in both H5Seurat and H5AD formats to support cross-platform analysis in R and Python environments.

# Understanding Harmony
```bash
p1 <- DimPlot(
  seurat_processed,
  group.by = "Sample",
  shuffle = TRUE,
  pt.size = 0.5
) + NoLegend()


p2 <- DimPlot(
  harmony_processed,
  group.by = "Sample",
  shuffle = TRUE,
  pt.size = 0.5
) + NoLegend()

ggsave(file.path(plotDir, "before_harmony.png"), plot = p1, width = 8, height = 7, dpi = 300)
ggsave(file.path(plotDir, "after_harmony.png"), plot = p2, width = 8, height = 7, dpi = 300)

p3 <- DimPlot(harmony_processed, group.by = "seurat_clusters", label = TRUE)
ggsave(file.path(plotDir, "clusters.png"), plot = p3, width = 6, height = 5, dpi = 300)
```
<img width="763" height="820" alt="image" src="https://github.com/user-attachments/assets/b9d4028a-d7f6-4f7e-b8a3-729db8151fa3" />   

To correct for potential sample-specific batch effects, Harmony integration was applied using “Sample” as the batch variable.  
UMAP visualization before integration showed partial segregation of cells by sample, indicating technical variation influencing clustering.  
After Harmony correction, cells from different samples were well mixed across clusters, suggesting successful removal of batch effects while preserving biological structure.    
Subsequent graph-based clustering identified ~20 distinct cell populations, representing transcriptionally defined cell types or states within the dataset.  
These clusters form the basis for downstream biological interpretation, including cell type annotation and condition-specific differential expression analysis.    

# Dimention reduction plot on PCA
```bash
##Visualize uncorrected cell clusters on UMAP colored by Sample and Condition, and save the plots
p <- DimPlot(seurat_processed, reduction = "umap", group.by = "Sample") + ggtitle("Uncorrected") +  NoLegend()
  ggtitle("Uncorrected")
ggsave(filename = file.path(plotDir, "BENCHMARKING_GSE183276_RawPCA_sample.png"), width = 7, height = 6, dpi = 600)
p

p <-DimPlot(seurat_processed, reduction = "umap", group.by = "Condition") +
  ggtitle("Uncorrected")
ggsave(filename = file.path(plotDir, "BENCHMARKING_GSE183276_RawPCA_condition.png"), width = 8, height = 6, dpi = 600)
p
```
<img width="428" height="622" alt="image" src="https://github.com/user-attachments/assets/72e271e2-cc2c-42f8-866e-b78b5abae970" />   

I visualized the uncorrected UMAP embeddings colored by both sample identity and biological condition (Tumor vs Normal) to assess the presence of batch effects and underlying biological structure. The sample-wise distribution showed substantial overlap across samples, suggesting minimal batch-driven separation. Additionally, condition-based visualization indicated that biological differences were already captured in the uncorrected space. Despite the apparent low batch effect, Harmony integration was performed to remove any subtle technical variation and ensure robust downstream clustering and biological interpretation.    

# visualizing post-integration (Harmony-corrected) UMAP
```bash
#harmony_processed <- RunUMAP(harmony_processed, reduction = "harmony", dims = 1:30)
p2 <-DimPlot(harmony_processed, reduction = "umap", group.by = "Condition") + ggtitle("Harmony") 
ggsave(filename = file.path(plotDir, "BENCHMARKING_GSE183276_Harmony_condition.png"), width = 8, height = 6, dpi = 600)
p2

p <- DimPlot(harmony_processed, reduction = "umap", group.by = "Sample") + ggtitle("Harmony") +  NoLegend()
ggsave(filename = file.path(plotDir, "BENCHMARKING_GSE183276_Harmony_sample.png"), width = 7, height = 6, dpi = 600)
p

umap_hcoords <- Embeddings(harmony_processed, "umap")
write.csv(umap_hcoords, file="GSE183276_umap_coordinates.csv")
```
<img width="430" height="606" alt="image" src="https://github.com/user-attachments/assets/be4c451f-b38e-4543-9e8f-383fc0eb269a" />

Harmony-corrected UMAP embeddings were visualized to assess the effectiveness of batch correction and preservation of biological structure. Condition-based visualization (Tumor vs Normal) demonstrated that biological differences were retained after integration. Sample-wise visualization showed improved mixing of cells across samples, indicating successful removal of batch effects. Notably, the overall UMAP structure remained similar before and after integration, suggesting that the dataset exhibited minimal batch-driven variation and that Harmony correction did not distort the underlying biological signal. UMAP coordinates were also exported to enable reproducibility and downstream analyses.    

<img width="569" height="601" alt="image" src="https://github.com/user-attachments/assets/36ca1c4b-4917-4c9a-9dad-6e9b2f5ff479" />    

UMAP coordinates derived from Harmony-corrected embeddings were exported for each cell. These coordinates represent a two-dimensional projection (UMAP_1 and UMAP_2) of the high-dimensional gene expression space, enabling visualization of transcriptional similarity between cells. Cells with similar expression profiles are positioned closer together, while distinct populations are separated in the embedding. The two axes do not correspond to specific genes but capture the major sources of variation in the dataset in a reduced-dimensional space. Each row corresponds to a single cell identified by its unique barcode, facilitating reproducible downstream analysis and integration with external tools such as Scanpy.    

# KNN-Based Batch Mixing Function Definition
```bash
##what we did above with umap harmony and all, was just visualising it but here in this chunk we are measuring it by KNN. UMAP provides qualitative visualization, but to objectively evaluate batch correction, we computed KNN-based batch mixing scores before and after Harmony.

# Compute the fraction of nearest neighbors from the same batch before and after Harmony integration, then combine results for visualization as a bar plot
library(ggplot2)   
library(FNN)    # for nearest neighbor calculations
# Computes the mean fraction of k-nearest neighbors (KNN) that belong to the same batch
# Arguments:
#   seu       : Seurat object
#   batch_var : Column in metadata specifying batch/sample
#   reduction : Dimensionality reduction to use (PCA or Harmony)
#   dims      : Dimensions to use from reduction
#   k         : Number of neighbors for KNN
compute_knn_batch_mixing <- function(
    seu,
    batch_var,
    reduction = "pca",
    dims = 1:50,
    k = 20
){
  
# Check if the batch column exists in metadata(We need this to know which batch each cell belongs to)
  if (!batch_var %in% colnames(seu@meta.data)) {
    stop(paste0("Metadata column '", batch_var, "' not found in seu@meta.data"))
  }
  
   # # Check if the required dimensional reduction (PCA/Harmony) exists. If not, create it from scratch
  if (!reduction %in% Reductions(seu)) {
    message("[INFO] PCA not found. Running Normalize → HVG → Scale → PCA")
    seu <- NormalizeData(seu, verbose = FALSE)     # Normalize gene expression (make cells comparable)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)    # Select important genes (highly variable genes)
    seu <- ScaleData(seu, verbose = FALSE)     # Scale data (center and standardize expression values)
    seu <- RunPCA(seu, npcs = max(dims), verbose = FALSE)   # Run PCA to reduce dimensionality (convert genes → coordinates)
  }
  
# Extract coordinates of each cell from PCA or Harmony space(Each cell becomes a point in multi-dimensional space)
  emb <- Embeddings(seu, reduction = reduction)[, dims, drop = FALSE]  
# Find k nearest neighbors for each cell(who are the closest cells in expression space)
   nn  <- FNN::get.knn(emb, k = k)$nn.index
# Get batch label for each cell(which sample/batch each cell belongs to)  
   labs <- seu@meta.data[[batch_var]]
  
   
# For each cell:Check how many of its neighbors belong to the same batch
  same <- sapply(seq_len(nrow(nn)), function(i){
# Compare neighbor batch labels with the cell’s own batch. Calculate fraction of neighbors from same batch
     mean(labs[nn[i,]] == labs[i], na.rm = TRUE)
  })
# Create a dataframe with: batch label + same-batch fraction for each cell  
  df <- data.frame(batch = labs, frac_same_batch = same)
  
  # Compute average same-batch fraction for each batch(final summary per batch)
  aggregate(frac_same_batch ~ batch, df, mean)
}

# Simple bar plot showing mean same-batch fraction per batch
plot_knn_batch_mixing <- function(mix_df){
  ggplot(mix_df, aes(x = batch, y = frac_same_batch)) +
    geom_col(fill = "steelblue") +
    labs(
      title = "Batch Mixing",
      x = "Batch",
      y = "Mean same-batch fraction (higher = stronger batch effect)"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p_batch <- plot_knn_batch_mixing(mix_df)

ggsave(
  filename = file.path(plotDir, "BatchMixing.png"),
  plot = p_batch,
  width = 14,
  height = 5,
  dpi = 600
)
```
After preprocessing and generating the initial PCA-based representation of the data, I observed potential batch-driven separation across samples. To address this, we applied Harmony integration using the Sample metadata as the batch variable, allowing us to correct for technical variation while preserving underlying biological signals. Following integration, I recomputed the neighborhood graph and performed clustering using the Harmony-corrected embeddings to ensure that downstream structure reflected biology rather than batch. We then visualized the data using UMAP both before (PCA-based) and after (Harmony-based) correction, coloring cells by sample and condition. Which, I have shown in the previous step.     
While the UMAP suggested improved mixing across batches, visual inspection alone can be subjective and potentially misleading. Therefore, to make a more objective decision about the effectiveness of batch correction, we implemented a KNN-based batch mixing metric. This approach quantifies, for each cell, the fraction of its nearest neighbors that belong to the same batch, providing a direct numerical measure of batch effect. We computed this metric in PCA space (representing the uncorrected structure) and in Harmony space (representing the corrected structure), keeping parameters consistent to ensure a fair comparison. The results were then combined and visualized as a comparative bar plot, enabling us to assess whether Harmony meaningfully reduced batch-driven clustering across samples.

<img width="1718" height="638" alt="image" src="https://github.com/user-attachments/assets/f4289e84-8059-4e04-9796-62c9e04b3996" />    

Bar plot showing the mean fraction of same-batch nearest neighbors for each sample in PCA space (before batch correction). Higher values indicate stronger batch effects, meaning cells preferentially cluster with cells from the same sample rather than mixing across samples. Several samples exhibit high same-batch fractions, suggesting substantial batch-driven structure in the data prior to integration. To get a comparative analysis, I plotted the before and after results of Harmony integration with computing the fraction of its nearest neighbors that belong to the same batch, providing a direct numerical measure of batch effect.    

# Before and after Harmony + final plot
```bash
# -------------------------------
# Compute batch mixing before and after Harmony
# -------------------------------

# Define which metadata column represents batch/sample(this tells the function how cells are grouped)
batch_column <- "Sample"

# Before batch correction
mix_df <- compute_knn_batch_mixing(
  seu       = seurat_processed,             
  batch_var = batch_column,
  reduction = "pca",   ##PCA is built from:raw gene expression (after normalization). no batch correction
  dims      = 1:50,
  k         = 20
)

# After Harmony batch correction
mix_df_harmony <- compute_knn_batch_mixing(
  seu       = harmony_processed,
  batch_var = batch_column,
  reduction = "harmony",   ##Harmony takes PCA and adjusts it. Cells are repositioned after removing batch effect
  dims      = 1:50,
  k         = 20
)
#so pca is the original structure and harmony the corrected structure.

# Stack the two results together(BUT this alone does not tell us which is before/after)
mix_df_combined <- rbind(mix_df, mix_df_harmony)

# Combine results for plotting
mix_df_combined1 <- rbind(
  transform(mix_df, Status = "Before"),
  transform(mix_df_harmony, Status = "Harmony")
)

##Visualize the effect of Harmony batch correction by plotting mean fraction of same-batch neighbors

batch_cor_plot <- ggplot(mix_df_combined1, aes(x = batch, y = frac_same_batch, fill = Status)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    title = "Batch Mixing",
    x = "Samples",
    y = "Mean same-batch fraction (higher = stronger batch effect)"
  ) +
  scale_fill_manual(
    values = c(
      "Before" = "#00CED1",     # turquoise
      "Harmony" = "#9370DB"     # new purple shade (example 3rd color)
    )
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

batch_cor_plot

ggsave(filename = file.path(plotDir,"batch_correction_barplot_combined.png"), plot = batch_cor_plot, width = 23, height = 8, dpi = 300)
```
After defining the KNN-based batch mixing metric, I applied it to quantitatively evaluate the effectiveness of Harmony integration. Specifically, then computed the fraction of same-batch nearest neighbors for each cell in both the uncorrected (PCA) and corrected (Harmony) embedding spaces, ensuring that identical parameters (dimensions and number of neighbors) were used for a fair comparison. This allowed us to directly assess how cell neighborhood composition changes after batch correction. The results were then aggregated at the sample level and combined into a single dataframe with labels indicating “Before” and “Harmony”, enabling side-by-side comparison. By visualizing these values in a grouped bar plot, It could be clearly observe, the reduction in same-batch neighbor fractions across samples. This comparison was crucial for decision-making, as it provided objective evidence that Harmony successfully reduced batch-driven clustering observed earlier in PCA space. At the same time, any samples that retained relatively higher values indicated residual batch effects, guiding further evaluation if needed, but as per the graph there were none. 

<img width="1419" height="519" alt="image" src="https://github.com/user-attachments/assets/f87d950c-8f17-4d04-a674-cc7650c77e29" />    

Bar plot comparing the mean fraction of same-batch nearest neighbors for each sample before (PCA) and after (Harmony) batch correction. A decrease in same-batch fraction following Harmony indicates improved batch mixing and successful integration, while persistently higher values suggest residual batch effects.

# Cell type Annotation (Cell Typist)
After batch correction and dimensionality reduction, automated cell type annotation was performed using CellTypist (vX.X).  
- Input: Batch-corrected AnnData object (post-Harmony / integration)  
- Model used: Immune_All_Low.pkl (or whichever you used)  
- Mode: majority voting (or probabilistic)  
- Output: Predicted cell labels stored in `.obs['cell_type']`  
This step assigns biologically meaningful identities to clusters based on reference transcriptomic profiles.
So, the following protocol starts on Jupyter Notebook and Anaconda prompt.

**You can download Anaconda from here**   [Click to view website](https://www.anaconda.com/download)

**Open Anaconda Prompt and type**
```bash
# Create a new conda environment named "celltypist_env" with Python 3.9
conda create -n celltypist_env python=3.9
```
```bash
# Activate the newly created environment
conda activate celltypist_env
```
```bash
# Install required packages for scRNA-seq analysis and notebook support
pip install celltypist scanpy jupyter ipykernel
```
```bash
# Register the environment as a Jupyter kernel for use in notebooks
python -m ipykernel install --user --name celltypist_env --display-name "Python (celltypist_env)"
```
```bash
# open Jupyter Notebook
jupyter notebook
```
<img width="602" height="447" alt="image" src="https://github.com/user-attachments/assets/a9af4fd6-3da1-42d5-89ea-264ff85a296a" />    

The notebook interface will open in a web browser, where you can create and run notebooks using the configured environment. Select the environment you just created as the kernel. 

 [Click to view python script file](https://github.com/Bidya122/scRNA-seq-Analysis/blob/main/02_scRNAseq_Bidya.ipynb)
```bash
import celltypist  ## Import celltypist to perform automated cell type annotation using pre-trained models
import scanpy as sc   ## Import scanpy, the primary toolkit for analyzing single-cell RNA-seq data in Python

print("CellTypist and Scanpy loaded successfully!")
```
```bash
# Import pandas to handle data manipulation and analysis, 
# especially for managing cell metadata (obs) and gene information (var) as DataFrames
import pandas as pd
```
```bash
import warnings   ## Import the warnings module to manage and silence non-critical system alerts

# Suppress PerformanceWarnings from pandas, which often occur when highly 
# fragmented DataFrames are created during large-scale single-cell data processing
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
```
```bash
 # Import the os module to interact with the operating system, 
# allowing for file path management and directory navigation
import os
```
```bash
##ad is reserved for AnnData (the library that handles the data structure Scanpy uses).
import scanpy as ad ## Import scanpy using the 'ad' alias
```
```bash
# Set the current working directory to the specified folder for saving plots
os.chdir("D:/Bidya Work/single/GSE183276/plots/celltypist_plots")
os.getcwd()
```
This section initializes the analysis environment by importing the essential Python libraries used in single-cell RNA sequencing (scRNA-seq) analysis.    

**🔬 Core Analysis Libraries**    
*<u>scanpy</u>*        
A comprehensive toolkit for analyzing single-cell RNA-seq data, including preprocessing, clustering, and visualization.    
Decision: Chosen as the primary framework because it efficiently handles large-scale single-cell datasets and integrates multiple analysis steps in one pipeline. 
Biological relevance: Enables identification of cell populations and gene expression patterns across individual cells.    
*<u>celltypist</u>*         
A machine learning-based tool for automated cell type annotation.    
Decision: Used to avoid manual annotation and leverage pre-trained reference models for faster, standardized results.    
Biological relevance: Assigns biological identities (e.g., immune cell types) to clusters derived from expression data.    
**📊 Data Handling & Computation**    
*<u>pandas</u>*    
Used for structured data manipulation (e.g., tables, metadata handling).    
Decision: Provides flexible and intuitive data structures for handling annotated biological datasets.    
*<u>numpy</u>*    
Supports numerical computations and array operations.    
Decision: Essential for efficient handling of large gene expression matrices.    
*<u>os</u>*    
Handles file paths and directory operations.    
Decision: Ensures portability and easy management of input/output files across systems.    
**📈 Visualization Libraries**    
*<u>matplotlib</u>*    
A foundational plotting library for creating static visualizations.    
Decision: Used for customizable, publication-quality plots.    
*<u>seaborn</u>*    
Built on top of matplotlib, providing more aesthetically pleasing and statistically informative visualizations.    
Decision: Simplifies the creation of complex plots like heatmaps and distributions.    
**⚙️ Utility**    
*<u>warnings</u>*    
Controls and suppresses unnecessary warning messages.    
Decision: Keeps output clean and focused, especially when working with large pipelines that may generate non-critical warnings.    

```bash
# Load the AnnData object from a specific file path
# This file likely contains pre-calculated PCA and UMAP coordinates from a Seurat-to-Scanpy conversion
GSE183276_raw_adata = ad.read_h5ad("D:/Bidya Work/single/GSE183276/output/03_GSE183276_seurat_pca_umap.h5ad")
GSE183276_raw_adata
```
<img width="1099" height="212" alt="image" src="https://github.com/user-attachments/assets/1a5c0df6-2043-4602-b13f-d92ffe2a0bb6" />    

By using this step, I loaded the single-cell RNA-seq dataset into an AnnData object using Scanpy. The dataset is stored in .h5ad format, which efficiently contains gene expression data along with metadata and analysis results. The file already includes precomputed PCA and UMAP embeddings, allowing immediate visualization and downstream analysis without repeating computationally intensive steps.

```bash
# Importing the batch-corrected single-cell RNA-seq data (GSE183276) 
# The 'harmony_corrected' suffix indicates that batch effects between samples have been removed.
GSE183276_harmony_adata = ad.read_h5ad("D:/Bidya Work/single/GSE183276/output/04_GSE183276_harmony_corrected.h5ad")
GSE183276_harmony_adata
```
<img width="1100" height="231" alt="image" src="https://github.com/user-attachments/assets/82adc0b8-1773-45ce-ae46-3eded9f8635d" />    

This step loads the batch-corrected single-cell RNA-seq dataset into an AnnData object. The data has been processed using Harmony, a method for correcting batch effects across multiple samples. This dataset has been batch-corrected using Harmony to remove technical variation between samples on R previously (refer the pipeline). This allows cells to cluster based on true biological similarities rather than experimental differences, improving the accuracy of downstream analysis.

```bash
import pandas as pd
import numpy as np

# 1. Load the external UMAP coordinates ( exported from Seurat in R)
# Setting index_col=0 ensures the cell barcodes are used as the index for matching
umap_df = pd.read_csv(
    "D:/Bidya Work/single/GSE183276/GSE183276_umap_coordinates.csv",
    index_col=0
)

# 2. Verify dimensions: ensure the number of cells (5405) and coordinates (2) match expectations
print(umap_df.shape)   # should be (5405, 2)

# 3. Align and inject coordinates into the AnnData object
# We use .loc[GSE183276_harmony_adata.obs_names] to ensure the UMAP rows 
# match the exact order of cells in our AnnData object.
GSE183276_harmony_adata.obsm['X_umap'] = umap_df.loc[
    GSE183276_harmony_adata.obs_names
].values

# 4. Visualize the data using the newly injected Seurat UMAP coordinates
sc.pl.umap(
    GSE183276_harmony_adata,
    color='seurat_clusters',
    show=True
)
```
<img width="763" height="454" alt="image" src="https://github.com/user-attachments/assets/1db538a4-fd43-47a2-a10a-a75112c3a881" />   
I took 5,405 cells → reduced them into 2D → and found 20 biologically meaningful groups.     
This visualization reveals the cellular heterogeneity within the dataset. Distinct clusters likely correspond to different:    
Cell types (e.g., immune cells, epithelial cells)    
Cellular states (e.g., activated vs resting cells)    
Biological conditions or subpopulations    
Each point corresponds to a single cell profiled using single-cell RNA sequencing. Cells with similar gene expression profiles cluster together in the 2D space, Distinct colors represent different Seurat clusters, A total of 20 clusters are observed, indicating strong cellular heterogeneity within the dataset.    

Well-separated clusters indicate strong transcriptional differences, while closely positioned clusters may represent related cell types or transitional states.
UMAP coordinates generated in Seurat (R) were imported and added to the Scanpy AnnData object. The coordinates were aligned using cell barcodes to ensure correct mapping between datasets. Instead of recomputing UMAP in Scanpy, previously generated coordinates were reused to maintain consistency with earlier analysis. UMAP projection reduces high-dimensional gene expression data into a 2D space, allowing visualization of cellular relationships. UMAP was stored in .obsm['X_umap'], following Scanpy’s standard data structure. Proper alignment ensures that  biologically meaningful clusters (e.g., cell types or states) are accurately represented.

```bash
#In an AnnData object, the .X attribute is the "active" matrix. Most functions look at .X by default.
#The .layers slot: Usually acts as a container for different versions of your data (e.g., raw, normalized, logcounts).
#The .copy() requirement: Without .copy(), Python might just create a "view" (a shortcut). If you then modify .X, you might accidentally modify your saved logcounts layer too. Copying keeps them separate and safe.
GSE183276_harmony_adata.X = GSE183276_harmony_adata.layers['logcounts'].copy()
```
In this step, I set the main data matrix to use the log-normalized gene expression values. The AnnData object can store different versions of the data (raw, normalized, etc.). I selected the log-normalized version because it is more stable and suitable for analysis.    
I also used .copy() to make sure this change does not accidentally affect the original stored data.    

```bash
#You are using the batch-corrected Harmony embeddings to group cells and then applying a pre-trained model to "name" those groups.
import scanpy as sc
import celltypist

# 1. Define cell-to-cell relationships based on corrected data
# We use 'HARMONY' instead of the default PCA to ensure the graph as it reflects the batch-corrected space.
sc.pp.neighbors(GSE183276_harmony_adata, use_rep='HARMONY')

# 2. Group cells into high-resolution clusters
# A resolution of 10 is very high (over-clustering); this creates many small 
# sub-groups to help CellTypist make more precise 'majority voting' decisions.
sc.tl.leiden(GSE183276_harmony_adata, resolution=10, key_added='celltypist_clusters')

# 3. Predict cell types using the CellTypist model
# We use a specific Kidney model and 'majority_voting', which assigns a single 
# cell type to an entire cluster based on the most frequent prediction.
ct_pred = celltypist.annotate(
    GSE183276_harmony_adata,
    model='D:/Bidya Work/single/GSE183276/Adult_Human_Kidney.pkl',
    majority_voting=True,
    over_clustering='celltypist_clusters'
)

# 4. View the results
print(ct_pred.predicted_labels.head())
```
<img width="1103" height="152" alt="image" src="https://github.com/user-attachments/assets/ae049fc8-5824-412b-9f78-50a52271c62b" />    
In this step, I first grouped similar cells together using the batch-corrected Harmony data. This ensures that cells are grouped based on real biological similarity, not technical differences. Then, I split the cells into many small groups to make the predictions more accurate.    
Finally, I used a pre-trained model (CellTypist) to assign a cell type label (like kidney cell types) to each group based on the most common prediction in that group.

<img width="601" height="232" alt="image" src="https://github.com/user-attachments/assets/d359544a-35e9-4010-a5fd-2260ac844710" />      
This table shows the final kidney cell type assigned to each cell based on cluster-level majority voting.

```bash
# Count and display the total number of clusters generated by the high-resolution Leiden clustering.
# This helps verify if the 'over-clustering' worked as intended for CellTypist.
print(f"Number of clusters created: {len(GSE183276_harmony_adata.obs['celltypist_clusters'].unique())}")
```
This created 81 clusters as this step quantifies the number of Leiden clusters generated after applying high-resolution clustering. The goal is to intentionally over-cluster the data prior to automated cell type annotation using CellTypist ensuring rare and subtle cell states are preserved and improving the accuracy and robustness. 

```bash
# Transfer the 'majority voting' predictions from the CellTypist result to the AnnData object.
# We convert to string (.astype(str)) to ensure compatibility with Scanpy's plotting functions.
GSE183276_harmony_adata.obs['majority_voting'] = (
    ct_pred.predicted_labels['majority_voting']
    .astype(str)
)
```
I tranferred the CellTypist majority-vote predictions into the AnnData object as the final per-cell annotations, providing biologically interpretable labels for downstream visualization and analysis in Scanpy. 

```bash
# Calculate and display the top 5 most frequent cell types identified by majority voting.
# This provides a quick census of the major cell populations in the dataset
print(GSE183276_harmony_adata.obs['majority_voting'].value_counts().head())
```
This step provides a high-level census of predicted cell identities, summarizing the major cellular composition of the dataset after CellTypist annotation.    
<img width="408" height="130" alt="image" src="https://github.com/user-attachments/assets/71cf29e7-8281-4ad1-b03f-caa83399e6f1" />    

```bash
import os
import matplotlib.pyplot as plt
import scanpy as sc

# Create output directory if it doesn't exist
out_dir = "D:/Bidya Work/single/GSE183276/plots/celltypist_plots/"
os.makedirs(out_dir, exist_ok=True)

# Iterate through each unique group in the 'Condition' column
for cond in GSE183276_harmony_adata.obs['Condition'].unique():

    # Generate a UMAP plot for a subset of the data belonging to the current condition
    adata_subset = GSE183276_harmony_adata[
        GSE183276_harmony_adata.obs['Condition'] == cond
    ]

    fig = sc.pl.umap(
        adata_subset,
        color='majority_voting',   # Color cells by CellTypist predicted labels
        title=f'Condition: {cond}',   # Dynamically set the title
        legend_loc='right margin',
        show=False,         # Don't show immediately so we can save it first
        return_fig=True     # Capture the figure object for saving
    )

    # Save the figure for each condition
    fig.savefig(
        os.path.join(out_dir, f"umap_celltypes_ct_{cond}.png"),
        dpi=300,
        bbox_inches='tight'
    )

    # Show the plot
    plt.show()

    # Close figure to free memory
    plt.close(fig)
```
To investigate how predicted cell types are distributed across experimental conditions, I generated UMAP embeddings separately for each condition using CellTypist majority-voted labels. For each condition, the dataset was subset based on the Condition metadata field, and cell identities were visualized in UMAP space using the majority_voting annotations as the color key. This approach enables direct comparison of cellular composition and structure across biological conditions. 

<img width="1114" height="611" alt="image" src="https://github.com/user-attachments/assets/f741e0c4-8c17-4cb6-9171-c826a46186a8" />    

<img width="1129" height="612" alt="image" src="https://github.com/user-attachments/assets/25493e74-a564-4989-b4f4-ceb5af48a324" />    

<img width="1132" height="617" alt="image" src="https://github.com/user-attachments/assets/e1ec8fe9-486a-42fd-a5b7-9bca93bd9ba7" />    

<img width="1105" height="618" alt="image" src="https://github.com/user-attachments/assets/e0a11d94-9909-4e5c-9201-a3f9cee0a17d" />    

Across all four conditions, major cell populations (e.g., tubular, endothelial, and immune cells) form well-separated and consistent clusters in UMAP space. This indicates that the integration and annotation pipeline preserved biologically meaningful structure across datasets. Immune populations (e.g., B cells, dendritic cells) appear: More prominent in AKI and DKD, Less abundant in Healthy. This pattern is consistent with immune activation or infiltration during kidney injury and disease progression.  
The UMAP analysis reveals that kidney disease conditions (AKI and DKD) are associated with shifts in cellular composition particularly increased immune cell representation and changes in epithelial cell states, while the overall cellular architecture remains conserved.

NOTE: In this workflow, both Harmony-integrated and Scanpy-standard AnnData objects were explored. However, all downstream analyses (UMAP visualization, CellTypist annotation, and compositional analysis) are fully compatible with the Harmony-integrated object. The use of a separate raw_adata object was primarily for format standardization and pipeline testing, and was not required for biological interpretation.

```bash
# Verify that the variable is a valid AnnData object.
# This confirms the data was loaded correctly into the Scanpy-compatible format.
print(type(GSE183276_harmony_adata))
```
Verified that the dataset is a valid AnnData object to ensure compatibility with Scanpy-based single-cell analysis workflows. <img width="306" height="33" alt="image" src="https://github.com/user-attachments/assets/c2a94b2d-322d-4135-8558-71bb9d56c798" />    

```bash
##This block of code is a troubleshooting and refinement step. You are renaming your coordinates to follow Scanpy’s naming conventions 
##and then re-running the clustering to ensure the automated annotation tool (CellTypist) has the best possible input.
import scanpy as sc
import celltypist

# 1. Standardize coordinate naming
# Scanpy functions look for 'X_pca' by default. This clones the existing 'PCA' 
# slot into 'X_pca' to ensure compatibility with downstream tools.
GSE183276_raw_adata.obsm['X_pca'] = GSE183276_raw_adata.obsm['PCA']

# 2. Re-calculate clusters for annotation
# We build a new neighborhood graph from the PCA and use high-resolution 
# Leiden clustering (over-clustering) to create fine-grained groups for prediction.
sc.pp.neighbors(GSE183276_raw_adata, use_rep='X_pca')
sc.tl.leiden(GSE183276_raw_adata, resolution=10, key_added='ct_clusters')

# 3. Perform automated cell type annotation
# CellTypist uses the model to predict labels and then 'smoothes' those 
# labels by applying a majority vote within each 'ct_cluster'.
GSE183276_predictions3 = celltypist.annotate(
    GSE183276_raw_adata, 
    model="D:/Bidya Work/single/GSE183276/Adult_Human_Kidney.pkl",
    majority_voting=True,
    over_clustering='ct_clusters'
)

# 4. Finalize the data object
# Cast the predictions to strings and save them into the metadata (.obs) 
# for easy plotting and data export.
GSE183276_raw_adata.obs["majority_voting"] = GSE183276_predictions3.predicted_labels["majority_voting"].astype(str)

print("Annotation complete!")
```
<img width="568" height="169" alt="image" src="https://github.com/user-attachments/assets/a083bb38-8ceb-4a8a-af08-311ac70d497d" />    

The PCA renaming step was necessary to make sure Scanpy functions work correctly without errors or inconsistencies. I used high-resolution clustering because it creates smaller, more detailed groups of cells. This helps separate closely related or rare cell populations, which might get merged if clustering is too coarse. Using CellTypist allowed me to automatically assign biological cell type labels instead of manually identifying them, which saves time and improves consistency. The majority voting step helps reduce noise from individual cell predictions and produces more stable and reliable labels.   

The initial CellTypist annotation was performed on earlier clustering results. However, since annotation quality depends on how well cells are grouped, I re-ran the process after improving the clustering. In this step, I used high-resolution Leiden clustering to create smaller and more homogeneous clusters. This improves the accuracy of CellTypist predictions, especially when using majority voting, which relies on cluster-level agreement. As a result, the second annotation provides more reliable and biologically meaningful cell type labels.I basically refined my clustering to create detailed groups of cells and then used CellTypist to assign reliable cell type labels, making the dataset ready for meaningful biological analysis.    

```bash
##This code performs a compositional analysis, which is a fancy way of saying it calculates which cell types are increasing or decreasing between your experimental conditions.
import matplotlib.pyplot as plt
import seaborn as sns

# Step 1: Quantify the number of cells for each cell type within each condition
# This creates a summary table of raw counts
n_cells_condition = (
    GSE183276_raw_adata.obs
    .groupby(["Condition", "majority_voting"])
    .size()
    .reset_index(name="count")
)

# Step 2: Normalize counts to percentages (Proportions)
# This allows for fair comparison even if one condition has more total cells than another
n_cells_condition["total"] = n_cells_condition.groupby("Condition")["count"].transform("sum")
n_cells_condition["proportion"] = (n_cells_condition["count"] / n_cells_condition["total"]) * 100

# Step 3: Sort cell types by abundance
# We find the average percentage of each cell type across all conditions to order the plot
avg_proportions = (
    n_cells_condition
    .groupby("majority_voting")["proportion"]
    .mean()
    .sort_values(ascending=False)
)
ordered_celltypes = avg_proportions.index.tolist()

# Step 4: Generate the Bar Plot
# 'dodge=True' places the condition bars side-by-side for easy comparison
plt.figure(figsize=(40, 10))
ax = sns.barplot(
    data=n_cells_condition,
    x="majority_voting",
    y="proportion",
    hue="Condition",
    order=ordered_celltypes,
    dodge=True,
    width=0.9
)

# Step 5: Add percentage labels above the bars
# This makes the exact values readable without looking at the axis
for container in ax.containers:
    ax.bar_label(container, fmt='%.1f%%', label_type='edge', padding=2, fontsize=10, rotation = 90)

# Customize axes
plt.xticks(rotation=45, ha="right", fontsize=12, fontweight="bold")
plt.ylabel("Proportion (%)", fontsize=12, fontweight="bold")
plt.xlabel("Cell Type", fontsize=12, fontweight="bold")
plt.title("Cell Type Proportions", fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig("D:/Bidya Work/single/GSE183276/plots/celltypist_plots/Cellproportions_barplot.png", dpi=300, bbox_inches='tight')
plt.show()

n_cells_condition.to_csv("Cellproportions.csv", index=False)
```
I performed a compositional analysis to compare how different cell types are distributed across experimental conditions (Healthy, AKI, DKD, HCKD).    
First, I counted the number of cells belonging to each predicted cell type (majority_voting) within each condition. This gave a raw count of how many cells of each type are present per condition. Next, I converted these counts into percentages (proportions) so that conditions with different total cell numbers could be fairly compared.    
After that, I calculated the average proportion of each cell type across all conditions to determine an order for plotting. This helps visualize the most abundant cell types first. Finally, I created a grouped bar plot to compare cell type proportions across conditions and saved both the figure and the processed data table.   

<img width="1919" height="484" alt="image" src="https://github.com/user-attachments/assets/84249c3b-90b9-41b0-bf91-9647a9df6053" />    

<img width="499" height="652" alt="image" src="https://github.com/user-attachments/assets/07e62053-ba8a-465a-ba2b-9ff6161c66c5" />

Across the four conditions (Healthy, AKI, DKD, and HCKD), clear differences in cell type composition are observed. In Healthy samples, the distribution of major kidney cell populations is more balanced, with relatively stable proportions across epithelial and vascular compartments. In contrast, disease conditions (AKI and DKD) show noticeable shifts in cellular composition, particularly in epithelial populations such as proximal tubule and related tubular subtypes, which display altered proportions compared to Healthy. Immune-related cell types (such as B cells and dendritic cells) appear more prominent in AKI and DKD, suggesting increased immune involvement in disease states. Endothelial and vascular-associated populations also show variability across conditions, indicating possible structural or functional changes in the kidney microenvironment. HCKD generally shows an intermediate pattern between Healthy and disease conditions, suggesting a partial transition in cellular composition.

Overall, the results indicate that kidney disease conditions are associated with clear compositional shifts in key cell populations especially increased immune representation and altered epithelial cell proportions—highlighting disease-driven remodeling of the renal cellular landscape.    

```bash
# Create a new combined label by merging the predicted cell type with the condition.
# Resulting format: "CellType_Condition" (e.g., "Podocyte_Healthy").
# This is useful for performing differential expression between conditions within the same cell type.
GSE183276_raw_adata.obs['majority_voting'] = (
    GSE183276_raw_adata.obs['majority_voting'].astype(str)
    + "_" +
    GSE183276_raw_adata.obs['Condition'].astype(str)
)
```
I merged cell type and condition into a single label so I can directly compare the same cell type across different disease states.     

```bash
def find_unique_markers(
    adata, 
    groupby='cell_type', 
    method='wilcoxon', 
    pval_threshold=0.05, 
    logfc_threshold=0.25,
    top_n=3,
    min_cells_per_group=2
):
    """
    Mimics Seurat's FindAllMarkers + filters for unique DEGs per cluster.
    
    Parameters:
        adata : AnnData object
        groupby : column in adata.obs to group cells (e.g., clusters, cell_type)
        method : DEG test method ('wilcoxon', 't-test', 'logreg')
        pval_threshold : adjusted p-value threshold for significance
        logfc_threshold : log fold change threshold for filtering

    Returns:
        unique_degs_df : DataFrame of DEGs unique to each group
    """

   # 1. Quality Control: Remove groups with too few cells to ensure statistical power
    group_counts = adata.obs[groupby].value_counts()
    valid_groups = group_counts[group_counts >= min_cells_per_group].index.tolist()

  # 2. Subset to valid data: Use .copy() to avoid modifying the original AnnData
    adata_filtered = adata[adata.obs[groupby].isin(valid_groups)].copy()
    
   # 3. Differential Expression Analysis: Compute rankings for all genes across groups
    sc.tl.rank_genes_groups(adata_filtered, groupby=groupby, method=method)
    
   # 4. Extract Results: Convert the structured numpy array into a tidy Pandas DataFrame
    all_degs = sc.get.rank_genes_groups_df(adata_filtered, group=None)
    
   # 5. Significance Filtering: Apply adjusted p-value and log-fold change cutoffs
    filtered_degs = all_degs[
        (all_degs['pvals_adj'] < pval_threshold) &
        (abs(all_degs['logfoldchanges']) > logfc_threshold)
    ]

  # 6. Exclusivity Check: Identify genes that passed filtering in ONLY ONE group.
    # This prevents 'Pan-marker' genes (like general immune markers) from appearing.
    unique_genes = (
        filtered_degs.groupby('names')['group']
        .nunique()
        .reset_index()
        .query('group == 1')['names']
        .tolist()
    )

   # 7. Final Subset: Keep only the genes identified as truly unique to a cluster
    unique_degs_df = filtered_degs[filtered_degs['names'].isin(unique_genes)].copy()


  # 8. Representative Selection: Sort by effect size (logFC) and pick the top N per group
    top_unique_degs_df = (
        unique_degs_df
        .sort_values(['group', 'logfoldchanges'], ascending=[True, False])
        .groupby('group')
        .head(top_n)
        .reset_index(drop=True)
    )

    return top_unique_degs_df
```
I created a custom function to find marker genes for each cell group (such as cell types or clusters) using differential gene expression analysis. First, I filtered out groups with too few cells to make sure the results were statistically reliable. Then, I ran Scanpy’s differential expression test to identify genes that are differentially expressed across groups. After that, I applied filters based on statistical significance (adjusted p-value) and expression strength (log fold change) to keep only meaningful genes. To make the results more specific, I further removed genes that appeared as markers in more than one group, keeping only truly unique markers for each cell type. Finally, I selected the top few strongest markers per group to make the output easier to interpret. This was done to ensure that the final marker list is biologically meaningful, specific to each cell population, and easier to use for understanding cell identity and validating annotations.     

```bash
# 1. Identify the top 3 highly specific markers for each cell type
# Uses the Wilcoxon rank-sum test to find genes that are significant (p < 0.05) 
# and exclusive to only one group.unique_markers = find_unique_markers(
    GSE183276_raw_adata, 
    groupby='majority_voting',
    method='wilcoxon',
    pval_threshold=0.05, 
    logfc_threshold=0.25,
    top_n=3, 
    min_cells_per_group=2
)
# 2. Clean up group names using Regular Expressions (Regex)
# This removes redundant suffixes that might have been created during concatenation 
# (e.g., changing "Podocyte_Healthy_Healthy" back to "Podocyte_Healthy").
unique_markers['group'] = (
    unique_markers['group']
    .str.replace(r'_(HCKD|Healthy|AKI|DKD)_\1$', r'_\1', regex=True)
)

# 3. Validation and Export
# Print the top results, verify unique group names, and save to a CSV file.
print(unique_markers.head())
unique_markers['group'].drop_duplicates()

unique_markers.to_csv("Unique_cluster_markers.csv", index=False)
```
I identified the top 3 most specific marker genes for each cell type using a differential gene expression approach based on the Wilcoxon rank-sum test. The analysis was performed on the CellTypist-annotated labels (majority_voting) to find genes that are significantly enriched in one cell type compared to all others, using thresholds for statistical significance (adjusted p-value < 0.05) and expression strength (log fold change > 0.25). I also ensured that each gene is unique to a single cell type so that only truly specific markers are retained. After generating the results, I cleaned the group names using pattern matching to remove any duplicated condition suffixes that may have been introduced during label creation, making the group names consistent and readable. Finally, I reviewed the output and saved the final list of unique marker genes to a CSV file for further analysis and validation.

<img width="542" height="681" alt="image" src="https://github.com/user-attachments/assets/3cc71771-3f76-4b0d-a4d5-932c86a95aba" />    
This table shows the top 3 most specific marker genes identified for each cell type across different conditions (Healthy, AKI, DKD, HCKD). Each row represents a gene that is strongly and uniquely expressed in a particular cell population.    

```bash
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

# 1. Subset Data: Remove specific populations (e.g., Classical monocytes) 
# that might be redundant or distracting in this specific visualization.
excluded_prefixes = ["Classical monocytes"]
mask = ~GSE183276_raw_adata.obs['majority_voting'].str.startswith(tuple(excluded_prefixes))
adata_to_plot = GSE183276_raw_adata[mask].copy()

# 2. Filter Gene List: Ensure all selected markers actually exist in the matrix.
genes_to_plot = unique_markers["names"].tolist()  # top 3 per group
genes_to_plot = [g for g in genes_to_plot if g in adata_to_plot.var_names]

# 3. Data Transformation: Convert sparse matrix to dense array for calculation.
# This allows us to perform group-wise summation of gene expression.
X = adata_to_plot[:, genes_to_plot].X
if hasattr(X, "toarray"):  # sparse -> dense
    X = X.toarray()

# 4. Aggregation: Create a summary table of gene expression per cell-type group.
expr_df = pd.DataFrame(X, index=adata_to_plot.obs['majority_voting'], columns=genes_to_plot)
expr_per_group = expr_df.groupby(expr_df.index).sum()  # rows = groups, columns = genes

# 5. Stringent Filtering: Keep only genes with substantial total signal (sum >= 30).
# This removes "weak" markers that might be significant but aren't visually clear.
genes_filtered = expr_per_group.columns[(expr_per_group.max(axis=0) >= 30)].tolist()
expr_per_group = expr_per_group[genes_filtered]

# 6. Group Filtering: Ensure each group shown has a minimum total expression signal.
valid_groups = expr_per_group.index[(expr_per_group.sum(axis=1) > 20)].tolist()
adata_to_plot = adata_to_plot[adata_to_plot.obs['majority_voting'].isin(valid_groups)].copy()

# 7. Dot Plot Visualization:
# standard_scale='var' normalizes color from 0 to 1 for each gene, making it
# easier to see where each gene is uniquely "turned on."
sc.pl.dotplot(
    adata_to_plot,
    var_names=genes_filtered,
    groupby='majority_voting',
    standard_scale='var',
    show=False,
    figsize=(20, 20),
    dendrogram=False
)
plt.savefig("D:/Bidya Work/single/GSE183276/plots/celltypist_plots/Cluster_Markergenes.png", dpi=300, bbox_inches='tight')
plt.show()
```
I created a dot plot to visualize how strongly the selected marker genes are expressed across different CellTypist-annotated cell types.
First, I removed a few cell populations (such as Classical monocytes) that were not relevant for this specific visualization and could make the plot harder to interpret. Then, I selected the list of top marker genes identified earlier and ensured that only genes present in the dataset were used. Next, I extracted the gene expression matrix for these markers and converted it into a usable format for analysis. I then grouped the expression values by cell type and calculated the total expression of each gene within each group. This helped summarize how strongly each marker gene is expressed in each cell population.    
To improve clarity, I filtered out genes with weak overall expression and also removed cell groups with very low total expression. This ensured that only meaningful and visually interpretable signals were included in the final plot. Finally, I used a dot plot to visualize the results, where each dot represents the expression level of a gene in a specific cell type. The values were normalized to make it easier to compare patterns across genes and cell types.    

<img width="949" height="911" alt="image" src="https://github.com/user-attachments/assets/4cead082-d64f-4d8f-89fe-26743f91afe7" />

This dot plot visualizes the expression of selected top marker genes across all CellTypist-annotated cell types and conditions. Each row represents a cell type–condition combination (e.g., PT_AKI, B cells_DKD), and each column represents a marker gene identified from the previous differential expression analysis. The size of each dot indicates the fraction of cells expressing a gene in that group, while the color intensity represents the average expression level.

```bash
GSE183276_raw_adata.write("D:/Bidya Work/single/GSE183276/output/GSE183276_celltypist.h5ad")
```
I saved the final processed AnnData object (GSE183276_raw_adata) as an .h5ad file using Scanpy’s built-in write function. This file contains all important results from the analysis pipeline, including normalized data, PCA embeddings, clustering results, CellTypist annotations, and metadata such as cell type labels and conditions. I saved the final processed dataset so I can reload it later without repeating the entire analysis and easily continue working from this point. 

Following cell type annotation and marker gene visualization, downstream pathway enrichment analysis was performed to better understand the biological functions associated with the identified marker genes. This step helped me to interpret how different cell populations contribute to disease-associated molecular pathways and cellular processes.  
The marker genes identified from the differential expression analysis were therefore used for pathway and functional enrichment analysis in R.    
[Click to view script file](https://github.com/Bidya122/scRNA-seq-Analysis_kidney_disease/blob/main/03_scRNAseq_Bidya.Rmd)

# Installation of other packages for this workflow
```bash
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("reactome.db")
```
I did not have to download each and every package separately, since the main stack of tools required for the analysis had already been installed during the initial setup of the workflow. Only the remaining packages needed for downstream analysis were installed here.   

# Setting working directories and the needed libraries
```bash
setwd("D:/Bidya Work/single/GSE183276")
inputDir <- "D:/Bidya Work/single/GSE183276/input"
outputDir <- "D:/Bidya Work/single/GSE183276/output"
plotDir <- "D:/Bidya Work/single/GSE183276/plots"
library(zellkonverter)
library(SingleCellExperiment)
library(Seurat)
library(MAST)
library(tidyverse)
```

# Uploading the required files/datasets
```
sce <- readH5AD("D:/Bidya Work/single/GSE183276/output/GSE183276_celltypist.h5ad") ## readH5AD() converts it into a SingleCellExperiment (SCE) object in R

####you will see some messages like:
#ℹ Using stored X_name value 'counts'
#<sys>:0: FutureWarning: Use varm (e.g. `k in adata.varm` or `adata.varm.keys() | {'u'}`) instead of AnnData.varm_keys, AnnData.varm_keys is deprecated and will be removed in the future.
#<sys>:0: FutureWarning: Use obsm (e.g. `k in adata.obsm` or `adata.obsm.keys() | {'u'}`) instead of AnnData.obsm_keys, AnnData.obsm_keys is deprecated and will be removed in the future.

seurat_obj <- as.Seurat(sce, counts = "counts", data = "logcounts")  # Convert the SingleCellExperiment object into a Seurat object

saveRDS(seurat_obj, "D:/Bidya Work/single/GSE183276/output/GSE183276_with_celltypist.rds") # Save the Seurat object as an .rds file. saveRDS() stores the object in compressed R format, so you can reload it later without repeating conversion steps
```
I loaded the processed .h5ad file generated from the Scanpy workflow into R using readH5AD(), which converted the dataset into a SingleCellExperiment object. I then converted this object into a Seurat object using as.Seurat() so the dataset could be used for downstream analysis and visualization in R. Finally, I saved the Seurat object as an .rds file using saveRDS() to allow easy reloading of the processed dataset in future sessions without repeating the conversion steps.

# Cleaning and standardizing Annotations
```bash
colnames(seurat_obj@meta.data)
# Create a contingency table between:
# rows   = disease condition
# columns = predicted cell types from CellTypist
table(seurat_obj$Condition, seurat_obj$majority_voting)  

#This cleans cell type names so disease status and cell identity are separated properly.
seurat_obj$majority_voting <- gsub("_(AKI|DKD|HCKD|Healthy)$", "", seurat_obj$majority_voting) # Remove disease suffixes from CellTypist labels

table(seurat_obj$Condition, seurat_obj$majority_voting)  # Recheck the table after cleaning

# If suffixes are still present in some labels, run the gsub() command one more time.
```
<img width="1258" height="278" alt="image" src="https://github.com/user-attachments/assets/e6aed588-e00e-412d-a13c-5a1a06cb211a" />    

<img width="1194" height="632" alt="image" src="https://github.com/user-attachments/assets/6c095239-c140-4f5c-afd6-f5bc7bc8c35a" />  


I first examined the metadata columns in the Seurat object and generated a contingency table to compare disease conditions with the CellTypist-predicted cell type labels. During this step, I observed that some cell type labels contained disease-specific suffixes such as _AKI, _DKD, _HCKD, and _Healthy. To ensure that cell identity and disease condition remained separate metadata fields, I cleaned the majority_voting labels using gsub() to remove these suffixes from the cell type names as shown. The contingency table was then rechecked to confirm that the labels had been cleaned correctly before proceeding with downstream analysis. 

# Preparaing cell type lables
```bash
DefaultAssay(seurat_obj) <- "originalexp"  ## The default assay tells Seurat which expression matrix should be used for downstream functions by default.
# Here, we are choosing the original expression data.

celltypes <- sort(unique(seurat_obj$majority_voting)) # Extract all unique cell type names from the metadata column "majority_voting"
celltypes

celltypes <- gsub("-", "_", celltypes)  # Replace hyphens (-) with underscores (_)
celltypes
```
<img width="1145" height="75" alt="image" src="https://github.com/user-attachments/assets/c599c5cc-a30b-42c0-8093-26ff9f583906" />

I set the default assay of the Seurat object to originalexp so that downstream analyses would use the original expression matrix by default. Next, I extracted all unique CellTypist-predicted cell type labels from the majority_voting metadata column and sorted them for easier handling in downstream analyses. Finally, I replaced hyphens (-) with underscores (_) in the cell type names to maintain consistent formatting and avoid potential issues during file naming, plotting, or automated analysis steps.    

# Cell Type-Specific Differential Expression Analysis Using MAST 
```bash
#DKD vs Healthy: What changes in diabetic kidney disease?
#AKI vs Healthy: What changes during acute kidney injury?
#HCKD vs Healthy: What changes in chronic kidney disease?
#DKD vs AKI: What makes diabetic kidney disease different from acute kidney injury?


run_mast_per_celltype_DKD_vs_Healthy <- function(    ##unction to perform cell type–specific differential expression analysis
# between DKD and Healthy samples using the MAST statistical method
  seu,
  celltype_col = "majority_voting",     # apna annotation column
  condition_col = "Condition",   
  outdir = "D:/Bidya Work/single/GSE183276/output/DEG_MAST_DKD_vs_Healthy",
  min_cells_per_group = 20   # Minimum number of cells required in each condition to perform DEG analysis reliably
) {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  meta <- seu@meta.data   # Extract metadata from Seurat object
  celltypes <- sort (unique(meta[[celltype_col]]))    # Get all unique cell types present in the dataset and sort them alphabetically

  for (ct in celltypes) {     # Loop through each cell type one by one; Take each value from celltypes one by one, and temporarily store it in a variable called ct

    cat("\n▶ Processing cell type:", ct, "\n") ##cat() means concatenate and print.

    cells_ct <- rownames(meta)[meta[[celltype_col]] == ct]      # Identify cells belonging only to the current cell type
    cond_tab <- table(meta[cells_ct, condition_col]) # counts ALL conditions that exist within that cell type.
    print(cond_tab)
    
    ##now we are seeing if the two required conditions by us are present in all the cell types or not

    # Check whether both groups exist ---If either DKD or Healthy cells are absent, skip this cell type because comparison is impossible
    if (!all(c("DKD", "Healthy") %in% names(cond_tab))) {
      cat("  ❌ DKD Or Healthy missing → skipping\n")
      next
    }

    # # ---- Minimum cell count check --Skip analysis if either condition has fewer cells than the defined threshold. 
    if (any(cond_tab[c("DKD", "Healthy")] < min_cells_per_group)) {
      cat("  ❌ Not enough cells → skipping\n") 
      next
    }

   
    seu_ct <- subset(seu, cells = cells_ct)      # Create a new Seurat object containing only cells from the current cell type
    DefaultAssay(seu_ct) <- "originalexp"     # Set assay to original expression matrix
    seu_ct <- NormalizeData(seu_ct, verbose = FALSE)       # Normalize expression data Required before running DEG analysis

    Idents(seu_ct) <- condition_col    # Set cell identities based on condition. This tells Seurat which groups to compare

    # # ---- Run MAST differential expression analysis ----
     # ident.1 = DKD group #ident.2 = Healthy group
    deg <- FindMarkers(
      seu_ct,
      ident.1 = "DKD",
      ident.2 = "Healthy",
      test.use = "MAST",
      logfc.threshold = 0,   # -> tests all genes regardless of fold change
      min.pct = 0.1,        # -> gene must be expressed in at least 10% of cells
      latent.vars = c("nCount_RNA", "percent.mt")   # latent.vars-> adjusts for technical confounders:
    #    nCount_RNA = sequencing depth
    #    percent.mt = mitochondrial gene percentage

    )

    deg$gene <- rownames(deg)    # Add gene names as a separate column

    # 🔥 CRITICAL FIX (THIS LINE SOLVES YOUR ERROR)
    ct_clean <- gsub("[^a-zA-Z0-9]", "_", ct)    # Prevents file path and naming errors Clean cell type names for safe file creation
    # Replaces spaces/special characters with underscores
    # Prevents file path and naming errors

    out_file <- file.path(outdir, paste0("DEG_", ct_clean, "_DKD_vs_Healthy.csv"))

    write.csv(deg, out_file, row.names = FALSE)

    cat("  ✔ Saved:", out_file, "\n")
  }
}
```

Among the different disease comparisons available in the dataset, I initially focused on the DKD vs Healthy comparison to identify transcriptional changes specifically associated with diabetic kidney disease. Comparing diseased samples directly against healthy controls provides a clearer baseline for understanding disease-related molecular alterations within each cell type.
Before performing differential expression analysis, several quality checks and filtering steps were applied to ensure reliable statistical results. For each CellTypist-annotated cell population, I first verified that both DKD and Healthy cells were present. Cell types lacking one of the comparison groups were excluded since meaningful differential analysis would not be possible. I also applied a minimum cell count threshold of 20 minimum cells per group to ensure that each condition contained enough cells for statistically robust testing.
Only cell populations that passed these criteria were included in the downstream MAST-based differential expression analysis. The resulting differentially expressed genes were then saved separately for each cell type for further functional and pathway enrichment analysis. Differential expression analysis was performed using the MAST statistical method through Seurat’s FindMarkers() function. Although methods such as Wilcoxon rank-sum testing can also be used for differential expression analysis, MAST was selected because it is specifically designed for single-cell RNA-seq data.
Single-cell datasets are typically sparse and contain a large number of zero expression values (“dropouts”), along with substantial cell-to-cell variability. MAST accounts for these characteristics more effectively by modeling both gene detection and expression levels, making it more suitable for identifying biologically meaningful differential expression patterns in single-cell data. In addition, technical confounding factors such as sequencing depth (nCount_RNA) and mitochondrial gene percentage (percent.mt) were included as latent variables during the analysis to reduce technical bias. 

# Run the DEG function created above
```bash
run_mast_per_celltype_DKD_vs_Healthy(
  seu = seurat_obj,
  celltype_col = "majority_voting",   
  condition_col = "Condition"
)
```
After defining the differential expression analysis function, I executed it on the Seurat object to perform cell type–specific differential expression analysis between DKD and Healthy samples. The analysis used the majority_voting metadata column for CellTypist-based cell type annotations and the Condition column to define the disease groups for comparison. The function automatically processed each cell type individually and generated separate differential expression result files for downstream pathway and functional enrichment analysis. 

<img width="895" height="477" alt="image" src="https://github.com/user-attachments/assets/a5d308bc-d3c1-475b-8ab1-67b0579b8501" />    

<img width="776" height="689" alt="image" src="https://github.com/user-attachments/assets/762eba37-215e-465b-832a-dec66eee0226" />

During execution, the function processed each CellTypist-annotated cell type individually and checked whether sufficient numbers of DKD and Healthy cells were available for reliable differential expression analysis. Cell types with too few cells in either comparison group were automatically skipped to avoid statistically unreliable results. Only cell populations that passed the minimum cell count threshold were analyzed using the MAST method, and the resulting differential expression outputs were saved as separate CSV files for downstream pathway enrichment analysis.
After applying the filtering criteria and minimum cell count threshold, a total of eight cell populations contained sufficient DKD and Healthy cells for reliable differential expression analysis. Differentially expressed gene (DEG) results were successfully generated for the following cell types: C-TAL, CNT, DCT, EC-GC, EC-PTC, Podocytes, PT, and VSMC/P. Separate DEG result files were saved for each cell population for downstream functional enrichment and pathway analysis.

# Summerization of DEG Results
```bash
library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(enrichplot)

deg_path <- "D:/Bidya Work/single/GSE183276/output/DEG_MAST_DKD_vs_Healthy"   ##Define folder containing DEG result files generated from MAST analysis

deg_files <- list.files(    ### list.files() scans the folder and retrieves: files beginning with "DEG_" and ending with ".csv"
  deg_path,
  pattern = "^DEG_.*\\.csv$",
  full.names = TRUE
)

deg_files

# --------------------------------------------------
# Fix gene column name ONLY for DEG files
# --------------------------------------------------
##Some exported CSV files may accidentally store gene names in unnamed column ("X" or blank). This loop standardizes the first column name to "genes"

for (file in deg_files) {
  df <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  
  if (colnames(df)[1] %in% c("", "X")) { ## # Check whether first column name is missing or "X"
    colnames(df)[1] <- "genes"
  }
  
  write.csv(df, file, row.names = FALSE)
}
# --------------------------------------------------
# Create empty summary table for storing
# DEG statistics across all cell types
# --------------------------------------------------

deg_summary <- data.frame(
  CellType = character(),     # Cell type name
  Upregulated = integer(),       # Number of significantly upregulated genes
  Downregulated = integer(),     # Number of significantly downregulated genes
  Total_DEGs = integer(),     # Total significant DEGs
  stringsAsFactors = FALSE
)
# --------------------------------------------------
# Loop through each DEG file
# and calculate DEG statistics
# --------------------------------------------------

for (file in deg_files) {      # Read DEG result file
  deg_data <- read.csv(file)
  
  if (!all(c("avg_log2FC", "p_val_adj") %in% colnames(deg_data))) {       # Ensure required DEG columns exist, avg_log2FC = fold change, p_val_adj = adjusted p-value (FDR corrected)
    warning(paste("Skipping (missing columns):", basename(file)))
    next
  }
# Count significantly UPREGULATED genes

  up   <- sum(deg_data$avg_log2FC > 0 & deg_data$p_val_adj < 0.05)  # log2FC > 0  → higher in DKD adjusted p-value < 0.05 → statistically significant
  down <- sum(deg_data$avg_log2FC < 0 & deg_data$p_val_adj < 0.05)    # Count significantly DOWNREGULATED genes. log2FC < 0 → lower in DKD / higher in Healthy
  
  celltype <- gsub("^DEG_(.*?)_DKD_vs_Healthy\\.csv$", "\\1",     # Extract cell type name from file name basename(file))
  
  deg_summary <- rbind(      # Add DEG statistics into summary table
    deg_summary,
    data.frame(
      CellType = celltype,
      Upregulated = up,
      Downregulated = down,
      Total_DEGs = up + down
    )
  )
}
# --------------------------------------------------
# Sort summary table by total DEG count
# Highest DEG-containing cell types appear first
# --------------------------------------------------

deg_summary <- deg_summary %>% arrange(desc(Total_DEGs))

print(deg_summary)

write.csv(
  deg_summary,
  "D:/Bidya Work/single/GSE183276/output/DEG_MAST_DKD_vs_Healthy/DEGsummary_table.csv",
  row.names = FALSE
)
```
<img width="810" height="282" alt="image" src="https://github.com/user-attachments/assets/c7062608-151f-4d6a-a673-afb80cff1bbb" />    

I performed a cell type–resolved differential expression summary analysis to quantify transcriptional changes between DKD and Healthy kidney samples using previously generated MAST-based DEG results. The main objective was to move from gene-level differential expression outputs to a higher-level view of how strongly each renal cell type is affected in disease.    
All DEG result files were systematically aggregated and standardized to ensure consistent formatting across cell types before analysis. Only results with validated MAST output structure were included to maintain reliability. For each cell type, I focused on significantly differentially expressed genes based on a consistent threshold (adjusted p-value < 0.05), while also distinguishing the direction of change to capture both upregulated and downregulated transcriptional shifts in DKD. The key decision in this step was to prioritize a uniform statistical cutoff and apply it across all cell types, enabling fair comparison of DEG burden rather than gene-specific interpretation. This allowed identification of cell populations with the strongest transcriptional perturbations in DKD. Cell type identities were programmatically extracted from filenames to ensure traceability and reproducibility without manual annotation.    
Finally, I generated a consolidated summary table ranking cell types based on total differential expression load. This provides a comparative landscape of disease impact across renal compartments and serves as a foundation for downstream visualization and pathway-level interpretation.  

From the results, epithelial and tubular compartments show the strongest transcriptional perturbation, with EC_PTC exhibiting the highest number of DEGs, followed by C_TAL and CNT, indicating these segments are major contributors to DKD-associated molecular changes. In contrast, cell types like DCT and Podocytes show minimal differential expression under the same thresholds, suggesting either lower transcriptional responsiveness or more stable gene regulation in the dataset context.    
Importantly, separating upregulated and downregulated genes helps capture not just disease intensity but also directionality of transcriptional remodeling, highlighting whether DKD is associated with activation or suppression of gene programs in each cell type.    

# Cell type–resolved Reactome pathway enrichment analysis using ranked gene lists (GSEA framework)
```bash
all_pathways_df <- data.frame()  # Create empty dataframe to store Reactome pathway enrichment results from ALL cell types together

for (file_path in deg_files) {

  file_name <- basename(file_path)     # DEG_Podocyte_DKD_vs_Healthy.csv
  message("\n==============================")     # Print separator for cleaner console output
  message("▶ Processing file: ", file_name)    # Print currently processed DEG file

  ### Extract cell type name from file name,  # Removes: "DEG_" and "_DKD_vs_Healthy.csv", leaving only cell type name
  celltype <- sub("^DEG_(.*?)_DKD_vs_Healthy\\.csv$", "\\1", file_name)

  if (celltype == file_name) {    # If extraction failed, celltype remains identical to original file name
    warning("⚠ Celltype extraction FAILED for file: ", file_name)
    next
  }

  celltype_clean <- gsub("[^a-zA-Z0-9]", "_", celltype)    # Clean special characters from cell type name
  condition <- "DKD_vs_Healthy"   # Store biological comparison label

  message("  ✔ Cell type identified: ", celltype)   # Print extracted cell type

  ### ---- Read DEG file ----
  res <- tryCatch(     # tryCatch prevents entire pipeline from crashing if one file is corrupted or unreadable
    read.csv(file_path, stringsAsFactors = FALSE),    # Read DEG CSV file
    error = function(e) {    # Error handling function
      message(" Failed to read file: ", e$message)  # Print error message
      return(NULL)   # Return NULL instead of crashing
    }
  )
  if (is.null(res)) next    # Skip file if reading failed

  message("  ✔ DEG rows read: ", nrow(res))      # Print number of DEG rows loaded

    # Verify required DEG columns exist
  # gene = gene symbols
  # avg_log2FC = fold-change ranking metric
  if (!all(c("gene", "avg_log2FC") %in% colnames(res))) {
    warning("  ⚠ Missing required columns in: ", file_name)
    next
  }

 # Convert Gene SYMBOLs → ENTREZ IDs
    # Reactome enrichment requires ENTREZ IDs instead of gene symbols
  ncbi_map <- suppressMessages(
    clusterProfiler::bitr(
      res$gene,                     # Input gene symbols
      fromType = "SYMBOL",          # Current gene ID type
      toType = "ENTREZID",          # Desired gene ID type
      OrgDb = org.Hs.eg.db         # Human annotation database
    )
  )

  message("  ✔ Genes mapped to ENTREZ: ", nrow(ncbi_map))    # Print number of successfully mapped genes
  
   # --------------------------------------------------
  # Merge DEG table with ENTREZ mapping
  # -------------------------------------------------- 

  res_mapped <- res %>%
    left_join(ncbi_map, by = c("gene" = "SYMBOL")) %>%     # Match DEG genes with ENTREZ IDs
    filter(!is.na(ENTREZID)) %>%          # Remove genes without ENTREZ mapping
    distinct(ENTREZID, .keep_all = TRUE)

  message(" Genes after filtering: ", nrow(res_mapped))   # Print remaining usable genes

  # --------------------------------------------------
  # Create ranked gene list for GSEA
  # --------------------------------------------------
  # GSEA requires: named numeric vector
    # Values = fold changes
    # Names = ENTREZ IDs


  gene_list <- res_mapped$avg_log2FC
  names(gene_list) <- res_mapped$ENTREZID      # Assign ENTREZ IDs as vector names
  
  # Sort genes from highest positive logFC to most negative logFC
  # Highly positive genes: upregulated in DKD
  # Highly negative genes: downregulated in DKD
  gene_list <- sort(gene_list, decreasing = TRUE)

  message("  ✔ Ranked gene list length: ", length(gene_list))     # Print gene list size
  
 # GSEA becomes unreliable with very small gene sets
  if (length(gene_list) < 20) {
    message("  ⏭ Skipping GSEA: too few genes")
    next
  }

 # --------------------------------------------------
  # Run Reactome pathway GSEA
  # --------------------------------------------------
   # GSEA asks: "Are genes from specific pathways enriched at the top or bottom of the ranked DEG list?"
#Positive NES: pathway activated/upregulated in DKD
# Negative NES: pathway suppressed/downregulated in DKD

  message(" Running Reactome GSEA...")

  gsea_res <- tryCatch(
    gsePathway(
      geneList = gene_list,     # Ranked gene list
      organism = "human",       # Human Reactome pathways
      eps = 0,                  # More precise p-value estimation
      verbose = FALSE             # Suppress verbose output
    ),
    error = function(e) {     # Handle GSEA errors safely
      message("  GSEA FAILED: ", e$message)
      return(NULL)
    }
  )

  if (is.null(gsea_res) || nrow(gsea_res@result) == 0) {
    message("  ⏭ No significant pathways")
    next
  }

##### Convert ENTREZ IDs back into readable gene names
  gsea_res <- setReadable(
    gsea_res,    # GSEA object
    OrgDb = org.Hs.eg.db,  # Human annotation database
    keyType = "ENTREZID"     # Original ID type
  )
  
 ## Extract clean pathway result table
  pathways <- gsea_res@result %>%   # Remove incomplete/invalid pathways
    filter(
      !is.na(Description),
      !is.na(NES),
      !is.na(p.adjust),
      !is.na(setSize)
    )

  message(" Valid pathways: ", nrow(pathways)) # Print number of valid pathways

  if (nrow(pathways) == 0) next  # Skip empty pathway tables

  ##  # Add metadata columns
  pathways$CellType  <- celltype
  pathways$Condition <- condition

  all_pathways_df <- rbind(all_pathways_df, pathways)    # Combine current pathway results with master pathway dataframe

  ## ---- Save per-cell CSV ----
  out_csv <- file.path(    # Save pathway enrichment results for current cell type
    deg_path,
    paste0(celltype_clean, "_Reactome_GSEA_", condition, ".csv")
  )
  write.csv(pathways, out_csv, row.names = FALSE)

  message("  ✔ GSEA CSV written for ", celltype)
}
```
<img width="1048" height="513" alt="image" src="https://github.com/user-attachments/assets/8ba84191-851c-42e1-a6a2-6fe43f145a50" />

I performed a cell type–specific pathway enrichment analysis using Reactome Gene Set Enrichment Analysis (GSEA) to interpret the biological processes underlying differential gene expression in DKD vs Healthy kidney samples. Instead of focusing only on individual differentially expressed genes, I extended the analysis to pathway-level interpretation to understand coordinated functional changes across renal cell types. For each cell type, DEG results were first read and filtered to ensure the presence of required gene identifiers and fold-change values. Gene symbols were converted to Entrez IDs to ensure compatibility with Reactome pathway annotations, as pathway databases require standardized identifiers for enrichment mapping. Only successfully mapped genes were retained to maintain annotation accuracy.  
A ranked gene list was then constructed using log fold-change values, where genes were ordered from most upregulated to most downregulated in DKD. This ranking strategy was chosen because GSEA does not rely on hard significance cutoffs but instead evaluates whether genes from a pathway are systematically enriched at the top or bottom of the ranked list. This allows detection of subtle but coordinated pathway-level changes that may not be captured by threshold-based DEG filtering.   
Reactome GSEA was performed separately for each cell type to preserve cell type–specific biological signals. Robust error handling was included to ensure that individual file failures or insufficient gene sets did not interrupt the full pipeline. Finally, all pathway enrichment results were consolidated into a single dataset while also exporting individual cell type–specific outputs. This structure enables both global comparison across renal cell types and focused interpretation of pathway dysregulation within specific compartments.    

<img width="422" height="126" alt="image" src="https://github.com/user-attachments/assets/31b9cc08-0c67-4edd-a8d7-163cc3c4d007" /> 
<img width="478" height="120" alt="image" src="https://github.com/user-attachments/assets/9043d535-318a-4776-9b9d-bf67bc5d7439" />  
<img width="508" height="120" alt="image" src="https://github.com/user-attachments/assets/acd6a3bd-ce29-4b64-aba8-8ef64dc553b6" />
<img width="453" height="121" alt="image" src="https://github.com/user-attachments/assets/a1d21a2a-61c2-4730-8b32-9dddc5c539ef" /> 
<img width="357" height="73" alt="image" src="https://github.com/user-attachments/assets/1ebdb5d9-9604-4fa1-8c3c-e747e593addf" />

I performed Reactome-based GSEA across all cell type–specific DEG files to interpret functional pathway-level changes in DKD vs Healthy kidney samples. Although the analysis was applied uniformly to all cell types, only a subset of them produced significant pathway enrichment results after filtering for mapping success, gene list quality, and adjusted p-value thresholds. This was mainly due to differences in DEG burden across cell types and the loss of genes during Entrez ID conversion, which reduced the effective input size for some populations. As a result, only cell types with stronger and more coordinated transcriptional signals contributed meaningful pathway enrichment outputs. This indicates that pathway-level biological changes in DKD are not uniform across all renal compartments but are instead driven by specific cell populations showing higher transcriptional perturbation.    

<img width="1897" height="208" alt="image" src="https://github.com/user-attachments/assets/199b8ed9-4171-42c3-bab0-f426e8f5ea20" />    

# Reactome GSEA Pathway Visualization
```bash
gsea_png_files <- c()  # Create empty vector to store paths of generated GSEA pathway plot PNG files

# --------------------------------------------------
# Loop through each unique cell type present in
# the combined pathway enrichment dataframe
# --------------------------------------------------

for (ct in unique(all_pathways_df$CellType)) {

  df_ct <- all_pathways_df %>%     # Extract pathway enrichment results for current cell type only
    filter(CellType == ct)

  if (nrow(df_ct) == 0) {             # Skip plotting if no pathways exist
    message(" No pathways to plot")
    next
  }

  ## ---- Select top pathways ----
    # NES = Normalized Enrichment Score
  # Higher positive NES: pathways activated in DKD
  # Higher negative NES: pathways suppressed in DKD
 # abs(NES) selects strongest biological signals regardless of direction

  top_pathways <- df_ct %>%
    arrange(desc(abs(NES))) %>%     # Sort pathways by strongest enrichment magnitude
    slice_head(n = 10)     # Keep top 10 pathways

 if (nrow(top_pathways) == 0) next

  # --------------------------------------------------
  # Reorder pathway names according to NES values
  # Improves plot readability
  # --------------------------------------------------
  top_pathways$Description <- factor(   # Pathways with low NES appear at bottom, Pathways with high NES appear at top
    top_pathways$Description,
    levels = top_pathways$Description[order(top_pathways$NES)]
  )

  # --------------------------------------------------
  # Create pathway enrichment dot plot
  # --------------------------------------------------
 # x-axis: NES = enrichment direction + strength
 # y-axis: pathway names
 # point color: adjusted p-value significance
 # point size: pathway gene-set size
  p <- ggplot(
    top_pathways,
    aes(x = NES, y = Description, color = p.adjust, size = setSize)
  ) +
    geom_point() +
    theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7, face = "bold"),
    axis.text.x = element_text(size = 7, face = "bold")
  )
  ## ---- Save PNG ----
  ct_clean <- gsub("[^a-zA-Z0-9]", "_", ct)     # Clean cell type name for safe file naming
 
  png_file <- file.path(
    deg_path,
    paste0(ct_clean, ".png")
  )

  png(png_file, width = 9, height = 6, units = "in", res = 300)
  print(p)
  dev.off()

  gsea_png_files <- c(gsea_png_files, png_file)
}
write.csv(
  all_pathways_df,
  file.path(deg_path, "All_CellTypes_Reactome_GSEA_results.csv"),
  row.names = FALSE
)
```
I performed pathway-level enrichment visualization separately for each identified cell type using the combined Reactome GSEA output dataframe. First, I iterated through all unique cell types present in the enrichment results and extracted cell type–specific pathway enrichment profiles. To avoid generating empty or biologically irrelevant plots, I implemented conditional filtering steps that skipped cell types with no enriched pathways. For each cell type, pathways were ranked based on the absolute value of the Normalized Enrichment Score (NES), allowing prioritization of the strongest biological signals irrespective of enrichment direction. Positive NES values were interpreted as pathways activated in DKD, whereas negative NES values represented suppressed pathways. I selected the top 10 pathways with the highest absolute NES values to focus on the most biologically significant enrichment signatures while maintaining plot interpretability.
To improve visualization clarity, pathway descriptions were reordered according to NES values so that negatively enriched pathways appeared at the bottom and positively enriched pathways appeared at the top of the plot. I then generated GSEA dot plots using ggplot2, where the x-axis represented NES values, the y-axis represented Reactome pathway names, point color encoded adjusted p-values (FDR-corrected significance), and point size reflected pathway gene-set size. 

<img width="1129" height="915" alt="image" src="https://github.com/user-attachments/assets/26421742-90c5-4d41-b74d-146c2c39c38f" />    

<img width="1106" height="871" alt="image" src="https://github.com/user-attachments/assets/970365e4-ad30-4916-8039-c7d7ca6f7dd7" />    

<img width="1099" height="878" alt="image" src="https://github.com/user-attachments/assets/3bf225d9-9e63-435a-b5e8-1e30ebb59676" />    

<img width="1097" height="873" alt="image" src="https://github.com/user-attachments/assets/dedd9443-eadc-4bcc-ba0f-988c28d7330b" />     

<img width="1094" height="877" alt="image" src="https://github.com/user-attachments/assets/6f9f3efd-008b-4761-b360-b15933452dc8" />     

<img width="1886" height="450" alt="image" src="https://github.com/user-attachments/assets/b76ce705-3881-4ba6-95db-fa85f45b5584" />     

Reactome pathway enrichment analysis using ranked gene lists revealed pronounced cell type–specific functional remodeling in diabetic kidney disease (DKD) compared with healthy controls. Across tubular epithelial populations, including the cortical thick ascending limb (C_TAL), connecting tubule (CNT), distal convoluted tubule (DCT), and proximal tubule (PT), distinct patterns of transcriptional dysregulation were observed. C_TAL and CNT populations demonstrated marked suppression of neurotrophic signaling programs, including NGF-stimulated transcription, Signaling by NTRK1 (TRKA), and broader NTRK signaling pathways, together with reduced activation of nuclear transcription factor–associated processes. These findings suggest widespread attenuation of growth factor–responsive transcriptional networks and stress-adaptive signaling in tubular epithelial compartments during DKD progression.    
In contrast, CNT cells additionally exhibited significant enrichment of extracellular matrix (ECM) organization and ECM degradation pathways, indicating concurrent activation of tissue remodeling and fibrotic processes. Similarly, DCT cells displayed selective enrichment of ECM-associated pathways, suggesting more subtle but directionally consistent matrix remodeling within distal nephron segments. In proximal tubule (PT) cells, enrichment of Golgi-associated vesicle biogenesis pathways indicated alterations in intracellular trafficking, vesicular transport, and protein-processing machinery in response to diabetic stress.    
Among all analyzed populations, vascular smooth muscle/perivascular cells (VSMC_P) exhibited some of the most pronounced transcriptional alterations. These changes were characterized by suppression of pathways involved in elastic fibre formation, molecules associated with elastic fibres, and insulin-like growth factor (IGF) transport and uptake, alongside dysregulation of glycosaminoglycan metabolism and ECM-associated processes. Collectively, these findings indicate substantial disruption of vascular structural integrity, extracellular matrix homeostasis, and growth factor signaling within vascular compartments.    
Overall, pathway-level analysis demonstrated that DKD induces highly compartmentalized molecular reprogramming across kidney cell populations, characterized by suppression of neurotrophic and transcription-associated signaling in tubular epithelial cells alongside extensive extracellular matrix remodeling and structural degeneration in vascular and stromal compartments. These results highlight distinct yet coordinated injury-associated functional programs contributing to DKD pathogenesis.    

Overall, these findings show that DKD causes strong cell type–specific molecular changes across the kidney, marked by reduced signaling activity in tubular cells and increased extracellular matrix remodeling and structural damage in vascular compartments. These findings improve our understanding of how different kidney cell populations respond to diabetic injury at the pathway level.     

This scRNA-seq analysis identified distinct cell-type-specific transcriptional alterations associated with kidney disease progression. Differential expression and pathway enrichment analyses revealed that proximal tubule (PT), cortical thick ascending limb (C_TAL), connecting tubule (CNT), distal convoluted tubule (DCT), and vascular smooth muscle cell (VSMC) populations exhibited significant molecular changes under disease conditions. Enrichment results highlighted dysregulation of pathways related to cellular stress responses, ion transport, inflammation, metabolic dysfunction, and vascular remodeling, suggesting that multiple nephron segments and vascular compartments contribute to disease pathology. In particular, PT cells demonstrated strong injury-associated signatures, while TAL-, CNT-, and DCT-associated populations reflected alterations in electrolyte transport and tubular homeostasis. Additionally, VSMC-related pathway changes indicated potential vascular involvement during kidney injury progression.    

Overall, this study demonstrates the utility of single-cell RNA sequencing in uncovering cell-type-specific molecular heterogeneity that may not be detectable using bulk transcriptomic approaches, providing insights into the complex cellular mechanisms underlying kidney disease.    














