#Protein and RNAseq DE analyses and CpG differential methylation analysis

# Tumour vs normal comparison
This notebook compares tumour vs normal samples in the three data types:
1) Protein data  
2) RNAseq data   
3) DNA Methylation data  

## DE DA CpG analysis
```
library("limma")
library("dplyr")
library(matrixStats)

source("helper_functions.R")

##Chdir to source dir if you aren't running from the relative location.
GSEA_file_dir <- '../data/raw_downloads/supps/GSEA/'
figure_dir <- '../figures/'
input_data_dir <- '../data/sircle/F1_DE_input_TvN/'
output_data_dir <- '../data/sircle/F2_DE_output_TvN/'

# Setup the pathways for using GSEA
KEGG <- read.csv(file.path(GSEA_file_dir,  "c2.cp.kegg.v6.2.symbols.csv"))
Reactome <- read.csv(file.path(GSEA_file_dir,  "c2.cp.reactome.v6.2.symbols.csv"))
Biocarta <- read.csv(file.path(GSEA_file_dir, "c2.cp.biocarta.v6.2.symbols.csv"))
Hallmarks <- read.csv(file.path(GSEA_file_dir,  "h.all.v6.2.symbols.csv"))
GO_BP <- read.csv(file.path(GSEA_file_dir, "c5.go.bp.v7.2.symbols.csv"))
GO_CC <- read.csv(file.path(GSEA_file_dir, "c5.go.cc.v7.2.symbols.csv"))
GO_MF <- read.csv(file.path(GSEA_file_dir, "c5.go.mf.v7.2.symbols.csv"))
Metabolic <- read.csv(file.path(GSEA_file_dir, "41467_2016_BFncomms13041_MOESM340_ESM.csv"))
Correction <- read.csv(file.path(GSEA_file_dir, "41467_2016_BFncomms13041_MOESM340_ESM.csv"))

#Run the GSEA analysis
##1."KEGG", "Reactome", "Biocarta", "Hallmarks"
pathways <- rbind(KEGG, Reactome, Biocarta, Hallmarks)
pathway_list <- list()
for (pathway in unique(pathways$term)) {
  pathway_list[[pathway]] <- as.character(pathways[pathways$term == pathway, 1])
}
```


## Perform differential analysis on methylation, RNAseq and Protein data
```
# Do RNAseq analysis on different patient subsets, so we'll need to read in the sample DF
# and also the RNAseq data.
# Only need to run the TvN comparison because we will be doing the late and early stage comparisons after we do the 
# tests.
runDEAll('all_patients_ccRCC', input_data_dir, output_data_dir, GSEA_file_dir, paired=FALSE, output_figure_dir=figure_dir, pathway_list=pathway_list)
```

