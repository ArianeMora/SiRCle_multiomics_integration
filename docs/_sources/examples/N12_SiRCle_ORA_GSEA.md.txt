# ORA for SiRCle clusters

## ORA setup
Run ORA on each of the regulatory clusters. The files are saved in the current working directory.
```
library(org.Hs.eg.db)
library(clusterProfiler)
library(svglite)
library(sircle)
# For running the immune cell shit install.packages("getopt")
source("helper_functions.R")

GSEA_file_dir <- '../data/raw_downloads/supps/GSEA/'
figure_dir <- '../figures/'
input_data_dir <- '../data/sircle/F1_DE_input_TvN/'
output_data_dir <- '../data/sircle/F2_DE_output_TvN/'
 
#Import the pathways:
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


### ORA on the Tumour vs normal 

```
input_data_dir <- '../data/S050_CCRCC_Clark_Cell2019/sircle/all_patients/all_cpgs/'
figure_dir <- '../figures/sircle/all_patients/all_cpgs/'

regLabel <- 'Regulation_Grouping_2'
test_title <- 'RCM_P0.5-R1.0-M0.1'

# Read in the regulatory clustering output
sircleFileName <- paste0(input_data_dir, "RCM_all_patients_ccRCC_P0.5-R1.0-M0.1-GENES.csv")

# Plot the size of each regulatory cluster as a circle
sirclePlot(sircleFileName, output_dir=figure_dir, regLabels=regLabel, title=paste(test_title), fileType="pdf")

# Run ORA on each of the clusters
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, regLabel, emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, output_dir = figure_dir)

```