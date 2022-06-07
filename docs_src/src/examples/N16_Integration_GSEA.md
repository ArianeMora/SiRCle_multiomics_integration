# GSEA on the integrated VAE value

## ORA setup
Run ORA on each of the regulatory clusters. The files are saved in the current working directory.
```
library(org.Hs.eg.db)
library(clusterProfiler)
library(svglite)
library(sircle)
source("helper_functions.R")

GSEA_file_dir <- "../data/S050_CCRCC_Clark_Cell2019/supps/GSEA/"
 
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


## Look at the changes between patient groups using the VAE

```
input_data_dir <- '../data/sircle/F6_integrative_comparisons/'
figure_dir <- '../figures/'

# Use the RCM df as the reference (this is how we get the background genes for the clusters for ORA)
refDf <- read.csv(paste0('../data/sircle/F3_regulatory_clustering/', "RCM_all_patients_ccRCC_P0.5-R1.0-M0.1-GENES.csv"))

reg_grps <- c("MDS", "MDS_TMDE", "MDE", "MDE_TMDS", "TMDE", "TMDS", "TPDE", "TPDE_TMDS",  "TPDS", "TPDS_TMDE")
sircleFileName <- paste0(input_data_dir, 'mean_Integrated_comparison_Stage IV-Stage I.csv')
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
vis_df <- read.csv(sircleFileName)
for (r in reg_grps) {

  grpGenes <- subset(vis_df, vis_df[["Regulation_Grouping_2"]] == r)
  # Run GSEA on the consistently changed genes between late and early stage using the t-stat as the rank
  GSEA_results <- GSEASetupSircle(grpGenes, paste0("G4_Stage IV-Stage I ", r), gene_name="external_gene_name", stat='Integrated.diff..Stage.IV.Stage.I.',   GSEA_file_dir=GSEA_file_dir, output_dir=figure_dir, pathway_list=pathway_list)
}
  
sircleFileName <- paste0(input_data_dir, 'mean_Integrated_comparison_old-young.csv')
vis_df <- read.csv(sircleFileName)
for (r in reg_grps) {

  grpGenes <- subset(vis_df, vis_df[["Regulation_Grouping_2"]] == r)
  # Run GSEA on the consistently changed genes between late and early stage using the t-stat as the rank
  GSEA_results <- GSEASetupSircle(grpGenes, paste0("G4_old-young ", r), gene_name="external_gene_name", stat='Integrated.diff..old.young.',   GSEA_file_dir=GSEA_file_dir, output_dir=figure_dir, pathway_list=pathway_list)
}

## Lastly look at mutations
sircleFileName <- paste0(input_data_dir, 'mean_Integrated_comparison_PBRM1-BAP1.csv')
vis_df <- read.csv(sircleFileName)
for (r in reg_grps) {

  grpGenes <- subset(vis_df, vis_df[["Regulation_Grouping_2"]] == r)
  # Run GSEA on the consistently changed genes between late and early stage using the t-stat as the rank
  GSEA_results <- GSEASetupSircle(grpGenes, paste0("G4_PBRM1-BAP1 ", r), gene_name="external_gene_name", stat='Integrated.diff..PBRM1.BAP1.',   GSEA_file_dir=GSEA_file_dir, output_dir=figure_dir, pathway_list=pathway_list)
}

```




