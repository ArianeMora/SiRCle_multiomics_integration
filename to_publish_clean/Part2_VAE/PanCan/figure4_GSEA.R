
library(EnhancedVolcano)
library(org.Hs.eg.db)#For Human load:
library(clusterProfiler)#To run ORA
library(enrichplot)#For emapplot, dotplot,...
library(ggplot2)#For saving and legends
library(tidyverse) # used for data manipulation
#Establish a couple of functions needed to perform the GSEA (made by Aurelien)
library(fgsea)
library(GSEABase)



GSEASetupSircle <- function (GSEA_data, title, gene_name="external_gene_name", stat="stat", GSEA_file_dir=NULL, output_dir=NULL, 
                             pathway_list=NULL) {
  # Set output dir to be the current directory if it is null otherwise set 
  if (is_null(output_dir)) {
    output_dir <- getwd()
  }
  if (is_null(GSEA_file_dir)) {
    GSEA_file_dir <- getwd()  # Set it by default to be the current working directory
  }
  
  #Prepare data matrix for GSEA:
  Column_Gene <- as.character(GSEA_data[[gene_name]])
  Column_tval <- as.numeric(GSEA_data[[stat]])
  MyData_Extracted <- data.frame(cbind(Column_Gene, Column_tval), stringsAsFactors = F)
  MyData_Extracted$Column_tval <- as.numeric(MyData_Extracted$Column_tval)
  
  t_val <- as.numeric(MyData_Extracted$Column_tval)#Make the data into a vector
  names(t_val) <- MyData_Extracted$Column_Gene
  
  if (is_null(pathway_list)) {
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
    pathways <- rbind(Hallmarks, KEGG, Reactome, Biocarta, Hallmarks)
    pathway_list <- list()
    for (pathway in unique(pathways$term)) {
      pathway_list[[pathway]] <- as.character(pathways[pathways$term == pathway, 1])
    }
  }
  
  filename <- file.path(output_dir, gsub(" ", "_", title))
  
  gsea_result1 <- fgsea(pathways = pathway_list, stats = t_val, nperm = 10000, minSize=4)
  write_csv(gsea_result1, paste0(filename, "_GSEA_Pathways.csv"))
  
  ##2."GO_BP", "GO_CC", "GO_MF"
  #pathways <- rbind(GO_BP, GO_CC, GO_MF)
  #pathway_list <- list()
  #for(pathway in unique(pathways$term))
  #{pathway_list[[pathway]] <- as.character(pathways[pathways$term == pathway, 1])}
  
  #gsea_result2 <- fgsea(pathways = pathway_list, stats = t_val, nperm = 10000, minSize=4)
  #write_csv(gsea_result2, paste0(filename,  "_GSEA_GO-terms.csv"))
  
  #3 Metabolic enzymes
  pathways <- Metabolic
  pathway_list <- list()
  for(pathway in unique(pathways$term))
  {pathway_list[[pathway]] <- as.character(pathways[pathways$term == pathway, 1])}
  
  gsea_result3 <- fgsea(pathways = pathway_list, stats = t_val, nperm = 10000, minSize=4)
  write_csv(gsea_result3, paste0(filename,  "_GSEA_MetabolicPathways.csv"))
  
  #### Promiscuity Corrected metabolic enzymes
  # Promiscuity correction:
  #Prepare data matrix for GSEA:
  #Column_Gene <- as.character(GSEA_data[[gene_name]])
  #Column_tval <- as.numeric(GSEA_data[[stat]])
  #MyData_Extracted <- data.frame(cbind(Column_Gene, Column_tval), stringsAsFactors = F)
  #MyData_Extracted$Column_tval <- as.numeric(MyData_Extracted$Column_tval)
  
  #t_val <- as.numeric(MyData_Extracted$Column_tval)#Make the data into a vector
  #names(t_val) <- MyData_Extracted$Column_Gene
  
  #Correction <- read.csv(file.path(GSEA_file_dir, "41467_2016_BFncomms13041_MOESM341_ESM.csv"))
  #MyData_Extracted_Corr <- merge(x=MyData_Extracted, by.x="Column_Gene", y=Correction, by.y="gene")
  #MyData_Extracted_Corr$Column_tval_Corr <- (MyData_Extracted_Corr$Column_tval)/(MyData_Extracted_Corr$pathway)
  
  # matrix for GSEA:
  #t_val <- as.numeric(MyData_Extracted_Corr$Column_tval_Corr)#Make the data into a vector
  #names(t_val) <- MyData_Extracted_Corr$Column_Gene
  
  #3 Metabolic enzymes
  #pathways <- Metabolic  # Ariane update
  #pathway_list <- list()
  #for(pathway in unique(pathways$term))
  #{pathway_list[[pathway]] <- as.character(pathways[pathways$term == pathway, 1])}
  
  #gsea_result3_Corrected <- fgsea(pathways = pathway_list, stats = t_val, nperm = 10000, minSize=4)
  #write_csv(gsea_result3_Corrected, paste0(filename, "_GSEA_MetabolicPathways_PromiscuityCorrected.csv"))
  
  # Return each of the GSEA results
  #return(list("pathways" = gsea_result1, "GO" = gsea_result2, "MetabolicPathways" = gsea_result3, "MetabolicPathwaysCorrected" =  gsea_result3_Corrected))
}


#  ----------------------------------------------------------------------------------------------
#                  File and naming setup
#  ----------------------------------------------------------------------------------------------

GSEA_file_dir <- 'Required_Refs/GSEA/' # Will need to change the path if you're on windows
figure_dir <- 'Output_Figures'
input_data_dir <- 'Output_Data'
output_data_dir <- 'Output_Data/ORA/'
cancer <- 'PanCan'

#Import the pathways:
KEGG <- read.csv(file.path(GSEA_file_dir,  "c2.cp.kegg.v6.2.symbols.csv"))
Reactome <- read.csv(file.path(GSEA_file_dir,  "c2.cp.reactome.v6.2.symbols.csv"))
Biocarta <- read.csv(file.path(GSEA_file_dir, "c2.cp.biocarta.v6.2.symbols.csv"))
Hallmarks <- read.csv(file.path(GSEA_file_dir,  "h.all.v6.2.symbols.csv"))
GO_BP <- read.csv(file.path(GSEA_file_dir, "c5.go.bp.v7.2.symbols.csv"))
GO_CC <- read.csv(file.path(GSEA_file_dir, "c5.go.cc.v7.2.symbols.csv"))
GO_MF <- read.csv(file.path(GSEA_file_dir, "c5.go.mf.v7.2.symbols.csv"))
Metabolic <- read.csv(file.path(GSEA_file_dir, "41467_2016_BFncomms13041_MOESM340_ESM.csv"))
#Correction <- read.csv(file.path(GSEA_file_dir, "41467_2016_BFncomms13041_MOESM340_ESM.csv"))

#Run the GSEA analysis
##1."KEGG", "Reactome", "Biocarta", "Hallmarks"
pathways <- rbind(KEGG)#, Reactome, Biocarta, Hallmarks)
pathway_list <- list()
for (pathway in unique(pathways$term)) {
  pathway_list[[pathway]] <- as.character(pathways[pathways$term == pathway, 1])
}

cancer <- 'PanCan'
refDf <- read.csv(paste0('Input_RCM/RCM_', cancer, '.csv'))

reg_grps <- c("MDS", "MDS_TMDE", "MDE", "MDE_TMDS", "TMDE", "TMDS", "TPDE", "TPDE_TMDS",  "TPDS", "TPDS_TMDE")


sircleFileName <- file.path(input_data_dir, paste0('stats_Late-Early_', cancer, '.csv'))
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
vis_df <- read.csv(sircleFileName)
for (r in reg_grps) {
  
  grpGenes <- subset(vis_df, vis_df[["RG2_Changes_filtered"]] == r)
  # Run GSEA on the consistently changed genes between late and early stage using the t-stat as the rank
  GSEA_results <- GSEASetupSircle(grpGenes, paste0("Late vs Early ", r), gene_name="id", stat='Integrated.diff..Late.Early.',   GSEA_file_dir=GSEA_file_dir, output_dir=figure_dir, pathway_list=pathway_list)
}
sircleFileName <- file.path(input_data_dir, paste0('stats_Late-Early_', cancer, '.csv'))
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
vis_df <- read.csv(sircleFileName)
GSEA_results <- GSEASetupSircle(vis_df, paste0("Late vs Early ", 'All'), gene_name="id", stat='Integrated.diff..Late.Early.',   GSEA_file_dir=GSEA_file_dir, output_dir=figure_dir, pathway_list=pathway_list)



# ------------- Now look at S4 vs S1
sircleFileName <- file.path(input_data_dir, paste0('stats_Stage IV-Stage I_', cancer, '.csv'))
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
vis_df <- read.csv(sircleFileName)
for (r in reg_grps) {
  
  grpGenes <- subset(vis_df, vis_df[["RG2_Changes_filtered"]] == r)
  # Run GSEA on the consistently changed genes between late and early stage using the t-stat as the rank
  
}

# Do it for all pathways simultaeneously
GSEA_results <- GSEASetupSircle(vis_df, paste0("Stage IV vs Stage I ", "All"), gene_name="id", stat='Integrated.diff..Stage.IV.Stage.I.',   GSEA_file_dir=GSEA_file_dir, output_dir=figure_dir, pathway_list=pathway_list)
