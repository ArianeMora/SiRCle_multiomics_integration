library(EnhancedVolcano)
library(org.Hs.eg.db)#For Human load:
library(clusterProfiler)#To run ORA
library(enrichplot)#For emapplot, dotplot,...
library(ggplot2)#For saving and legends
library(tidyverse) # used for data manipulation
#Establish a couple of functions needed to perform the GSEA (made by Aurelien)
library(fgsea)
library(GSEABase)

#  ----------------------------------------------------------------------------------------------
#                  File and naming setup
#  ----------------------------------------------------------------------------------------------

GSEA_file_dir <- 'Required_Refs/GSEA/' # Will need to change the path if you're on windows
figure_dir <- 'Output_Figures'
input_data_dir <- 'Output_Data'
output_data_dir <- 'Output_Data/ORA/'
cancer <- 'ClearCellRenalCellCarcinoma'

sircleFileName <- file.path(input_data_dir, paste0('sircle_PorMandR_', cancer, '.csv'))

# Run the ORA

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

sircleORAHuman <- function(filename, entrezId, title, regLabels="RegulatoryLabels", emptyRegLabel="", fileType="pdf",
                           minGSSize=10, qvalueCutoff=0.2, pvalueCutoff=0.05, showCatagory=30, outputFolder='') {
  ## ------------ Setup and installs ----------- ##
  packages <- c("org.Hs.eg.db", "clusterProfiler", "svglite", "enrichplot")
  install.packages(setdiff(packages, rownames(installed.packages())))
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(svglite)
  library(enrichplot)
  ## ------------ Run ----------- ##
  # open the data
  df <- read.csv(filename)
  allGenes <- as.character(df[[entrezId]]) #
  clusterGenes <- subset(df, ! df[[regLabels]] == emptyRegLabel)
  clusterGenes <- subset(clusterGenes, ! clusterGenes[[regLabels]] == "Not-Background")

  grps_labels <- unlist(unique(clusterGenes[regLabels]))
  for(g in grps_labels) {
    grpGenes <- subset(df, df[[regLabels]] == g)
    print(g)
    print(dim(grpGenes))
    clusterGo <- enrichGO(gene = as.character(grpGenes[[entrezId]]),
                          universe = allGenes,
                          keyType = "ENTREZID",
                          OrgDb = org.Hs.eg.db,
                          ont = "ALL",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 1.0,
                          pvalueCutoff = 1.0,
                          minGSSize = minGSSize,
                          readable = TRUE)
    # We have a cutoff of all of them, and then only visualise the ones that the user wants...

    clusterGoSummary <- data.frame(clusterGo)
    write.csv(clusterGoSummary, paste(outputFolder, 'ClusterGoSummary_', g, title, '.csv', sep=""))#Export the ORA results as .csv

    if (!(dim(clusterGoSummary)[1] == 0)) {#exclude df's that have no observations
      Dotplot <- dotplot(clusterGo, showCategory=showCatagory) +
        ggtitle(paste("Dotplot ", g, sep=""))
      ggsave(file=paste(outputFolder, "SiRCle-ORA_Dotplot_Human_", g, title, ".", fileType, sep=""), plot=Dotplot, width=10, height=8)
      x2 <- pairwise_termsim(clusterGo)

      Emapplot <- emapplot(x2, pie_scale=1.5, layout = "nicely")+
        ggtitle(paste("Emapplot ", g, sep=""))
      ggsave(file=paste(outputFolder, "SiRCle-ORA_Emapplot_Human_", g, title, ".", fileType, sep="" ), plot=Emapplot, width=10, height=8)

      Heatplot <- heatplot(clusterGo,showCategory=showCatagory) +
        theme(axis.text.x =element_text(size=5), axis.text.y =element_text(size=8,face="bold"), axis.title=element_text(size=12,face="bold"))+
        ggtitle(paste("Heatplot ", g, sep=""))
      ggsave(file=paste(outputFolder, "SiRCle-ORA_Heatplot_Human_", g,title,  ".", fileType, sep="" ), plot=Heatplot, width=10, height=8)

    }
  }
}


sircleORAHuman(sircleFileName, "entrezgene_id", '_RCM', 'RG2_Changes_filtered', emptyRegLabel="None",
               minGSSize=10, qvalueCutoff=0.1, outputFolder=output_data_dir)


sircleFileName <- file.path(input_data_dir, paste0('Protein_', cancer, '_ORA.csv'))
test_title <- '_protein'
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, 'ORA_label', emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, outputFolder=output_data_dir)

sircleFileName <- file.path(input_data_dir, paste0('RNA_', cancer, '_ORA.csv'))
test_title <- '_RNA'
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, 'ORA_label', emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, outputFolder=output_data_dir)

sircleFileName <- file.path(input_data_dir, paste0('DNAMethylation_', cancer, '_ORA.csv'))
test_title <- '_Methylation'
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, 'ORA_label', emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, outputFolder=output_data_dir)

#  ----------------------------------------------------------------------------------------------
#                 Do GSEA on the protein 
#  ----------------------------------------------------------------------------------------------

sircleFileName <- "Input_Methylation/ClearCellRenalCellCarcinoma_filtered_DCpG.csv"
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
vis_df <- read.csv(sircleFileName)

GSEA_results <- GSEASetupSircle(vis_df, "Pathways DMC", gene_name="gene_name", stat='beta_diff',   GSEA_file_dir=GSEA_file_dir, 
                                output_dir=figure_dir)

sircleFileName <- "Input_Protein/ClearCellRenalCellCarcinoma_filtered_DA_Protein.csv"
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
vis_df <- read.csv(sircleFileName)

GSEA_results <- GSEASetupSircle(vis_df, "Pathways Protein", gene_name="gene_name", stat='logFC_protein',   GSEA_file_dir=GSEA_file_dir, 
                                output_dir=figure_dir)


sircleFileName <- "Input_RNAseq/ClearCellRenalCellCarcinoma_filtered_DE_RNA.csv"
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
vis_df <- read.csv(sircleFileName)

GSEA_results <- GSEASetupSircle(vis_df, "Pathways RNA", gene_name="gene_name", stat='logFC_rna',   GSEA_file_dir=GSEA_file_dir, 
                                output_dir=figure_dir)


#  ----------------------------------------------------------------------------------------------
#                 Now look at the Late vs Early dataset and save these
#  ----------------------------------------------------------------------------------------------
output_data_dir <- 'Output_Data/ORA_EarlyLate/' # _ORA_EarlyLate

sircleFileName <- file.path(input_data_dir, paste0('Protein_', cancer, '_ORA_EarlyLate.csv'))
test_title <- '_protein'
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, 'ORA_label', emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, outputFolder=output_data_dir)

sircleFileName <- file.path(input_data_dir, paste0('RNA_', cancer, '_ORA_EarlyLate.csv'))
test_title <- '_RNA'
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, 'ORA_label', emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, outputFolder=output_data_dir)

sircleFileName <- file.path(input_data_dir, paste0('DNAMethylation_', cancer, '_ORA_EarlyLate.csv'))
test_title <- '_Methylation'
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, 'ORA_label', emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, outputFolder=output_data_dir)



#  ----------------------------------------------------------------------------------------------
#                 Now look at the Late vs Early dataset and save these
#  ----------------------------------------------------------------------------------------------
output_data_dir <- 'Output_Data/ORA_Early/' # _ORA_EarlyLate

sircleFileName <- file.path(input_data_dir, paste0('Protein_', cancer, '_ORA_Early.csv'))
test_title <- '_protein'
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, 'ORA_label', emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, outputFolder=output_data_dir)

sircleFileName <- file.path(input_data_dir, paste0('RNA_', cancer, '_ORA_Early.csv'))
test_title <- '_RNA'
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, 'ORA_label', emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, outputFolder=output_data_dir)

sircleFileName <- file.path(input_data_dir, paste0('DNAMethylation_', cancer, '_ORA_Early.csv'))
test_title <- '_Methylation'
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, 'ORA_label', emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, outputFolder=output_data_dir)



#  ----------------------------------------------------------------------------------------------
#                 Now look at the Late vs Early dataset and save these
#  ----------------------------------------------------------------------------------------------
output_data_dir <- 'Output_Data/ORA_Late/' # _ORA_EarlyLate

sircleFileName <- file.path(input_data_dir, paste0('Protein_', cancer, '_ORA_Late.csv'))
test_title <- '_protein'
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, 'ORA_label', emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, outputFolder=output_data_dir)

sircleFileName <- file.path(input_data_dir, paste0('RNA_', cancer, '_ORA_Late.csv'))
test_title <- '_RNA'
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, 'ORA_label', emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, outputFolder=output_data_dir)

sircleFileName <- file.path(input_data_dir, paste0('DNAMethylation_', cancer, '_ORA_Late.csv'))
test_title <- '_Methylation'
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, 'ORA_label', emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, outputFolder=output_data_dir)

#  ----------------------------------------------------------------------------------------------
#                 Do GSEA for each as well so that we can see if there are any coordinated changes
#  ----------------------------------------------------------------------------------------------


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
cancer <- 'ClearCellRenalCellCarcinoma'

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

sircleFileName <- "Input_Methylation/ClearCellRenalCellCarcinoma_filtered_DCpG_EarlyLate.csv"
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
vis_df <- read.csv(sircleFileName)

GSEA_results <- GSEASetupSircle(vis_df, "Early vs Late DMC", gene_name="gene_name", stat='beta_diff',   GSEA_file_dir=GSEA_file_dir, 
                                output_dir=figure_dir)

sircleFileName <- "Input_Protein/ClearCellRenalCellCarcinoma_filtered_DA_Protein_EarlyLate.csv"
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
vis_df <- read.csv(sircleFileName)

GSEA_results <- GSEASetupSircle(vis_df, "Early vs Late Protein", gene_name="gene_name", stat='logFC_protein',   GSEA_file_dir=GSEA_file_dir, 
                                output_dir=figure_dir)


sircleFileName <- "Input_RNAseq/ClearCellRenalCellCarcinoma_filtered_DE_RNA_EarlyLate.csv"
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
vis_df <- read.csv(sircleFileName)

GSEA_results <- GSEASetupSircle(vis_df, "Early vs Late RNA", gene_name="gene_name", stat='logFC_rna',   GSEA_file_dir=GSEA_file_dir, 
                                output_dir=figure_dir)