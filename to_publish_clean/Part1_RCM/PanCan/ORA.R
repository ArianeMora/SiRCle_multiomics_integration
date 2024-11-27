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
cancer <- 'PanCan'

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
