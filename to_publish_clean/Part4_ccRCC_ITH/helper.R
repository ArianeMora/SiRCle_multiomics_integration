#### Visualisation and differential analysis helper functions

library(EnhancedVolcano)
library(org.Hs.eg.db)#For Human load:
library(clusterProfiler)#To run ORA
library(enrichplot)#For emapplot, dotplot,...
library(ggplot2)#For saving and legends
library(tidyverse) # used for data manipulation
#Establish a couple of functions needed to perform the GSEA (made by Aurelien)
library(fgsea)
library(GSEABase)
library(ggrepel)
library(missMethyl)
library(minfi)

# 

#### Volcano plot for visualising the data
VolcanoPlotSircle <- function(input_data, gene_name, pvalue="padj", log_fc="Log2FC", p_cutoff=0.05, 
                              logfc_cutoff=1.0, subtitle="subtitle", title="T", col_custom = NULL, x_label=NULL, y_label=NULL,
                              colours=c("black", "cyan4", "cyan4", "red3"), output_dir=NULL, point_size=4, label_size=2.0) {
  
  # Set output dir to be the current directory if it is null otherwise set 
  if (is_null(output_dir)) {
    output_dir <- getwd()
  }
  if (is_null(xlab)) {
    x_label <- bquote(~Log[2]~ "FC")
  }
  if (is_null(ylab)) {
    y_label <- bquote(~-Log[10]~p.adj)
  }
  
  ylim <- c(0, ceiling(max(-log10(input_data[[pvalue]]))))  # make the max the biggest value we see
  xlim <- c(floor(min(log2(input_data[[log_fc]]))), ceiling(max(log2(input_data[[log_fc]]))))  # make the max the biggest value we see
  VolcanoPlot<- EnhancedVolcano (input_data,
                                 lab = input_data[[gene_name]],
                                 x = log_fc, #Log2FC
                                 y = pvalue, #p-value or q-value
                                 xlab = x_label,
                                 ylab = y_label,#(~-Log[10]~adjusted~italic(P))
                                 pCutoff = p_cutoff,
                                 FCcutoff = logfc_cutoff, # Cut off Log2FC, automatically 2
                                 pointSize = point_size,
                                 labSize = label_size,
                                 titleLabSize = 16,
                                 col=colours,
                                 colCustom=col_custom,
                                 colAlpha = 0.5,
                                 title=title,
                                 xlim=xlim,
                                 ylim=ylim,
                                 subtitle = subtitle,
                                 caption = paste0("total = ", nrow(input_data)),
                                 cutoffLineType = "dashed",
                                 cutoffLineCol = "black",
                                 cutoffLineWidth = 0.5,
                                 legendPosition = 'right',
                                 legendLabSize = -1,
                                 legendIconSize = -1
  )
  ggsave(file=file.path(output_dir, gsub(" ", "_", paste("VolcanoPlot_", title,  "_", subtitle, ".pdf", sep=""))), plot=VolcanoPlot, width=10, height=8)
}


GoSircle <- function(df, gene_name="entrezgene_id", title="GoResults", key_type="ENTREZID",
                     logfc='log2FC', pvalue='adj.P.Val', 
                     logfc_cutoff=0.5, p_cutoff=0.05, q_cutoff=0.1, output_dir=NULL) {
  # Set output dir to be the current directory if it is null otherwise set 
  if (is_null(output_dir)) {
    output_dir <- getwd()
  }
  
  # If the fold change is greater than 0 then we want a positive cutoff.
  if (logfc_cutoff > 0) {
    selected_genes <- subset(df, df[[pvalue]] < p_cutoff & df[[logfc]]  > logfc_cutoff)
  } else {
    selected_genes <- subset(df, df[[pvalue]] < p_cutoff & df[[logfc]] < logfc_cutoff)  # Otherwise we're testing for a negative cutoff
  }
  
  filename <- gsub(" ", "_", paste0(title, "_LogFC", logfc_cutoff, "_pvalue", p_cutoff))
  go_res <- enrichGO(gene = unlist(selected_genes[[gene_name]]), 
                     universe = df[[gene_name]],
                     keyType = key_type,
                     OrgDb = org.Hs.eg.db, 
                     ont = "ALL",
                     pAdjustMethod = "BH", 
                     qvalueCutoff = q_cutoff, 
                     readable = TRUE)
  go_res_csv <- data.frame(go_res)
  write_csv(go_res_csv, file.path(output_dir, paste0(filename, ".csv")))
  
  ouput_filename <- file.path(output_dir, filename)
  x <- tryCatch(
    {
      Dotplot <- dotplot(go_res, showCategory=30) + ggtitle(title)+ theme(axis.text.x =element_text(size=8,face="bold"), axis.text.y =element_text(size=8,face="bold"), axis.title=element_text(size=5))
      ggsave(paste0(ouput_filename, "_Dotplot.pdf"), plot=Dotplot, width=10, height=8)
      #https://github.com/YuLab-SMU/enrichplot/issues/22
      
      Heatplot <- heatplot(go_res,showCategory=30) + theme(axis.text.x =element_text(size=2),
                                                           axis.text.y =element_text(size=5,face="bold"), axis.title=element_text(size=12,face="bold")) + ggtitle(title)
      ggsave(paste0(ouput_filename, "_Heatplot.pdf"), plot=Heatplot, width=10, height=8)
      
      x2 <- pairwise_termsim(go_res)
      Emapplot<- emapplot(x2, pie_scale=1.5, showCategory=30, layout = "nicely", cex_label_category=0.4, min_edge=0.2)+ ggtitle(title)
      ggsave(paste0(ouput_filename, "_Emapplot.pdf"), plot=Emapplot, width=10, height=8)
    },
    error = function(e){
      print(e)
    }
  )
}


GSEAPathwaysVisSircle <- function (GSEA_result, logfc_cutoff=0.5, p_cutoff=0.25, title="Volcano", 
                                   subtitle="GSEA", output_dir=NULL, plot_type="") {
  # Set output dir to be the current directory if it is null otherwise set 
  if (is_null(output_dir)) {
    output_dir <- getwd()
  }
  
  if (plot_type == "pathways") {
    Volcano1 <- separate(GSEA_result, "pathway", into = c("signature", "rest"), sep= "_", remove=FALSE)%>% 
      mutate(colour = case_when(signature =="KEGG" ~ 'blue',
                                signature =="BIOCARTA" ~ 'gold4',
                                signature =="HALLMARK" ~ 'deeppink4',
                                signature =="REACTOME" ~ 'seagreen4',
                                TRUE ~ 'Not_Detected'))
    #Prepare new colour scheme:
    keyvals <- ifelse(
      Volcano1$colour == "blue", "blue",
      ifelse(Volcano1$colour == "gold4", "gold4",
             ifelse(Volcano1$colour == "deeppink4", "deeppink4",
                    ifelse(Volcano1$colour == "seagreen4", "seagreen4",
                           "black"))))
    keyvals[is.na(keyvals)] <- 'black'
    names(keyvals)[keyvals == 'blue'] <- "KEGG"
    names(keyvals)[keyvals == 'gold4'] <- "BIOCARTA"
    names(keyvals)[keyvals == 'deeppink4'] <- "HALLMARK"
    names(keyvals)[keyvals == 'seagreen4'] <- "REACTOME"
    names(keyvals)[keyvals == 'black'] <- 'X'
    
    VolcanoPlotSircle(Volcano1, "pathway", pvalue="padj", log_fc="NES", x_label=bquote(~Log[2]~ "NES"), p_cutoff=p_cutoff, col_custom = keyvals,
                      logfc_cutoff=logfc_cutoff, subtitle=subtitle, title=title, output_dir=output_dir, point_size=6.0, label_size=3.0)
  } else {
    VolcanoPlotSircle(GSEA_result, "pathway", pvalue="padj", log_fc="NES", p_cutoff=p_cutoff, x_label=bquote(~Log[2]~ "NES"),
                      logfc_cutoff=logfc_cutoff, subtitle=subtitle, title=title, output_dir=output_dir, point_size=6.0, label_size=3.0)
  }
}


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
  pathways <- rbind(GO_BP, GO_CC, GO_MF)
  pathway_list <- list()
  for(pathway in unique(pathways$term))
  {pathway_list[[pathway]] <- as.character(pathways[pathways$term == pathway, 1])}
  
  gsea_result2 <- fgsea(pathways = pathway_list, stats = t_val, nperm = 10000, minSize=4)
  write_csv(gsea_result2, paste0(filename,  "_GSEA_GO-terms.csv"))
  
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
  Column_Gene <- as.character(GSEA_data[[gene_name]])
  Column_tval <- as.numeric(GSEA_data[[stat]])
  MyData_Extracted <- data.frame(cbind(Column_Gene, Column_tval), stringsAsFactors = F)
  MyData_Extracted$Column_tval <- as.numeric(MyData_Extracted$Column_tval)
  
  t_val <- as.numeric(MyData_Extracted$Column_tval)#Make the data into a vector
  names(t_val) <- MyData_Extracted$Column_Gene
  
  Correction <- read.csv(file.path(GSEA_file_dir, "41467_2016_BFncomms13041_MOESM341_ESM.csv"))
  MyData_Extracted_Corr <- merge(x=MyData_Extracted, by.x="Column_Gene", y=Correction, by.y="gene")
  MyData_Extracted_Corr$Column_tval_Corr <- (MyData_Extracted_Corr$Column_tval)/(MyData_Extracted_Corr$pathway)
  
  # matrix for GSEA:
  t_val <- as.numeric(MyData_Extracted_Corr$Column_tval_Corr)#Make the data into a vector
  names(t_val) <- MyData_Extracted_Corr$Column_Gene
  
  #3 Metabolic enzymes
  pathways <- Metabolic  # Ariane update
  pathway_list <- list()
  for(pathway in unique(pathways$term))
  {pathway_list[[pathway]] <- as.character(pathways[pathways$term == pathway, 1])}
  
  gsea_result3_Corrected <- fgsea(pathways = pathway_list, stats = t_val, nperm = 10000, minSize=4)
  write_csv(gsea_result3_Corrected, paste0(filename, "_GSEA_MetabolicPathways_PromiscuityCorrected.csv"))
  
  # Return each of the GSEA results
  #return(list("pathways" = gsea_result1, "GO" = gsea_result2, "MetabolicPathways" = gsea_result3, "MetabolicPathwaysCorrected" =  gsea_result3_Corrected))
}


library(DESeq2)
library(edgeR)

pairedPatientDE <- function(data_file, sample_file, output_file, paired=FALSE) {
  cat(paste("Differential Expression analysis for: \n", data_file, "\n"))
  
  counts <- read.csv(data_file, header = TRUE, sep = ",")
  rownames(counts) <- counts$gene_id
  
  experimental_design <- read.csv(sample_file)
  other_info <- counts[, c('gene_name', 'gene_id')] # Keep the metadata info
  rownames(other_info) <- rownames(counts)
  # Let's make sure our count data is in matrix format and is only the numeric columns i.e. everything but the genes
  count_matrix <- as.matrix(counts[,experimental_design$FullLabel])
  
  count_matrix[is.na(count_matrix)] = 0
  # We now set the row names to be the gene IDs
  rownames(count_matrix) <- rownames(counts) 
  
  # Separate out each of the sample names to be the different experiment conditions
  condition_id <- as.factor(experimental_design$CondId)
  case_id <- as.factor(experimental_design$SafeCases)
  
  # Before getting the means we want to normalise the counts (i.e. so our mean isn't stupidly big)
  dge <- DGEList(counts=count_matrix)
  dge <- calcNormFactors(dge, method="TMM")
  tmm <- cpm(dge)
  
  sample_df = data.frame(case_id = case_id, condition_id=condition_id)
  if (paired == TRUE) {
    dds_mat <- DESeqDataSetFromMatrix(countData = count_matrix,
                                      colData = sample_df,
                                      design = ~case_id + condition_id) # Have patient case as a factor (i.e. a batch correction)
    
  } else {
    dds_mat <- DESeqDataSetFromMatrix(countData = count_matrix,
                                      colData = sample_df,
                                      design = ~condition_id) # Just do on condition
    
  }
  
  dds <- estimateSizeFactors(dds_mat)
  
  num_samples_meeting_criteria <- 4  # be strict and enforce that at least half the samples need to meet the criteria (i.e. one full condition)
  num_counts_in_gene <- 10  # They need at least 10 counts
  
  keep <- rowSums(counts(dds_mat) >= num_counts_in_gene) >= num_samples_meeting_criteria
  dds <- dds_mat[keep,] # Only keep the rows with this criteria
  # Now we want to also filter out our other info using the same filter
  other_info_filtered <- other_info[keep,]
  counts_filtered <- counts[keep,]
  tmm_filtered <- tmm[keep, ]
  
  # Log2 the TMM counts for better variance & mean estimation.
  tmm_filtered <- log2(tmm_filtered + 1)
  # Let's print the number of rows
  cat(paste("Dataset dimensions: ", nrow(dds), ncol(dds), "\n"))
  # Run DEseq2
  dds <- DESeq(dds)
  resultsNames(dds)
  # Build results table
  res <- results(dds, independentFiltering=FALSE)
  other_info_filtered <- other_info_filtered[rownames(res), c('gene_name', 'gene_id')]
  
  # Ensure ordering is correct 
  tmm_filtered <- tmm_filtered[rownames(res), colnames(tmm_filtered)]
  other_info_filtered$logFC_rna <- res$log2FoldChange
  other_info_filtered$stat_rna <- res$stat
  other_info_filtered$pvalue_rna <- res$pvalue
  other_info_filtered$padj_rna <- res$padj
  other_info_filtered$lfcSE_rna <- res$lfcSE
  other_info_filtered$baseMean_rna <- res$baseMean
  other_info_filtered$var_rna <- matrixStats::rowVars(as.matrix(tmm_filtered))
  
  # Add in mean info
  all_rna_df <- cbind(other_info_filtered, tmm_filtered)
  
  write.csv(all_rna_df, file = output_file)
  return(all_rna_df)
}

pairedPatientDA <- function(data_file, sample_file, output_file, gene='Gene', FullLabel='FullLabel', CondId='SampleID', CaseId='SafeCases', 
                            paired=FALSE) {
  cat(paste("Differential Abundence analysis for: \n", data_file, "\n"))
  # Do protein analysis on the stage 4 sample
  prot_data <- read.csv(data_file)
  gene_names <- prot_data[[gene]]
  experimental_design <- read.csv(sample_file)
  data_columns <- experimental_design[[FullLabel]] # get LFQ column numbers
  
  prot_data_mat <- prot_data[, data_columns]
  prot_data_mat[is.na(prot_data_mat)] <- 0
  
  rownames(prot_data_mat) <- gene_names
  rownames(prot_data) <- gene_names
  
  # Get the column of interest from the experimental design
  cond <- as.factor(experimental_design[[CondId]])
  case_id <- as.factor(experimental_design[[CaseId]])
  
  # Create a mean and variance columns that we'll add to the final dataframe for both normal and tumour 
  # We use this later in the SIRCLE model
  tumors <- subset(experimental_design, CondId == 1)
  normals <- subset(experimental_design, CondId == 0)
  tumor_mean <- rowMeans(prot_data_mat[, tumors[[FullLabel]]])
  
  ##Go through each row and determine if a value is zero
  prot_data_mat = prot_data_mat[tumor_mean != 0, ]
  ##Subset as usual
  normal_mean <- rowMeans(prot_data_mat[, normals[[FullLabel]]])
  prot_data_mat = prot_data_mat[normal_mean != 0, ]
  prot_data_mat[is.na(prot_data_mat)] <- 0
  
  if (paired == TRUE) {
    design <- model.matrix(~case_id + cond) # We want to perform it on each patient (i.e. since we have matched samples)
  } else {
    design <- model.matrix(~cond) # We didn't have matching patient info
  }
  # Limma is good for detecting differentially abundent proteins
  fit <- lmFit(prot_data_mat, design)
  fit2 <- eBayes(fit, robust=TRUE) # Use robust:  https://www.biostars.org/p/496806/, https://support.bioconductor.org/p/118495/
  
  # Keep in mind the condition of interest (tumour vs normal) is always the final column in the design
  fit2_adj <- topTable(fit2, coef=length(colnames(design)), adjust="fdr", sort.by="none", number=1000000) 
  
  # Create the dataframe to return
  all_prot_df <- prot_data[rownames(fit2_adj), colnames(prot_data)]
  prot_data_mat <- prot_data_mat[rownames(fit2_adj), colnames(prot_data_mat)]
  # Add in the statistics from the DA analysis
  all_prot_df$logFC_protein <- fit2_adj$logFC
  all_prot_df$stat_protein <- fit2_adj$t
  all_prot_df$pvalue_protein <- fit2_adj$P.Value
  all_prot_df$padj_protein <- fit2_adj$adj.P.Val
  all_prot_df$B_protein <- fit2_adj$B
  all_prot_df$mean_protein <- fit2_adj$AveExpr
  
  write.csv(all_prot_df, output_file, row.names = FALSE)
  return(all_prot_df)
}

library(DMRcate)
library(limma)
library(edgeR)
#library("regRCPqn") # https://github.com/regRCPqn/regRCPqn install using devtools

#https://rdrr.io/bioc/lumi/src/R/methylation_preprocessing.R
# convert beta-value to m-value
beta2m <- function(beta) {
  m <- log2(beta/(1-beta))
  return(m)
}

# convert m-value to beta-value 
m2beta <- function(m) {
  beta <- 2^m/(2^m + 1)
  return(beta)
}

pairedPatientDMC <- function(data_file, sample_file, output_file, label, array_type='EPIC', project_dir='', paired=FALSE) {
  
  cat(paste("Differential Methylation analysis for: \n", data_file, "\n"))
  
  #### Import data ####
  cpg_raw <- read.csv(data_file, header = TRUE, sep = ",")
  sample_df <- read.csv(sample_file)
  
  #### Change rownames ####
  rowNames <- unlist(cpg_raw['Locus'])
  cpg_data <- cpg_raw[, sample_df$FullLabel]
  rownames(cpg_data) <- rowNames
  rownames(cpg_raw) <- rowNames
  #### QC/Filtering
  # First remove all NA values
  cpg_data[is.na(cpg_data)] <- 0
  summary(cpg_data)
  
  # Convert data to a matrix & plot the means for the rows
  cpg_data_m <- as.matrix(cpg_data)
  row_means_data <- rowMeans(cpg_data_m)
  #hist(row_means_data)
  #nrow(cpg_data_m)
  
  # Remove rows with a very small amount of overall DNA methylation
  cpg_data_filtered <- cpg_data_m[row_means_data > 0.05, ]
  cpg_raw_filtered <- cpg_raw[row_means_data > 0.05, ] # Make sure we apply the same filter to our CpG names etc
  nrow(cpg_data_filtered)
  row_means_data <- rowMeans(cpg_data_filtered)
  # hist(row_means_data)
  
  # Remove rows with overly high DNA methylation
  cpg_data_filtered <- cpg_data_filtered[row_means_data < 0.95, ]
  cpg_raw_filtered <- cpg_raw_filtered[row_means_data < 0.95, ] # Make sure we apply the same filter to our CpG names etc
  
  nrow(cpg_data_filtered)
  
  # The function model.matrix is used to generate the design matrix
  cases <- as.factor(sample_df$SafeCases)
  cond_id <- as.factor(sample_df$CondId)
  
  if (paired == FALSE) {
    design = model.matrix(~cond_id)  # Can't do paired in DNA methylation
  } else {
    design = model.matrix(~cases + cond_id) # This way we're actually using the TCGA or the CPTAC as the batch
  }
  
  # Before running limma we want to do two steps, (following the steps of Miss Methyl)
  # 1) convert to M values 
  # 2) perform limma
  # References: 
  # https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
  # https://bioconductor.org/packages/release/bioc/vignettes/missMethyl/inst/doc/missMethyl.html 
  # Add a very small amount so that if we have 0's we don't get an issue with NAs and log2
  cpg_data_M_filtered <- cpg_data_filtered #log2((cpg_data_filtered + 0.000001) /(1-cpg_data_filtered + 0.000001))
  # Normalise M values using https://github.com/regRCPqn/regRCPqn, https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0229763
  # Install: https://github.com/regRCPqn/regRCPqn
  coef_id <- length(colnames(design))
  fit <- lmFit(cpg_data_M_filtered, design)
  fit2 <- eBayes(fit, robust=TRUE)
  
  # Don't sort it otherwise we won't be able to match the data
  fit2_adj <- as.data.frame(topTable(fit2, coef=coef_id, adjust="fdr", sort.by="none", number=1000000))
  all_cpg_df <- cpg_raw[rownames(fit2_adj), c('Locus')]
  
  # Add in the statistics from the CpG analysis
  all_cpg_df$beta_logFC_meth <- fit2_adj$logFC
  all_cpg_df$beta_stat_meth <- fit2_adj$t
  all_cpg_df$beta_pvalue_meth <- fit2_adj$P.Value
  all_cpg_df$beta_padj_meth <- fit2_adj$adj.P.Val
  all_cpg_df$beta_B_meth <- fit2_adj$B
  all_cpg_df$beta_mean_cpg <- fit2_adj$AveExpr
  all_cpg_df$beta_var_cpg <- matrixStats::rowVars(as.matrix(cpg_data_M_filtered))
  write.csv(fit2_adj, paste0(project_dir, label, "_Beta_DMCLimma.csv"))

  ### Do the DMR comparison on the beta values
  out <- tryCatch(
    expr = {
      myAnnotation <- cpg.annotate(object = cpg_data_M_filtered, datatype = "array", what = "Beta", 
                                   arraytype = c(array_type), 
                                   analysis.type = "differential", design = design, coef = coef_id)
      
      DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
      
      results.ranges <- extractRanges(DMRs, genome='hg38')
      
      #### Write DMR results and information to a csv file ####
      write.csv(results.ranges, paste0(project_dir, label, "_Beta_DMRcate.csv"))
    },
    error = function(e){
      message('Caught an error!')
      print(e)
    },
    warning = function(w){
      message('Caught an warning!')
      print(w)
    },
    finally = {
      message('All done, quitting.')
    }
  )
  
  # Do the same with the M values 
  cpg_data_M_filtered <- log2(cpg_data_filtered/(1-cpg_data_filtered))

  # Normalise M values using https://github.com/regRCPqn/regRCPqn, https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0229763
  # Install: https://github.com/regRCPqn/regRCPqn
  # cpg_data_M_filtered <- regRCPqn(cpg_data_M_filtered, "", "", save_ref=FALSE) # Normalising data 
  coef_id <- length(colnames(design))
  fit <- lmFit(cpg_data_M_filtered, design)
  fit2 <- eBayes(fit, robust=TRUE)
  fit2_adj <- as.data.frame(topTable(fit2, coef=coef_id, adjust="fdr", sort.by="none", number=1000000))
  
  # Add in the statistics from the CpG analysis
  all_cpg_df$M_logFC_meth <- fit2_adj$logFC
  all_cpg_df$M_betadiff_meth <- (2^fit2_adj$logFC)/(2^fit2_adj$logFC + 1)
  all_cpg_df$M_stat_meth <- fit2_adj$t
  all_cpg_df$M_pvalue_meth <- fit2_adj$P.Value
  all_cpg_df$M_padj_meth <- fit2_adj$adj.P.Val
  all_cpg_df$M_B_meth <- fit2_adj$B
  all_cpg_df$M_mean_cpg <- fit2_adj$AveExpr
  all_cpg_df$M_var_cpg <- matrixStats::rowVars(as.matrix(cpg_data_M_filtered))
  
  write.csv(fit2_adj, paste0(project_dir, label, "_M_DMCLimma.csv"))
  
  ### Do the DMR comparison
  out <- tryCatch(expr = {
    myAnnotation <- cpg.annotate(object = cpg_data_M_filtered, datatype = "array", what = "M", 
                                 arraytype = c(array_type), 
                                 analysis.type = "differential", design = design, coef = coef_id)
    
    DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
    
    results.ranges <- extractRanges(DMRs, genome='hg38')
    
    #### Write DMR results and information to a csv file ####
    write.csv(results.ranges, paste0(project_dir, label, "_M_DMRcate.csv"))
  },
  error = function(e){
    message('Caught an error!')
    print(e)
  },
  warning = function(w){
    message('Caught an warning!')
    print(w)
  },
  finally = {
    message('All done, quitting.')
  })
  return(all_cpg_df)
}


### Running each of the plots for the DE analyses
runDEplots <- function(data, title, gene_name, entrez_id, logfc_cutoff, p_cutoff, q_cutoff, logfc,  pvalue, stat, output_dir) {
  # Volcano of the data
  # Now run the charts for this (make sure we don't have any NA entrez IDs)
  
  VolcanoPlotSircle(data, gene_name, pvalue=pvalue, log_fc=logfc, p_cutoff=p_cutoff, 
                    logfc_cutoff=logfc_cutoff, title=title, subtitle="Data", output_dir=output_dir)
  
  # Do enrichment for positive and negative logFC
  GoSircle(data, gene_name=entrez_id, title=paste(title, 'Positive'), key_type="ENTREZID", logfc=logfc, pvalue=pvalue, 
           logfc_cutoff=logfc_cutoff, p_cutoff=p_cutoff, q_cutoff=q_cutoff, output_dir=output_dir)
  
  GoSircle(data, gene_name=entrez_id, title=paste(title, 'Negative'), key_type="ENTREZID", logfc=logfc, pvalue=pvalue, 
           logfc_cutoff=-1 * logfc_cutoff, p_cutoff=p_cutoff, q_cutoff=q_cutoff, output_dir=output_dir)
  
}


runDEAll <- function(test_title, input_data_dir, output_data_dir, GSEA_file_dir, paired=FALSE, output_figure_dir="", pathway_list=NA) {
  cat('-------------------------------------------------------------------------------------\n')
  cat(paste0("Running... \n", test_title, "\n"))
  prot_df <- read.csv(file.path(input_data_dir, paste0('prot_data_', test_title, '_sircle.csv')))
  
  prot_df <- pairedPatientDA(file.path(input_data_dir, paste0('prot_data_', test_title, '_sircle.csv')),
                             file.path(input_data_dir, paste0('prot_sample_data_', test_title, '_sircle.csv')),
                             file.path(output_data_dir, paste0('prot_DE_', test_title, '_sircle.csv')), paired=paired)
  
  # Constant between all analyese
  gene_name <- "external_gene_name"
  entrez_id <- "entrezgene_id"
  logfc_cutoff <- 0.5
  p_cutoff <- 0.05
  q_cutoff <- 0.1
  
  # Specific for the protein DE
  title <- paste0(test_title, ' Protein')
  logfc <- "logFC_protein"
  pvalue <- "padj_protein"
  stat <- "stat_protein"
  
  # Now run the charts for this (make sure we don't have any NA entrez IDs)
  prot_df <- subset(prot_df, !is.na(prot_df$entrezgene_id))
  
  runDEplots(prot_df, title, gene_name, entrez_id, logfc_cutoff, p_cutoff, q_cutoff, logfc, pvalue, stat, output_figure_dir)
  
  GSEA_results <- GSEASetupSircle(prot_df, title, gene_name=gene_name, stat=stat, GSEA_file_dir=GSEA_file_dir, output_dir=output_figure_dir, pathway_list=pathway_list)
  
  ##########   RNAseq analysis
  
  rna_df <- pairedPatientDE(file.path(input_data_dir, paste0('rna_data_', test_title, '_sircle.csv')),
                            file.path(input_data_dir, paste0('rna_sample_data_', test_title, '_sircle.csv')),
                            file.path(output_data_dir, paste0('rna_DE_', test_title, '_sircle.csv')), paired=paired)
  
  title <- paste0(test_title, ' RNAseq')
  logfc <- "logFC_rna"
  pvalue <- "padj_rna"
  gene_name <- "external_gene_name"
  stat <- 'stat_rna'
  logfc_cutoff <- 1.0
  # put all the GSEA results in one directory
  output_dir <- file.path(output_data_dir, 'Figures', 'RNAseq')
  rna_df <- subset(rna_df, !is.na(rna_df$entrezgene_id))
  rna_df[is.na(rna_df)] <- 0 # Replace NAs
  rna_de_df <- subset(rna_df, !is.na(rna_df$entrezgene_id))
  rna_de_df <- subset(rna_de_df, !rna_de_df$padj_rna == 0)
  rna_de_df <- subset(rna_de_df, !rna_de_df$pvalue_rna == 0)
  x <- tryCatch(
    {
      runDEplots(rna_de_df, title, gene_name, entrez_id, logfc_cutoff, p_cutoff, q_cutoff, logfc, pvalue, stat, output_figure_dir)
    },
    error = function(e){
      print(e)
    }
  )
  # GSEA section (need to use external gene name here)
  GSEA_results <- NULL
  x <- tryCatch(
    {
      GSEA_results <- GSEASetupSircle(rna_df, title, gene_name=gene_name, stat=stat, GSEA_file_dir=GSEA_file_dir, output_dir=output_figure_dir, pathway_list=pathway_list)
    },
    error = function(e){
      print(e)
    }
  )
  ##########   Methylation analysis
  dmc_df <- pairedPatientDMC(file.path(input_data_dir, paste0('cpg_data_', test_title, '_sircle.csv')),
                             file.path(input_data_dir, paste0('cpg_sample_data_', test_title, '_sircle.csv')),
                             file.path(output_data_dir, paste0('cpg_DE_', test_title, '_sircle.csv')), paired=FALSE)
  
  # Constant between all analyese
  gene_name <- "external_gene_name"
  entrez_id <- "entrezgene_id"
  logfc_cutoff <- 0.1
  p_cutoff <- 0.05
  q_cutoff <- 0.1
  
  # Specific for the protein DE
  title <- paste0(test_title, ' Methylation')
  logfc <- "logFC_meth"
  pvalue <- "padj_meth"
  stat <- "stat_meth"
  
  dmc_df <- subset(dmc_df, !is.na(dmc_df$entrezgene_id))
  dmc_df[is.na(dmc_df)] <- 0 # Replace NAsx
  dmc_de_df <- subset(dmc_df, !is.na(rna_df$entrezgene_id))
  dmc_de_df <- subset(dmc_de_df, !rna_de_df$padj_rna == 0)
  dmc_de_df <- subset(dmc_de_df, !rna_de_df$pvalue_rna == 0)
  x <- tryCatch(
    {
      runDEplots(dmc_de_df, title, gene_name, entrez_id, logfc_cutoff, p_cutoff, q_cutoff, logfc, pvalue, stat, output_figure_dir)
    },
    error = function(e){
      print(e)
    }
  )
  x <- tryCatch(
    {
      ## GSEA section
      GSEA_results <- GSEASetupSircle(dmc_df, title, gene_name=gene_name, stat=stat, GSEA_file_dir=GSEA_file_dir, output_dir=output_figure_dir, pathway_list=pathway_list)
    },
    error = function(e){
      print(e)
    }
  )
  cat(paste0("Done! \n", test_title, "\n\n"))
  cat('-------------------------------------------------------------------------------------\n\n')
}

runDEAll_CPG <- function(test_title, input_data_dir, output_data_dir, GSEA_file_dir, paired=FALSE, output_figure_dir="", pathway_list=NA) {
  
  ##########   Methylation analysis
  dmc_df <- pairedPatientDMC(file.path(input_data_dir, paste0('cpg_data_', test_title, '_sircle.csv')),
                             file.path(input_data_dir, paste0('cpg_sample_data_', test_title, '_sircle.csv')),
                             file.path(output_data_dir, paste0('cpg_DE_', test_title, '_sircle.csv')), paired=paired)
  
  # Constant between all analyese
  gene_name <- "external_gene_name"
  entrez_id <- "entrezgene_id"
  logfc_cutoff <- 0.1
  p_cutoff <- 0.05
  q_cutoff <- 0.1
  
  # Specific for the protein DE
  title <- paste0(test_title, ' Methylation')
  logfc <- "logFC_meth"
  pvalue <- "padj_meth"
  stat <- "stat_meth"
  
  dmc_df <- subset(dmc_df, !is.na(dmc_df$entrezgene_id))
  dmc_df[is.na(dmc_df)] <- 0 # Replace NAsx
  x <- tryCatch(
    {
      runDEplots(dmc_df, title, gene_name, entrez_id, logfc_cutoff, p_cutoff, q_cutoff, logfc, pvalue, stat, output_figure_dir)
    },
    error = function(e){
      print(e)
    }
  )
  x <- tryCatch(
    {
      ## GSEA section
      GSEA_results <- GSEASetupSircle(dmc_df, title, gene_name=gene_name, stat=stat, GSEA_file_dir=GSEA_file_dir, output_dir=output_figure_dir, pathway_list=pathway_list)
    },
    error = function(e){
      print(e)
    }
  )
  cat(paste0("Done! \n", test_title, "\n\n"))
  cat('-------------------------------------------------------------------------------------\n\n')
}


