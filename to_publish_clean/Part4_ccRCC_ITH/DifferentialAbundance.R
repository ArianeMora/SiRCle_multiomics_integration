install.packages('tidyverse')
install.packages('remotes')
install.packages('impute')
install.packages('devtools')
BiocManager::install("impute")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")
BiocManager::install("EnhancedVolcano")
BiocManager::install("clusterProfiler")
BiocManager::install("fgsea")
BiocManager::install("enrichplot")
BiocManager::install("GSEABase")
BiocManager::install("missMethyl")
BiocManager::install("minfi")
BiocManager::install("DESeq2")
BiocManager::install("DMRcate")
install_github("WangLab-MSSM/DreamAI/Code")

  
library("limma")
library("dplyr")
library(matrixStats)
library(DMRcate)
source("helper.R")
library(DESeq2)
library(edgeR)

require("impute")

require("devtools")
require("remotes")
# You need to install via the below link

library("DreamAI")


project_dir <- ''
cancer <- 'ClearCellRenalCellCarcinoma-New'
data_file <- paste0(project_dir, cancer, '_filtered_Protein.csv')
sample_file <- paste0( project_dir, cancer, '_filtered_samples_Protein.csv') 
imputed_file <- paste0(project_dir, cancer, '_filtered_imputed_Protein.csv')

output_file <- paste0(project_dir, cancer, '_filtered_DA_Protein.csv')

cat(paste("Differential Protein analysis for: \n", data_file, "\n"))

CondId <- 'CondID'
FullLabel <- 'Sample'
CaseId <- 'SafeCases'
paired=FALSE
gene = 'gene_name'

cat(paste("Differential Abundence analysis for: \n", data_file, "\n"))
# Do protein analysis on the stage 4 sample
prot_data <- read.csv(data_file)
gene_names <- prot_data[[gene]]
experimental_design <- read.csv(sample_file)
data_columns <- experimental_design[[FullLabel]] # get LFQ column numbers
prot_data_mat <- prot_data[, data_columns]

rownames(prot_data_mat) <- gene_names
rownames(prot_data) <- gene_names

prot_num_data <- prot_data
prot_data_mat[prot_data_mat == 0] <- NA
imputed_data <- DreamAI(prot_data_mat, k = 10, maxiter_MF = 10, ntree = 100,
                        maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2),
                        gamma_ADMIN = NA, gamma = 50, CV = FALSE,
                        fillmethod = "row_mean", maxiter_RegImpute = 10,
                        conv_nrmse = 1e-06, iter_SpectroFM = 40, method = c("SpectroFM", "KNN"),
                        out = c("Ensemble"))


ens_data <- imputed_data$Ensemble
write.csv(ens_data, imputed_file)

# Pause here and go and check for outliers and filter!
# ---------------------------- GO BACK TO PYTHON

imputed_file <- paste0(project_dir, cancer, '_filtered_imputed-norm_Protein.csv')

prot_data <- read.csv(imputed_file)
gene_names <- prot_data[['gene_name']]
experimental_design <- read.csv(sample_file)
data_columns <- experimental_design[[FullLabel]] # get LFQ column numbers
prot_data_mat <- prot_data[, data_columns]

rownames(prot_data_mat) <- gene_names
rownames(prot_data) <- gene_names


# Get the column of interest from the experimental design
cond <- as.factor(experimental_design[[CondId]])
case_id <- as.factor(experimental_design[[CaseId]])
design <- model.matrix(~cond) # We didn't have matching patient info
# Limma is good for detecting differentially abundent proteins
fit <- lmFit(prot_data_mat, design)
fit2 <- eBayes(fit, robust=TRUE) # Use robust:  https://www.biostars.org/p/496806/, https://support.bioconductor.org/p/118495/

# Keep in mind the condition of interest (tumour vs normal) is always the final column in the design
fit2_adj <- topTable(fit2, coef=length(colnames(design)), adjust="fdr", sort.by="none", number=1000000) 

# Create the dataframe to return
all_prot_df <- data.frame(prot_data[rownames(fit2_adj), colnames(prot_data)])
prot_data_mat <- prot_data_mat[rownames(fit2_adj), colnames(prot_data_mat)]
# Add in the statistics from the DA analysis
all_prot_df$gene_name <- rownames(fit2_adj)
all_prot_df$logFC_protein <- fit2_adj$logFC
all_prot_df$stat_protein <- fit2_adj$t
all_prot_df$pvalue_protein <- fit2_adj$P.Value
all_prot_df$padj_protein <- fit2_adj$adj.P.Val
all_prot_df$B_protein <- fit2_adj$B
all_prot_df$mean_protein <- fit2_adj$AveExpr

write.csv(all_prot_df, output_file, row.names = FALSE)


