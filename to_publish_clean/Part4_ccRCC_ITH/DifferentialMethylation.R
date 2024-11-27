#install.packages('tidyverse')
#BiocManager::install("EnhancedVolcano")

#BiocManager::install("clusterProfiler")
#BiocManager::install("fgsea")
#BiocManager::install("enrichplot")
#BiocManager::install("GSEABase")
#BiocManager::install("missMethyl")
#BiocManager::install("minfi")
#BiocManager::install("DESeq2")

  
library("limma")
library("dplyr")
library(matrixStats)
library(DMRcate)
library(DESeq2)
library(edgeR)
library(data.table)

project_dir <- ''
cancer <- 'ClearCellRenalCellCarcinoma-New'
data_file <- paste0(project_dir, cancer, '_filtered_CpG.csv')
sample_file <- paste0( project_dir, cancer, '_filtered_samples_CpG.csv') 
output_file <- paste0(project_dir, cancer, '_filtered_DMC_CpG.csv')

cat(paste("Differential Methylation analysis for: \n", data_file, "\n"))

#### Import data ####
cpg_raw <- read.csv(data_file, header = TRUE, sep = ",")
sample_df <- read.csv(sample_file)

#### Change rownames ####
rowNames <- unlist(cpg_raw['Locus'])
sample_df <- sample_df[sample_df$Sample %in% names(cpg_raw), ]
cpg_data <- cpg_raw[, sample_df$Sample]
rownames(cpg_data) <- rowNames
rownames(cpg_raw) <- rowNames
#### QC/Filtering
# First remove all NA values
cpg_data[is.na(cpg_data)] <- 0
summary(cpg_data)

# Convert data to a matrix & plot the means for the rows
cpg_data_m <- as.matrix(cpg_data)

# The function model.matrix is used to generate the design matrix
cond_id <- as.factor(sample_df$CondID)
cases <- as.factor(sample_df$SafeCases)
design = model.matrix(~cond_id)  # Can't do paired in DNA methylation

# Before running limma we want to do two steps, (following the steps of Miss Methyl)
# 1) convert to M values 
# 2) perform limma
# References: 
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
# https://bioconductor.org/packages/release/bioc/vignettes/missMethyl/inst/doc/missMethyl.html 
# Add a very small amount so that if we have 0's we don't get an issue with NAs and log2
cpg_data_M_filtered <- log2((cpg_data_m) /(1-cpg_data_m))
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

write.csv(fit2_adj, output_file)

