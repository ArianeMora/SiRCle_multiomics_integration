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


project_dir <- ''
cancer <- 'ClearCellRenalCellCarcinoma-New'
data_file <- paste0(project_dir, cancer, '_filtered_RNA.csv')
sample_file <- paste0( project_dir, cancer, '_filtered_samples_RNA.csv') 
output_file <- paste0(project_dir, cancer, '_filtered_DE_RNA.csv')

cat(paste("Differential Expression analysis for: \n", data_file, "\n"))

counts <- read.csv(data_file, header = TRUE, sep = ",")
rownames(counts) <- counts$gene_id

experimental_design <- read.csv(sample_file)
other_info <- counts[, c('gene_name', 'gene_id')] # Keep the metadata info
rownames(other_info) <- rownames(counts)
# Let's make sure our count data is in matrix format and is only the numeric columns i.e. everything but the genes
#nn <- counts[, !(names(counts) %in% experimental_design$Sample)]
experimental_design <- experimental_design[experimental_design$Sample %in% names(counts), ]

count_matrix <- as.matrix(counts[,  experimental_design$Sample])

count_matrix[is.na(count_matrix)] = 0
# We now set the row names to be the gene IDs
rownames(count_matrix) <- rownames(counts) 

# Separate out each of the sample names to be the different experiment conditions
condition_id <- as.factor(experimental_design$CondID)
case_id <- as.factor(experimental_design$SafeCases)
sample_df = data.frame(case_id = case_id, condition_id=condition_id)

# Before getting the means we want to normalise the counts (i.e. so our mean isn't stupidly big)
dge <- DGEList(counts=count_matrix)
dge <- calcNormFactors(dge, method="TMM")
tmm <- cpm(dge)


dds_mat <- DESeqDataSetFromMatrix(countData = count_matrix,
                                  colData = sample_df,
                                  design = ~condition_id) # Just do on condition


dds <- estimateSizeFactors(dds_mat)

# Already did pre-filtering before
dds <- dds_mat
other_info_filtered <- other_info
counts_filtered <- counts
tmm_filtered <- tmm

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
