# Imputation script from: https://github.com/WangLab-MSSM/DreamAI
require("cluster")
require("survival")
require("randomForest")
require("missForest")
require("glmnet")
require("Rcpp")
require("foreach")
require("itertools")
require("iterators")
require("Matrix")
require("devtools")
#BiocManager::install("impute") #, version = "3.8")
require("impute")
require("remotes")
#install_github("WangLab-MSSM/DreamAI/Code")
library("DreamAI")
prot_data <- read.csv('../data/S050_CCRCC_Clark_Cell2019/downloaded_data/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB.tsv', sep='\t')
colnames(prot_data)
prot_num_data <- prot_data[, 5:length(colnames(prot_data))]
rownames(prot_num_data) <- prot_data$Index
imputed_data <- DreamAI(prot_num_data, k = 10, maxiter_MF = 10, ntree = 100,
                        maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2),
                        gamma_ADMIN = NA, gamma = 50, CV = FALSE,
                        fillmethod = "row_mean", maxiter_RegImpute = 10,
                        conv_nrmse = 1e-06, iter_SpectroFM = 40, method = c("KNN",
                                                                            "MissForest", "ADMIN", "Birnn", "SpectroFM", "RegImpute"),
                        out = c("Ensemble"))

ens_data <- imputed_data$Ensemble
# Save the imputed data to a csv for the remaining analysis
write.csv(ens_data, '../data/S050_CCRCC_Clark_Cell2019/sircle/prot_df_imputed.csv')

# Do the same for the new dataset
prot_data <- read.csv('../data/ccRCC proteome quantification matrix.csv', sep=',')
colnames(prot_data)
prot_num_data <- prot_data[, 3:length(colnames(prot_data))]
rownames(prot_num_data) <- prot_data$Genes.ID
imputed_data <- DreamAI(prot_num_data, k = 10, maxiter_MF = 10, ntree = 100,
                        maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2),
                        gamma_ADMIN = NA, gamma = 50, CV = FALSE,
                        fillmethod = "row_mean", maxiter_RegImpute = 10,
                        conv_nrmse = 1e-06, iter_SpectroFM = 40, method = c("KNN", "ADMIN", "Birnn", "SpectroFM", "RegImpute"),
                        out = c("Ensemble"))

ens_data <- imputed_data$Ensemble
# Save the imputed data to a csv for the remaining analysis
write.csv(ens_data, '../data/prot_df_imputed_new_dataset_more.csv')
