# Set up single cell files

## Look at immune response markers in single cell data

The single cell data comes from: 
https://www.sciencedirect.com/science/article/pii/S1535610821001653?via%3Dihub 
- Single-cell RNA-seq reveals the architecture of the ccRCC immune microenvironment  
- Multiregional immune phenotypes integrated with bulk RNA-seq and tumor pathology  
- TCR usage varies by phenotype and defines T cell differentiation trajectories  
- Signatures of tissue-resident T cells and TAMs predict clinical outcome  

#### Summary from single cell paper:  
Clear cell renal cell carcinomas (ccRCCs) are highly immune infiltrated, but the effect of immune heterogeneity on clinical outcome in ccRCC has not been fully characterized. Here we perform paired single-cell RNA (scRNA) and T cell receptor (TCR) sequencing of 167,283 cells from multiple tumor regions, lymph node, normal kidney, and peripheral blood of two immune checkpoint blockade (ICB)-naïve and four ICB-treated patients to map the ccRCC immune landscape. We detect extensive heterogeneity within and between patients, with enrichment of CD8A+ tissue-resident T cells in a patient responsive to ICB and tumor-associated macrophages (TAMs) in a resistant patient. A TCR trajectory framework suggests distinct T cell differentiation pathways between patients responding and resistant to ICB. Finally, scRNA-derived signatures of tissue-resident T cells and TAMs are associated with response to ICB and targeted therapies across multiple independent cohorts. Our study establishes a multimodal interrogation of the cellular programs underlying therapeutic efficacy in ccRCC.

#### Figure 2 
Mainly we're interested in the UMAP from Figure 2: Immune landscape of patients with ccRCC at single-cell resolution

(A) UMAP embedding of transcriptional profiles from all patients and samples (n = 167,283). Each dot represents a single cell, and colors represent clusters denoted by inferred cell type.  


#### Data details
Data were downloaded from https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?analysis=SRZ190804 

```
ccRCC_6pat_Seurat	23.4Gb	2021-03-26 10:51:14	c26b3e222bcf6d1f40ef9112864bc9f7
ccRCC_6pat_cell_annotations.txt	20.2Mb	2021-03-26 10:54:30	9c8f646ec8a1f54fa97c2f1524d505b3
ccRCC_TCRs.txt	3.7Mb	2021-03-26 10:54:31	204982b821836f39950627a642218f9f
ccRCC_regions.txt	348b	2021-03-26 10:54:32	af230043d451b477dc44386a0d26d075
```

## Format the single cell data from the paper
```

library(Seurat)
library(SeuratData)
library(loomR)
library(tidyverse)
devtools::install_github(repo = 'hhoeflin/hdf5r') # Update hdf5r for stability improvements
devtools::install_github(repo = 'mojaveazure/loomR', ref = 'develop') # Update loomR for up-to-date functions

install.packages('SeuratDisk')

singleCell <- readRDS('../data/raw_downloads/single_cell/ccRCC_6pat_Seurat')
output_folder <- '../data/raw_downloads/single_cell/'
kidney_cc <- UpdateSeuratObject(object = singleCell)

FeaturePlot(kidney_cc, features = c("MET", "BAP1", "SETD2"),  reduction = "reductions.tsne", cols = c("lightgrey", "darkred"), ncol = 3) & theme(plot.title = element_text(size = 10))

tsne_red <- kidney_cc@reductions[['tsne']]
rna_counts <- kidney_cc@assays[["RNA"]]@counts
library('Matrix')
writeMM(rna_counts, file=paste0(output_folder, 'rna_counts_ccRCC_6pat_Seurat.txt'))

gene_ids <- kidney_cc@assays[["RNA"]]@counts@Dimnames[[1]]
write.csv(gene_ids, paste0(output_folder, 'geneids_ccRCC_6pat_Seurat.csv'))

meta_data <- kidney_cc@meta.data
meta_data <- as.data.frame(meta_data)
write.csv(meta_data, 'metadata_ccRCC_6pat_Seurat.csv')

tsne_red <- kidney_cc@reductions[['tsne']]
tsne_red <- as.data.frame(tsne_red@cell.embeddings)
write.csv(tsne_red, paste0(output_folder, 'tsne_red_ccRCC_6pat_Seurat.csv'))

umap_red <- kidney_cc@reductions[['umap']]
umap_red <- as.data.frame(umap_red@cell.embeddings)
write.csv(umap_red, paste0(output_folder, 'umap_red_ccRCC_6pat_Seurat.csv'))

SaveH5Seurat(kidney_cc, filename = paste0(output_folder, "kidney_cc.h5Seurat"))
Convert(paste0(output_folder, "kidney_cc.h5Seurat"), dest = "h5ad")

```


