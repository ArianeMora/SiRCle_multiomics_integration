.. _about:

SUMMARY
=======

Preprint: `bioRxiv <X>`_

**********************************************************************************************************************
SiRCle (Signature Regulatory Clustering) model integration reveals mechanisms of phenotype regulation in renal cancer
**********************************************************************************************************************

Ariane Mora^1, Christina Schmidt^2,3, Brad Balderson1, Christian Frezza3#, Mikael Bodén1#

1) School of Chemistry and Molecular Biosciences, University of Queensland, Molecular Biosciences Building 76, St Lucia QLD 4072, Australia.
2) Medical Research Council Cancer Unit, University of Cambridge, Hutchison/MRC Research Centre, Box 197, Cambridge Biomedical Campus, Cambridge CB2 0X2, United Kingdom
3) CECAD Research Center, University Hospital Cologne, Joseph-Stelzmann-Str. 26, 50931 Cologne, Germany

^Joint first authors
#Joint last authors

Note Christina and Ariane are equal joint first authors and the authors may swap the order of their names as they so choose :)

Abstract
========
Clear cell renal cell carcinoma (ccRCC) tumours develop and progress via complex remodelling of the kidney epigenome, transcriptome, proteome, and metabolome. Given the subsequent tumour and inter-patient heterogeneity, drug-based treatments report limited success, calling for multi-omics studies to extract regulatory relationships, and ultimately, to develop targeted therapies. However, current methods are unable to extract nonlinear multi-omics perturbations.

Here, we present SiRCle (Signature Regulatory Clustering), a novel method to integrate DNA methylation, RNA-seq and proteomics data. Applying SiRCle to a case study of ccRCC, we disentangle the layer (DNA methylation, transcription and/or translation) where dysregulation first occurs and find the primary biological processes altered. Next, we detect regulatory differences between patient subsets by using a variational autoencoder to integrate omics’ data followed by statistical comparisons on the integrated space. In ccRCC patients, SiRCle allows to identify metabolic enzymes and cell-type-specific markers associated with survival along with the likely molecular driver behind the gene’s perturbations.


.. figure:: _static/summary.png
   :width: 800
   :align: center

sci-RegulatoryClusteringModel
-----------------------------

`Package <https://github.com/ArianeMora/scircm/>`_

The general table of how we define regulatory clusters.

.. code-block:: python

    """
        | Methylation      | RNAseq    | Proteomics | Regulation driver_1          | Regulation driver_2     | Regulation_Grouping1 | Regulation_Grouping2 | Regulation_Grouping3 |
        |------------------|-----------|------------|------------------------------|-------------------------|----------------------|----------------------|----------------------|
        | Hypermethylation | DOWN      | DOWN       | Methylation increase (MDS)   | None                    | MDS                  | MDS                  | MDS                  |
        | Hypermethylation | UP        | DOWN       | mRNA increase (TPDE)         | Protein decrease (TMDS) | TPDE+TMDS            | TPDE+TMDS            | TMDS                 |
        | Hypermethylation | UP        | UP         | mRNA increase (TPDE)         | None                    | TPDE                 | TPDE                 | TPDE                 |
        | Hypermethylation | DOWN      | UP         | Methylation increase (MDS)   | Protein increase (TMDE) | MDS+TMDE             | TMDE                 | TMDE                 |
        | Hypermethylation | No Change | UP         | mRNA increase (TPDE)         | Protein increase (TMDE) | TPDE+TMDE            | TMDE                 | TMDE                 |
        | Hypermethylation | No Change | DOWN       | mRNA increase (TPDE)         | Protein decrease (TMDS) | TPDE+TMDS            | TMDS                 | TMDS                 |
        | Hypermethylation | UP        | No Change  | mRNA increase (TPDE)         | Protein decrease (TMDS) | TPDE+TMDS            | TPDE+TMDS            | TMDS                 |
        | Hypermethylation | DOWN      | No Change  | Methylation increase (MDS)   | Protein increase (TMDE) | MDS+TMDE             | MDS+TMDE             | TMDE                 |
        | Hypermethylation | No Change | No Change  | Methylation increase (ncRNA) | None                    | MDS-ncRNA            | MDS_ncRNA            | MDS_ncRNA            |
        | Hypomethylation  | DOWN      | DOWN       | mRNA decrease (TPDS)         | None                    | TPDS                 | TPDS                 | TPDS                 |
        | Hypomethylation  | UP        | DOWN       | Methylation decrease (MDE)   | Protein decrease (TMDS) | MDE+TMDS             | TMDS                 | TMDS                 |
        | Hypomethylation  | UP        | UP         | Methylation decrease (MDE)   | None                    | MDE                  | MDE                  | MDE                  |
        | Hypomethylation  | DOWN      | UP         | mRNA decrease (TPDS)         | Protein increase (TMDE) | TPDS+TMDE            | TPDS+TMDE            | TMDE                 |
        | Hypomethylation  | No Change | UP         | mRNA decrease (TPDS)         | Protein increase (TMDE) | TPDS+TMDE            | TMDE                 | TMDE                 |
        | Hypomethylation  | No Change | DOWN       | mRNA decrease (TPDS)         | Protein decrease (TMDS) | TPDS+TMDS            | TMDS                 | TMDS                 |
        | Hypomethylation  | UP        | No Change  | Methylation decrease (MDE)   | Protein decrease (TMDS) | MDE+TMDS             | MDE+TMDS             | TMDS                 |
        | Hypomethylation  | DOWN      | No Change  | mRNA decrease (TPDS)         | Protein increase (TMDE) | TPDS+TMDE            | TPDS+TMDE            | TMDE                 |
        | Hypomethylation  | No Change | No Change  | Methylation decrease (ncRNA) | None                    | MDE+ncRNA            | MDE_ncRNA            | MDE_ncRNA            |
        | No Change        | DOWN      | UP         | mRNA decrease (TPDS)         | Protein increase (TMDE) | TPDS+TMDE            | TPDS+TMDE            | TMDE                 |
        | No Change        | UP        | DOWN       | mRNA increase (TPDE)         | Protein decrease (TMDS) | TPDE+TMDS            | TPDE+TMDS            | TMDS                 |
        | No Change        | DOWN      | DOWN       | mRNA decrease (TPDS)         | None                    | TPDS                 | TPDS                 | TPDS                 |
        | No Change        | UP        | UP         | mRNA increase (TPDE)         | None                    | TPDE                 | TPDE                 | TPDE                 |
        | No Change        | No Change | UP         | Protein increase (TMDE)      | None                    | TMDE                 | TMDE                 | TMDE                 |
        | No Change        | No Change | DOWN       | Protein decrease (TMDS)      | None                    | TMDS                 | TMDS                 | TMDS                 |
        | No Change        | UP        | No Change  | mRNA increase (TPDE)         | Protein decrease (TMDS) | TPDE+TMDS            | TPDE+TMDS            | TMDS                 |
        | No Change        | DOWN      | No Change  | mRNA decrease (TPDS)         | Protein increase (TMDE) | TPDS+TMDE            | TPDS+TMDE            | TMDE                 |
        | No Change        | No Change | No Change  | NoChange                     | NoChange                | NoChange             | NoChange             | NoChange             |
    """


Please post questions and issues related to sci-rcm on the `Issues <https://github.com/ArianeMora/scircm/issues>`_  section of the GitHub repository.

.. toctree::
   :caption: Reproducibility
   :maxdepth: 1

   examples/R4_Part1_RCM_Figure1
   examples/R4_Part1_RCM_Figure2
   examples/R4_Part1_RCM_Figure3
   examples/R4_Part1_RCM_Figure4
   examples/R4_Part1_RCM_Figure5-Late
   examples/R4_Part1_RCM_Figure5
   examples/R4_Part1_RCM_S.Fig1
   examples/R5_PanCan_Part1_RCM_Figure1
   examples/R5_PanCan_Part1_RCM_Figure2
   examples/R5_PanCan_Part1_RCM_Figure3
   examples/R5_PanCan_Part1_RCM_Figure4
   examples/R5_PanCan_Part1_RCM_S.Fig1
   examples/R6_ccRCC_Part2_VAE_Figure4
   examples/R6_ccRCC_Part2_VAE_Figure5
   examples/R6_ccRCC_Part2_VAE_Figure6
   examples/R7_PanCan_Part2_VAE_Figure4
   examples/R7_PanCan_Part2_VAE_Figure5
   examples/R7_PanCan_Part2_VAE_Figure6
   examples/R8_Comparison_Figure1
   examples/R8_Part4_ITH_Processing-Larger
   examples/R8_Part4_ITH_S.Fig1
   examples/R8_Part4_ITH_protein_analysis
   examples/RCM_all_patients-GENES_Rows
   examples/Benchmarking
   examples/N01_RNAProcessing
   examples/N02_ClinicalProcessing
   examples/N03_protein_imputation
   examples/N04_ProteinProcessing
   examples/N05_MethProcessing
   examples/N06_PhosphoProteinProcessing
   examples/N07_RNAProcessing
   examples/N08_DatasetGeneration
   examples/N09_TvN_DE_DA_CpG
   examples/N10_MethylationFilter
   examples/N11_SiRCle_RCM
   examples/N12_SiRCle_ORA_GSEA
   examples/N13_SiRCle_ORA_vis
   examples/N14_TF
   examples/N15_VAE_Integration
   examples/N16_Integration_GSEA
   examples/N17_SingleCell_setup
   examples/N17_singlecell_StageIV-StageI
   examples/N18_singlecell_PBRM1-BAP1
   examples/N20_Metabolomics_TvN
   examples/N21_Metabolomics_RCM
   examples/N22_Metabolomics_VAE