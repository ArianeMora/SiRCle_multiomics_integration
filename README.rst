***********************
**SiRCleMansucript:**
***********************

**********************************************************************************************************************
SiRCle (Signature Regulatory Clustering) model integration reveals mechanisms of phenotype regulation in renal cancer
**********************************************************************************************************************

Code
====
If you're interested in running SiRCLe: please see our R package: https://github.com/ArianeMora/SiRCleR
or python package for installation instructions, examples and documentation: https://github.com/ArianeMora/scircm

Paper
=====
We're in the process of submission but you can read our preprint here: https://www.biorxiv.org/content/10.1101/2022.07.02.498058v1
Always happy to hear feedback and suggestions :)

Information
===========
This site hosts the information associated with the paper: **SiRCle (Signature Regulatory Clustering) model integration reveals mechanisms of phenotype regulation in renal cancer**.
Here we provide the code and data used for all the analyses in the paper and link to the packages we developed as part of
producing the paper. Note the code in `manuscript_reproducibility` is from prior to our revision. While `to_publish_clean` is
 all the notebooks and HTML outputs from every figure in the paper. We also have all data, code and everything on Zenodo: https://zenodo.org/records/14176842

Links to analyses and data
--------------------------

- `Code and tutorials: <https://github.com/ArianeMora/scircm>`_
- `Analyses, reproducibility, and processed data <https://arianemora.github.io/SiRCle_multiomics_integration/>`_
- `Data and everything <https://zenodo.org/records/14176842>`_

Places where this (or a package we developed for this) has been presented
-------------------------------------------------------------------------

.. list-table::
   :widths: 15 30 15
   :header-rows: 1

   * - Date
     - Conference
     - Type
   * - 28 April 2021
     - Melbourne bioinformatics seminar series
     - Presentation
   * - 25 May 2021
     - `Vizbi <https://vizbi.org/Posters/2021/vD02>`_
     - Poster
   * - 15 - 17 Sep 2021
     - `Multiomics to Mechanisms: Challenges in Data Integration <https://www.embl.org/about/info/course-and-conference-office/events/ees21-09/>`_
     - Short talk
   * - 15 Jan 2022
     - Multi-Omics ONLINE - Webinar 2: Data integration and interpretation to unveil novel insights
     - Talk

Authors
=======

Ariane Mora^1, Christina Schmidt^2,3, Brad Balderson1, Christian Frezza3#, Mikael Bodén1#

1) School of Chemistry and Molecular Biosciences, University of Queensland, Molecular Biosciences Building 76, St Lucia QLD 4072, Australia.
2) Medical Research Council Cancer Unit, University of Cambridge, Hutchison/MRC Research Centre, Box 197, Cambridge Biomedical Campus, Cambridge CB2 0X2, United Kingdom
3) CECAD Research Center, University Hospital Cologne, Joseph-Stelzmann-Str. 26, 50931 Cologne, Germany

^Joint first authors equally contributed; the order is interchangeable and up to the authors discretion
#Joint last authors

Note Christina and Ariane are equal joint first authors and the authors may swap the order of their names as they so choose :)

Abstract
========
Clear cell renal cell carcinoma (ccRCC) tumours develop and progress via complex remodelling of the kidney epigenome, transcriptome, proteome, and metabolome. Given the subsequent tumour and inter-patient heterogeneity, drug-based treatments report limited success, calling for multi-omics studies to extract regulatory relationships, and ultimately, to develop targeted therapies. However, current methods are unable to extract nonlinear multi-omics perturbations.

Here, we present SiRCle (Signature Regulatory Clustering), a novel method to integrate DNA methylation, RNA-seq and proteomics data. Applying SiRCle to a case study of ccRCC, we disentangle the layer (DNA methylation, transcription and/or translation) where dysregulation first occurs and find the primary biological processes altered. Next, we detect regulatory differences between patient subsets by using a variational autoencoder to integrate omics’ data followed by statistical comparisons on the integrated space. In ccRCC patients, SiRCle allows to identify metabolic enzymes and cell-type-specific markers associated with survival along with the likely molecular driver behind the gene’s perturbations.

