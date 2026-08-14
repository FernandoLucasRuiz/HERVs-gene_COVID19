# HERV activation in COVID-19 lungs

This repository contains the analysis code associated with the manuscript:

**Severe pulmonary COVID-19 is marked by broad Human endogenous retrovirusesactivation and inflammatory tissue remodeling**

Fernando Lucas-Ruiz, Azahara Maria Garcia Serna, Clara Salas, Angel Esteban, Kevin De Man, Yann Heylen, Wim Vanden Berghe, Alberto Baroja-Mazo, Pablo Pelegrin and Santiago Cuevas.

## Overview

This study investigates endogenous retroviral activation in lung parenchyma from patients with fatal COVID-19 pneumonia and control lung tissue. The analyses integrate host-gene transcriptomics, locus-specific HERV quantification, cell-state deconvolution, immune-cell pseudobulk analysis, HERV-gene co-expression network analysis and exploratory DNA methylation profiling.

The main analyses include:

* host-gene differential expression from bulk RNA-seq;
* locus-specific HERV quantification using Telescope;
* functional enrichment of host genes and HERV-associated genes;
* cell-state deconvolution using a Human Lung Cell Atlas-derived CIBERSORTx reference;
* donor-level immune-cell pseudobulk analysis from the Human Lung Cell Atlas;
* integrated HERV-gene co-expression network analysis using GWENA;
* exploratory DNA methylation analysis from EPIC v2.0 arrays;

## Data availability

Raw human RNA-seq and DNA methylation array data are not included in this repository. These data will be deposited in an appropriate public or controlled-access repository under accession number `...............`.

## Software requirements

The analyses were performed in R version 4.4.0. The main R packages used include:

* DESeq2
* tximport
* AnnotationDbi
* org.Hs.eg.db
* limma
* edgeR
* clusterProfiler
* minfi
* GWENA
* ggplot2
* dplyr
* tidyverse
* uwot
* circlize
* igraph
* corrplot
* ggpubr
* ggrepel

External tools used include:

* FastQC
* Salmon
* Bowtie2
* Samtools
* Telescope
* CIBERSORTx


## Contact

For questions about the analyses, please contact:

Fernando Lucas-Ruiz
Biomedical Research Institute of Murcia (IMIB-Pascual Parrilla)
Email: fernando.lucas@um.es

Santiago Cuevas
Biomedical Research Institute of Murcia (IMIB-Pascual Parrilla)
Email: santiago.cuevas@imib.es
