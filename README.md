# Identifying Biomarkers for Drug Response in Canine Osteosarcoma

## Project Description

This project aims to identify potential biomarkers for drug response in canine osteosarcoma. By analyzing RNA-sequencing data from a panel of 13 canine osteosarcoma cell lines and their corresponding drug sensitivity profiles, we sought to uncover gene expression and pathway activity patterns associated with differential drug responses. Specifically, we focused on three drugs: Olaparib, Carboplatin, and Alisertib.

Through correlation analysis and statistical testing, we identified several genes and pathways that correlated with drug sensitivity. For example, we found that lower expression of genes like TMEM131L, APC, and DNAJB1 was associated with increased sensitivity to Olaparib. Additionally, pathway enrichment analysis revealed that higher Interferon Alpha Response scores were linked to Olaparib sensitivity. Similar findings were observed for Carboplatin and Alisertib, with different genes and pathways emerging as potential biomarkers.

The scripts and data files provided in this repository detail the computational methods used to conduct the analysis, including data preprocessing, differential expression analysis, pathway enrichment analysis, and correlation analysis.

## Purpose 

Identify potential new biomarkers for drug response in canine osteosarcoma to aid in the discovery of new treatments for both canine and human patients.


## Key Findings

* Olaparib: Lower expression of genes like TMEM131L, APC, and DNAJB1 was associated with increased sensitivity. Higher Interferon Alpha Response scores were linked to Olaparib sensitivity.
* Carboplatin: Sensitivity correlated with lower expression of APC and DNAJB1. Pathway enrichment analysis revealed a link to the Hallmark-Apical Junction pathway.
* Alisertib: Sensitivity correlated with LARP4B expression, while resistance was associated with higher Bile Acid Metabolism pathway scores.

## Contents
* Compare Normalization Methods
  * Scripts
  * Output
* GSVA Analysis
  * Scripts
  * Output
* Limma Analysis
  * Scripts
* Useful Scripts

## Packages Used 
### All scripts have needed packages listed at the top for easy install
* GSVA 1.52.3
* Limma 3.60.4
* GSEABase 1.66.0
* Flextable 0.9.6
* Gprofiler 0.2.3
* Mygene 1.40.0
* ComplexHeatmap 2.20.0
* Ggpubr 0.6.0.999
* Corrplot 0.94
* EDASeq 2.38.0
* Tidyverse 2.0.0
  * Tidyr 1.3.1
  * Ggplot2 3.5.1
  * Dplyr 1.1.4
* Conflicted 1.2.0
* pheatmap 1.0.12
* matrixStats 1.4.1

