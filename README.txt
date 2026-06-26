# Code Repository for: "Ochratoxin A increases the risk of clear cell renal cell carcinoma through HRAS protein: a comprehensive network analysis"

## Overview
This repository contains the full suite of author-generated computational scripts used for the data analysis, machine learning, immune infiltration evaluation, and molecular dynamics preparation in the manuscript. All codes are made publicly available to ensure transparency and reproducibility in accordance with PLOS ONE data-sharing policies.

## System Requirements
* **Computing Environment:** R version 4.5.1
* **Molecular Docking Software:** AutoDock Vina (vina.exe)
* **Main R Packages Required:** clusterProfiler, org.Hs.eg.db, ggplot2, ggvenn, limma, glmnet (for LASSO), e1071 (for SVM-RFE), pROC, CIBERSORT (or equivalent immune deconvolution algorithms), ComplexHeatmap.

## Analytical Pipeline & Script Descriptions

### Module 1: Network Pharmacology & Target Identification
* `01_Target_Database_Union.R`: Merges compound-related targets from ChEMBL, SEA, and Swiss Target Prediction databases.
* `02_OMIM_GeneCards_Union.R`: Merges disease-related targets from OMIM and GeneCards databases.
* `03_First_Venn_Intersection.R`: Identifies overlapping genes between compound targets and disease targets.
* `04_GO_Enrichment_Analysis.R`: Performs Gene Ontology (GO) functional enrichment on the intersected targets.
* `05_KEGG_Pathway_Analysis.R`: Evaluates KEGG pathway enrichment for the intersected targets.

### Module 2: Transcriptomics & Differential Expression (GEO: GSE36376)
* `06_GEO_Data_Processing.R`: Downloads and standardizes raw expression data (GSE36376).
* `07_Add_Sample_Type_Correction.R`: Annotates sample clinical traits and applies batch-effect corrections if necessary.
* `08_Differential_Expression.R`: Identifies Differentially Expressed Genes (DEGs) between control and disease samples.
* `09_Second_Venn_Intersection.R`: Intersects DEGs with the network pharmacology targets to narrow down core genes.
* `10_Extract_Intersection_Matrix.R`: Extracts the expression matrix for these candidate core genes.

### Module 3: Machine Learning Feature Selection
* `11_LASSO_Regression_Analysis.R`: Applies LASSO penalized regression to select diagnostic feature genes.
* `12_SVM_RFE_Analysis.R`: Implements Support Vector Machine-Recursive Feature Elimination (SVM-RFE) for feature selection.
* `13_Third_Venn_Intersection.R`: Intersects the results from LASSO and SVM-RFE to determine the final hub genes.
* `14_GEO_Gene_Validation_ROC.R`: Validates the hub genes' differential expression and plots Receiver Operating Characteristic (ROC) curves.

### Module 4: Immune Infiltration Analysis
* `15_Immune_Infiltration_Analysis.R`: Calculates immune cell fraction scores using the expression matrix.
* `16_Immune_Heatmap.R`: Visualizes immune cell infiltration landscapes via heatmaps.
* `17_Immune_Boxplot.R`: Generates boxplots to compare immune cell abundance between grouped/individual samples.
* `18_Gene_Immune_Correlation.R`: Conducts correlation analysis between the hub genes and immune cell infiltration scores.

### Module 5: Molecular Docking & Dynamics
* `19_AutoDock_Vina_Commands.txt`: Contains detailed command-line instructions (e.g., config data, log mapping) for executing small-molecule/protein molecular docking via AutoDock Vina.
* `20_Merge_Gro_Files.R`: Post-processing script to merge ligand and receptor `.gro` / `.top` files to prepare for subsequent molecular dynamics simulations.