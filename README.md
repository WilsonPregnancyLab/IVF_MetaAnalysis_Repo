# IVF_MetaAnalysis_Repo
 All project scripts and plots for IVF meta-analysis 2025 publication - order of read me is the order scipts should be run

## ExtractingMetadata_and_MakingSeriesDirectories.R 
Extract clinical and technical data for all samples in GEO data sets using GEOquery 

## Filtering_and_Normalization.R 
1. Complete sex prediction comparison quality control
2. load in all datasets and normalize using adjustedFunnorm
3. filtering step -> remove poor performing probes (bad p values, missing B values), SNP probes, crosshybridizing probes, and non-variable placental probes

## DNAMethylation_Analysis.R
1. Load in beta values for samples and separate probes by autosomes, X chromosomes, and Y chromosomes
2. Run linear models on M values - comparing IVF to spontnaoues in mixed fetal sex, male, and female strata 

## DNAMethylation_Plotting.R
1. Creates all necessary tables for calculating delta betas and formatting to include thresholds of biological significance cutoffs
2. Plots all CpGs in volcano plots for Autosome, X chromosome, and Y chromosome linear model results

## Linear_Model_Validation.R
1. Runs models outline in the sensitivity analysis: Excluding study with greatest imbalance, Robust = TRUE
2. Build table with top 15 significant CpGs of the linear model sensitivity analysis + original model 

## Cell_Deconvolution.R
1. Runs placental cell deconvolution for samples stratified by conception type and fetal sex and gives proportions of trophoblasts, syncytiotrophoblasts, stromal cells, endothelial cells, Hoffbauer cells, and nRBCs

## Comparing_CellType.R
1. Runs separate models to compare fetal sex in stratified conception types, and conception type stratified by fetal sex 
2. Runs a mixed effects model for both conception type and fetal sex 
3. Plots proportions for each cell type across conception type and fetal sex strata
