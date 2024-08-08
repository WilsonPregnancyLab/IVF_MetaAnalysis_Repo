# IVF_MetaAnalysis_Repo
 All project scripts and plots for IVF meta-analysis 2024 publication - order of read me is the order scipts should be run

## ExtractingMetadata_and_MakingSeriesDirectories.R 
Extract clinical and technical data for all samples in GEO data sets using GEOquery 

## Filtering_and_Normalization.R 
1. Complete sex prediction comparison quality control
2. load in all datasets and normalize using adjustedFunnorm
3. filtering step -> remove pood performing probes (bad p values, missing B values), SNP probes, crosshybridizing probes, and non-variable placental probes

## Methylation_Analyses.R
1. Comparison of DNA methylation in autosomal and X and Y chromosome data between IVF and spontaneous placentae in mixed and fetal sex stratified populations using linear modelling.
2. Plotting of results in volcano plots. 

## Variability_and_Genometracks.R
1. Calculate coefficients of variance for delta B values in IVF male, IVF female, Spontaneous male and Spontaneous female populations and plot densities to compare with KS test.
2. Create genome track figures for significant probe regions

## CellDeconvolution.R
1. Calculate placental cell proportions for samples
2. Compare placental cell proportions using an ANOVA and bonferonni test between IVF male, IVF female, Spontaneous male and Spontaneous female populations
3. Plot bar graph of cell proportions
