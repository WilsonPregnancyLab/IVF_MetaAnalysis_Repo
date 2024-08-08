# IVF_MetaAnalysis_Repo
 All project scripts and plots for IVF meta-analysis 2024 publication

## ExtractingMetadata_and_MakingSeriesDirectories.R 
Extract clinical and technical data for all samples in GEO data sets using GEOquery 

## Filtering_and_Normalization.R 
1. Complete sex prediction comparison quality control
2. load in all datasets and normalize using adjustedFunnorm
3. filtering step -> remove pood performing probes (bad p values, missing B values), SNP probes, crosshybridizing probes, and non-variable placental probes
