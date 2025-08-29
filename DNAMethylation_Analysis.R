#R version 4.4.2
library(minfi) #version 1.50.0
library(limma) #version 3.60.4

#Read in the filtered and normalized methylation data (RGset)
placmet_adjFunnorm_allfiltered <- readRDS ("./placmet_adjFunnorm_allfiltered.rds")
#Call X and Y chromosome probes to allow stratification of dataset into Autosomes, X chromosomes, and Y chromosomes
probeInfo <- as.data.frame(cbind(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations, 
                                 IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other, 
                                 IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Manifest)) 
probeInfo$probeID <- rownames(probeInfo)
chrXprobes <- subset(probeInfo, probeInfo$chr == "chrX") #19090 probes
chrYprobes <- subset(probeInfo, probeInfo$chr == "chrY") #537 probes
#Get Beta values from the RGset file 
placmet_adjFunnorm_filtbetas <- getBeta(placmet_adjFunnorm_allfiltered) 
#Convert Beta values to Mvalues
placmet_adjFunnorm_mvalues <- getM(placmet_adjFunnorm_allfiltered)
placmet_adjFunnorm_mvalues <- as.data.frame(placmet_adjFunnorm_mvalues)
#Autosomes Only 
placmet_adjFunnorm_mvalues_auto <- as.data.frame(
  placmet_adjFunnorm_mvalues[!rownames(placmet_adjFunnorm_mvalues) %in% 
                               c(chrXprobes$probeID, chrYprobes$probeID), ]) #296545 autosomal probes
allDNAm <- t(placmet_adjFunnorm_mvalues_auto)
allDNAm2 <- as.data.frame(apply(allDNAm, 2, as.numeric))
allDNAm2$Sample_Name <- rownames(allDNAm)
#Merge in metadata information for each sample 
pheno <- read.csv("./PlacMet_MetaData_outliersremoved.csv")
allmergedd <- merge(pheno, allDNAm2)
#view this table to make sure it looks okay - Each row = sample, first columns = Metadata, all other columns = each probe 
head(allmergedd[, 1:10], 5) 

#Set function to run the linear model and extract adjusted and unadjusted Pvalues
run_model <- function(data, design_formula, phenotype_cols, robust=FALSE, weights=NULL) {
  # Prepare the design matrix
  design <- model.matrix(design_formula, data = data)
  
  #make numeric data for limma
  Mval_numbers <- as.matrix(apply(data[, phenotype_cols], 2, as.numeric)) 
  
  # Fit the linear model
  if (is.null(weights)){
    fit <- lmFit(t(Mval_numbers), design)
  } else {
    fit <- lmFit(t(Mval_numbers), design, weights=weights)
  }
  fit <- eBayes(fit, robust = robust)
  
  # Extract Raw and Adjusted P-Values from the linear model outputs
  table <- topTable(fit, coef = "GroupSpontaneous", number = Inf, adjust = "fdr", sort.by = "none")
  raw_p <- table[, "P.Value"]
  adj_p <- table[, "adj.P.Val"]
  
  return(data.frame(Raw_P = raw_p, Adj_P = adj_p))
}
#set which columns are probes 
phenotype_cols <- 5:296549  # Replace with the actual column indices for CpG sites

# Run autosome models for the whole sample
whole_adj_p <- run_model(allmergedd, ~ Group + Fetal_Sex + Study, phenotype_cols)
# Run autosome models for males only
male_adj_p <- run_model(allmergedd[allmergedd$Fetal_Sex == "Male", ], ~ Group + Study, phenotype_cols)
# Run autosome models for females only
female_adj_p <- run_model(allmergedd[allmergedd$Fetal_Sex == "Female", ], ~ Group + Study, phenotype_cols)
# Combine results into a data frame
autosome_results <- data.frame(
  CpG = colnames(allmergedd[, phenotype_cols]),
  Whole_Raw_P = whole_adj_p$Raw_P,
  Whole_Adj_P = whole_adj_p$Adj_P,
  Male_Raw_P = male_adj_p$Raw_P,
  Male_Adj_P = male_adj_p$Adj_P,
  Female_Raw_P = female_adj_p$Raw_P,
  Female_Adj_P = female_adj_p$Adj_P
)
#checking how many significant hits found 
sum(autosome_results$Whole_Adj_P < 0.05, na.rm = TRUE) #127 CpGs
sum(autosome_results$Male_Adj_P < 0.05, na.rm = TRUE) #9 CpGs
sum(autosome_results$Female_Adj_P < 0.05, na.rm = TRUE) #0 CpGs

write.csv(autosome_results, "./autosome_results.csv")

#Running this Analysis in X chromosomes 
placmet_adjFunnorm_mvalues_X <- placmet_adjFunnorm_mvalues[rownames(placmet_adjFunnorm_mvalues) %in% chrXprobes$probeID,] #8779 X probes
XMval <- t(placmet_adjFunnorm_mvalues_X)
XMval2 <- as.data.frame(apply(XMval, 2, as.numeric))
XMval2$Sample_Name <- rownames(XMval)
allmergeddX <- merge(pheno, XMval2)
head(allmergeddX[, 1:10], 5) 
#set which columns are probes 
phenotype_cols <- 5:8783
#Run X chromosome models for males only 
maleX_adj_p <- run_model(allmergeddX[allmergeddX$Fetal_Sex == "Male", ], ~ Group + Study, phenotype_cols)
#Run X chromosome models for females only
femaleX_adj_p <- run_model(allmergeddX[allmergeddX$Fetal_Sex == "Female", ], ~ Group + Study, phenotype_cols)
#combine results into a data frame 
XChromosome_results <- data.frame(
  CpG = colnames(allmergeddX[, phenotype_cols]),
  MaleX_Raw_P = maleX_adj_p$Raw_P,
  MaleX_Adj_P = maleX_adj_p$Adj_P,
  FemaleX_Raw_P = femaleX_adj_p$Raw_P,
  FemaleX_Adj_P = femaleX_adj_p$Adj_P
)
#checking how many significant hits found
sum(XChromosome_results$MaleX_Adj_P < 0.05, na.rm = TRUE) #1 CpG
sum(XChromosome_results$FemaleX_Adj_P < 0.05, na.rm = TRUE) #0 CpGs

write.csv(XChromosome_results, "./XChromosome_results.csv")


#running this Analysis in Y Chomosome - males only 
placmet_adjFunnorm_mvalues_Y <- placmet_adjFunnorm_mvalues[rownames(placmet_adjFunnorm_mvalues) %in% chrYprobes$probeID,] #227 Y probes
YMval <- t(placmet_adjFunnorm_mvalues_Y)
YMval2 <- as.data.frame(apply(YMval, 2, as.numeric))
YMval2$Sample_Name <- rownames(YMval)
allmergeddY <- merge(pheno, YMval2)
head(allmergeddY[, 1:10], 5) 
#set which columns are probes 
phenotype_cols <- 5:231 
# Run Y Chromosome model for males only
maleY_adj_p <- run_model(allmergeddY[allmergeddY$Fetal_Sex == "Male", ], ~ Group + Study, phenotype_cols)
# Combine results into a data frame
YChromosome_results <- data.frame(
  CpG = colnames(allmergeddY[, phenotype_cols]),
  MaleY_Raw_P = maleY_adj_p$Raw_P,
  MaleY_Adj_P = maleY_adj_p$Adj_P
)
#checking how many significant hits found 
sum(YChromosome_results$MaleY_Adj_P < 0.05, na.rm = TRUE) #0 CpGs

write.csv(YChromosome_results, "./YChromosome_results.csv")


significant_X_results <- XChromosome_results_allinfo[
    XChromosome_results_allinfo$MaleX_Adj_P < 0.05 |
    XChromosome_results_allinfo$FemaleX_Adj_P < 0.05,
]
