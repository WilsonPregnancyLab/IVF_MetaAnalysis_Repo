#R Version - R4.4.1
#Packages
library(limma) #version 3.60.4 
library(lumi) #version 2.56.0 
library(stringr) #version 1.5.1
library(ggplot2) #version 3.5.1
library(gridExtra) #version 2.3
library(ggrepel) # version 0.9.5

#Creating groups for sex stratification
metadata <- read.csv("/workspace/lab/wilsonslab/lemairem/Placenta_Based_Studies/Met_array/rerun_June2024_outlier_exclusion/PlacMet_MetaData_outliersremoved.csv")
metadata$Sentrix_ID <- str_extract(metadata$Sample_Name, "(?<=_)[0-9]+(?=_)")
metadata$Sentrix_Position <- str_extract(metadata$Sample_Name, "R[0-9]{2}")
metadata$Group <- as.factor(metadata$Group)
metadata$Fetal_Sex <- as.factor(metadata$Fetal_Sex)
metadata$Sentrix_ID <- as.factor(metadata$Sentrix_ID)
metadata$Sentrix_Position <- as.factor(metadata$Sentrix_Position)
metadata$Study <- as.factor(metadata$Study)

males <- subset(metadata, metadata$Fetal_Sex == "Male")
males$Group <- as.factor(males$Group)
males$Fetal_Sex <- as.factor(males$Fetal_Sex)
males$Sentrix_ID <- as.factor(males$Sentrix_ID)
males$Sentrix_ID <- droplevels(males$Sentrix_ID)
males$Sentrix_Position <- as.factor(males$Sentrix_Position)
males$Study <- as.factor(males$Study)

females <- subset(metadata, metadata$Fetal_Sex == "Female")
females$Group <- as.factor(females$Group)
females$Fetal_Sex <- as.factor(females$Fetal_Sex)
females$Sentrix_ID <- as.factor(females$Sentrix_ID)
females$Sentrix_ID <- droplevels(females$Sentrix_ID) 
females$Sentrix_Position <- as.factor(females$Sentrix_Position)
females$Study <- as.factor(females$Study)

placmet_adjFunnorm_filtbetas <- getBeta(placmet_adjFunnorm_allfiltered)
placmet_filtfun_F <- placmet_adjFunnorm_filtbetas[, females$Sample_Name]
# 281 females - gives all the GSMs in either IVF or spont that are female
placmet_filtfun_M <- placmet_adjFunnorm_filtbetas[, males$Sample_Name]
# 294 males - gives all the GSMs in either IVF or spont that are male 

#Calculating the average Delta Beta values between control and IVF 
Spontaneousmetadata_all <- subset(metadata, metadata$Group == "Spontaneous")
Spontaneousmetadata_F <- subset(females, females$Group == "Spontaneous")
Spontaneousmetadata_M <- subset(males, males$Group == "Spontaneous")
# 479 spontaneous all (both males and females)
# 241 spontaneous females 
# 238 spontaneous males
IVFmetadata_all <- subset(metadata, metadata$Group == "IVF")
IVFmetadata_F <- subset(females, females$Group == "IVF")
IVFmetadata_M <- subset(males, males$Group == "IVF")
# 96 IVF all
# 40 IVF females 
# 56 IVF males

#Probe annotations 
probeInfo <- as.data.frame(cbind(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations, 
                                 IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other, 
                                 IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Manifest)) 
probeInfo$probeID <- rownames(probeInfo)
chrXprobes <- subset(probeInfo, probeInfo$chr == "chrX") #19090 probes
chrYprobes <- subset(probeInfo, probeInfo$chr == "chrY") #537 probes

# Total population Betas table 
# Autosomal Only Probes included for Beta Table
placmet_IVF_all_autosomes <-  as.data.frame(placmet_adjFunnorm_filtbetas[!rownames(placmet_adjFunnorm_filtbetas) %in% c(chrXprobes$probeID,chrYprobes$probeID), IVFmetadata_all$Sample_Name]) #dim: 296545, 96 
placmet_SPONT_all_autosomes <- as.data.frame(placmet_adjFunnorm_filtbetas[!rownames(placmet_adjFunnorm_filtbetas) %in% c(chrXprobes$probeID,chrYprobes$probeID), Spontaneousmetadata_all$Sample_Name]) #dim: 296545 479
# calculate All average betas
placmet_SPONT_all_autosomes$AvgBSPONT <- rowMeans(placmet_SPONT_all_autosomes)
placmet_SPONT_all_autosomes$ProbeSPONT <- rownames(placmet_SPONT_all_autosomes)
placmet_IVF_all_autosomes$AvgBIVF <- rowMeans(placmet_IVF_all_autosomes)
placmet_IVF_all_autosomes$ProbeIVF <- rownames(placmet_IVF_all_autosomes)
# Merge All Table 
placmet_AllAvgbetas_autosomes <- merge(placmet_SPONT_all_autosomes[,c("AvgBSPONT","ProbeSPONT")], placmet_IVF_all_autosomes[,c("AvgBIVF","ProbeIVF")], by = "row.names")
placmet_AllAvgbetas_autosomes$deltaB <- placmet_AllAvgbetas_autosomes$AvgBSPONT - placmet_AllAvgbetas_autosomes$AvgBIVF 
rownames(placmet_AllAvgbetas_autosomes) <- placmet_AllAvgbetas_autosomes$ProbeSPONT
write.table(placmet_AllAvgbetas_autosomes, sep = "\t", file = "./placmet_AllAvgbetas_autosomes.tsv")

# Male population Betas table 
#Autosomes
placmet_IVF_M_autosomes <- as.data.frame(placmet_filtfun_M[!rownames(placmet_filtfun_M) %in% c(chrXprobes$probeID,chrYprobes$probeID), IVFmetadata_M$Sample_Name]) #dim 296545, 56
placmet_SPONT_M_autosomes <- as.data.frame(placmet_filtfun_M[!rownames(placmet_filtfun_M) %in% c(chrXprobes$probeID,chrYprobes$probeID), Spontaneousmetadata_M$Sample_Name]) #dim 296545, 238
# calculate male autosome average betas
placmet_SPONT_M_autosomes$AvgBSPONT <- rowMeans(placmet_SPONT_M_autosomes)
placmet_SPONT_M_autosomes$ProbeSPONT <- rownames(placmet_SPONT_M_autosomes)
placmet_IVF_M_autosomes$AvgBIVF <- rowMeans(placmet_IVF_M_autosomes)
placmet_IVF_M_autosomes$ProbeIVF <- rownames(placmet_IVF_M_autosomes)
# Merge male autosomes Table 
placmet_MaleAvgbetas_autosomes <- merge(placmet_SPONT_M_autosomes[,c("AvgBSPONT","ProbeSPONT")], placmet_IVF_M_autosomes[,c("AvgBIVF","ProbeIVF")], by = "row.names")
placmet_MaleAvgbetas_autosomes$deltaB <- placmet_MaleAvgbetas_autosomes$AvgBSPONT - placmet_MaleAvgbetas_autosomes$AvgBIVF
rownames(placmet_MaleAvgbetas_autosomes) <- placmet_MaleAvgbetas_autosomes$ProbeSPONT
write.table(placmet_MaleAvgbetas_autosomes, sep = "\t", file = "./placmet_MaleAvgbetas_autosomes.tsv")
#Male X Chromosome 
placmet_IVF_M_X <- as.data.frame(placmet_filtfun_M[rownames(placmet_filtfun_M) %in% chrXprobes$probeID, IVFmetadata_M$Sample_Name]) #dim: 8779, 56
placmet_SPONT_M_X <- as.data.frame(placmet_filtfun_M[rownames(placmet_filtfun_M) %in% chrXprobes$probeID, Spontaneousmetadata_M$Sample_Name]) #dim: 8779, 238
# calculate male X average betas
placmet_SPONT_M_X$AvgBSPONT <- rowMeans(placmet_SPONT_M_X)
placmet_SPONT_M_X$ProbeSPONT <- rownames(placmet_SPONT_M_X)
placmet_IVF_M_X$AvgBIVF <- rowMeans(placmet_IVF_M_X)
placmet_IVF_M_X$ProbeIVF <- rownames(placmet_IVF_M_X)
# Merge male X table
placmet_MaleAvgbetas_X <- merge(placmet_SPONT_M_X[,c("AvgBSPONT","ProbeSPONT")], placmet_IVF_M_X[,c("AvgBIVF","ProbeIVF")], by = "row.names")
placmet_MaleAvgbetas_X$deltaB <- placmet_MaleAvgbetas_X$AvgBSPONT - placmet_MaleAvgbetas_X$AvgBIVF
rownames(placmet_MaleAvgbetas_X) <- placmet_MaleAvgbetas_X$ProbeSPONT
write.table(placmet_MaleAvgbetas_X, sep = "\t", file = "./placmet_MaleAvgbetas_X.tsv")
#Male Y Chromosome
placmet_IVF_M_Y <- as.data.frame(placmet_filtfun_M[rownames(placmet_filtfun_M) %in% chrYprobes$probeID, IVFmetadata_M$Sample_Name]) #dim: 227, 56
placmet_SPONT_M_Y <- as.data.frame(placmet_filtfun_M[rownames(placmet_filtfun_M) %in% chrYprobes$probeID, Spontaneousmetadata_M$Sample_Name]) #dim: 227, 238
# calculate male Y average betas
placmet_SPONT_M_Y$AvgBSPONT <- rowMeans(placmet_SPONT_M_Y)
placmet_SPONT_M_Y$ProbeSPONT <- rownames(placmet_SPONT_M_Y)
placmet_IVF_M_Y$AvgBIVF <- rowMeans(placmet_IVF_M_Y)
placmet_IVF_M_Y$ProbeIVF <- rownames(placmet_IVF_M_Y)
# Merge male Y table
placmet_MaleAvgbetas_Y<- merge(placmet_SPONT_M_Y[,c("AvgBSPONT","ProbeSPONT")], placmet_IVF_M_Y[,c("AvgBIVF","ProbeIVF")], by = "row.names")
placmet_MaleAvgbetas_Y$deltaB <- placmet_MaleAvgbetas_Y$AvgBSPONT - placmet_MaleAvgbetas_Y$AvgBIVF
rownames(placmet_MaleAvgbetas_Y) <- placmet_MaleAvgbetas_Y$ProbeSPONT
write.table(placmet_MaleAvgbetas_Y, sep = "\t", file = "./placmet_MaleAvgbetas_Y.tsv")

# Female population Betas table 
#Female Autosomes
placmet_IVF_F_autosomes <- as.data.frame(placmet_filtfun_F[!rownames(placmet_filtfun_F) %in% c(chrXprobes$probeID,chrYprobes$probeID), IVFmetadata_F$Sample_Name]) #dim 296545, 40
placmet_SPONT_F_autosomes <- as.data.frame(placmet_filtfun_F[!rownames(placmet_filtfun_F) %in% c(chrXprobes$probeID,chrYprobes$probeID), Spontaneousmetadata_F$Sample_Name]) #dim 296545, 241
# calculate female autosome average betas
placmet_SPONT_F_autosomes$AvgBSPONT <- rowMeans(placmet_SPONT_F_autosomes)
placmet_SPONT_F_autosomes$ProbeSPONT <- rownames(placmet_SPONT_F_autosomes)
placmet_IVF_F_autosomes$AvgBIVF <- rowMeans(placmet_IVF_F_autosomes)
placmet_IVF_F_autosomes$ProbeIVF <- rownames(placmet_IVF_F_autosomes)
# Merge female autosomes Table 
placmet_FemaleAvgbetas_autosomes <- merge(placmet_SPONT_F_autosomes[,c("AvgBSPONT","ProbeSPONT")], placmet_IVF_F_autosomes[,c("AvgBIVF","ProbeIVF")], by = "row.names")
placmet_FemaleAvgbetas_autosomes$deltaB <- placmet_FemaleAvgbetas_autosomes$AvgBSPONT - placmet_FemaleAvgbetas_autosomes$AvgBIVF
rownames(placmet_FemaleAvgbetas_autosomes) <- placmet_FemaleAvgbetas_autosomes$ProbeSPONT
write.table(placmet_FemaleAvgbetas_autosomes, sep = "\t", file = "./placmet_FemaleAvgbetas_autosomes.tsv")
#Female X Chromosome 
placmet_IVF_F_X <- as.data.frame(placmet_filtfun_F[rownames(placmet_filtfun_F) %in% chrXprobes$probeID, IVFmetadata_F$Sample_Name]) #dim: 8779, 40
placmet_SPONT_F_X <- as.data.frame(placmet_filtfun_F[rownames(placmet_filtfun_F) %in% chrXprobes$probeID, Spontaneousmetadata_F$Sample_Name]) #dim: 8779, 241
# calculate female X average betas
placmet_SPONT_F_X$AvgBSPONT <- rowMeans(placmet_SPONT_F_X)
placmet_SPONT_F_X$ProbeSPONT <- rownames(placmet_SPONT_F_X)
placmet_IVF_F_X$AvgBIVF <- rowMeans(placmet_IVF_F_X)
placmet_IVF_F_X$ProbeIVF <- rownames(placmet_IVF_F_X)
# Merge female X table
placmet_FemaleAvgbetas_X <- merge(placmet_SPONT_F_X[,c("AvgBSPONT","ProbeSPONT")], placmet_IVF_F_X[,c("AvgBIVF","ProbeIVF")], by = "row.names")
placmet_FemaleAvgbetas_X$deltaB <- placmet_FemaleAvgbetas_X$AvgBSPONT - placmet_FemaleAvgbetas_X$AvgBIVF
rownames(placmet_FemaleAvgbetas_X) <- placmet_FemaleAvgbetas_X$ProbeSPONT
write.table(placmet_FemaleAvgbetas_X, sep = "\t", file = "./placmet_FemaleAvgbetas_X.tsv")

#m value inputs for linear models 
#Whole population 
placmet_adjFunnorm_MVal <- beta2m(placmet_adjFunnorm_filtbetas)
write.table(placmet_adjFunnorm_MVal, sep = "\t", file = "./placmet_adjFunnorm_MVal.tsv")
#whole population Mvalues - Autosomes
Mval_whole_auto <- placmet_adjFunnorm_MVal[!rownames(placmet_adjFunnorm_MVal) %in% c(chrXprobes$probeID,chrYprobes$probeID),]

#Male Population Mvalues 
Mval_male <- beta2m(placmet_filtfun_M)
write.table(Mval_male, sep = "\t", file = "./MVal_male.tsv")
#male autosomes
Mval_male_auto <- Mval_male[!rownames(Mval_male) %in% c(chrXprobes$probeID,chrYprobes$probeID),]
#male XChr
Mval_male_X <- Mval_male[rownames(Mval_male) %in% chrXprobes$probeID,]
#male Ychr
Mval_male_Y <- Mval_male[rownames(Mval_male) %in% chrYprobes$probeID,]

#Female Population MVals 
Mval_female <- beta2m(placmet_filtfun_F)
write.table(Mval_female, sep = "\t", file = "./MVal_female.tsv")
#female autosomes 
Mval_female_auto <- Mval_female[!rownames(Mval_female) %in% c(chrXprobes$probeID,chrYprobes$probeID),]
#female XChr 
Mval_female_X <- Mval_female[rownames(Mval_female) %in% chrXprobes$probeID,]


#LINEAR MODELING 
#whole population autosomes 
SPONTvsIVF_wholemodel_auto <- model.matrix(~ Group + Fetal_Sex + Study, data = metadata)
SPONTvsIVF_wholefit_auto <- lmFit(Mval_whole_auto, SPONTvsIVF_wholemodel_auto)
SPONTvsIVF_wholefit_auto <- eBayes(SPONTvsIVF_wholefit_auto)
tt_SPONTvsIVF_whole_auto <- topTable(SPONTvsIVF_wholefit_auto, n = Inf, adjust = "fdr", coef = "GroupSpontaneous")
print(sum(tt_SPONTvsIVF_whole_auto$adj.P.Val < 0.05))
names(price_anno)[which(names(price_anno) == "SPOT_ID")] <- "probes"
tt_SPONTvsIVF_whole_auto$probes <- rownames(tt_SPONTvsIVF_whole_auto)
names(placmet_AllAvgbetas_autosomes)[which(names(placmet_AllAvgbetas_autosomes) == "ProbeSPONT")] <- "probes"
tt_SPONTvsIVF_whole_auto$probes <- as.factor(tt_SPONTvsIVF_whole_auto$probes)
placmet_AllAvgbetas_autosomes$probes <- as.factor(placmet_AllAvgbetas_autosomes$probes)
placmet_wholepop_auto <- merge(merge(tt_SPONTvsIVF_whole_auto, placmet_AllAvgbetas_autosomes[, c("deltaB","probes")], by = "probes"),
                               price_anno[,c("Closest_TSS_gene_name", "probes")], by = "probes")
write.table(tt_SPONTvsIVF_whole_auto, sep = "\t", file = "./tt_SPONTvsIVF_whole_auto_study.tsv")
write.table(placmet_wholepop_auto, sep = "\t", file = "./placmet_wholepop_auto_study.tsv")

#male autosomes
SPONTvsIVF_modelM <- model.matrix(~ Group + Study, data = males)
SPONTvsIVF_fitM_auto <- lmFit(Mval_male_auto, SPONTvsIVF_modelM_auto) 
SPONTvsIVF_fitM_auto <- eBayes(SPONTvsIVF_fitM_auto)
tt_SPONTvsIVF_M_all_auto <- topTable(SPONTvsIVF_fitM_auto, n = Inf, adjust = "fdr", coef = "GroupSpontaneous")
print(sum(tt_SPONTvsIVF_M_all_auto$adj.P.Val < 0.05)) #9
tt_SPONTvsIVF_M_all_auto$probes <- rownames(tt_SPONTvsIVF_M_all_auto)
names(placmet_MaleAvgbetas_autosomes)[which(names(placmet_MaleAvgbetas_autosomes) == "ProbeSPONT")] <- "probes"
tt_SPONTvsIVF_M_all_auto$probes <- as.factor(tt_SPONTvsIVF_M_all_auto$probes)
placmet_MaleAvgbetas_autosomes$probes <- as.factor(placmet_MaleAvgbetas_autosomes$probes)
placmet_M_fulldata_auto <- merge(merge(tt_SPONTvsIVF_M_all_auto, placmet_MaleAvgbetas_autosomes[, c("deltaB","probes")], by = "probes"),
                                 price_anno[,c("Closest_TSS_gene_name", "probes")], by = "probes")
write.table(tt_SPONTvsIVF_M_all_auto, sep = "\t", file = "./tt_SPONTvsIVF_M_all_auto_study.tsv")
write.table(placmet_M_fulldata_auto, sep = "\t", file = "./placmet_M_fulldata_auto_study.tsv")
maleauto_sig <- subset(placmet_M_fulldata_auto[placmet_M_fulldata_auto$adj.P.Val <0.05,])
write.csv(maleauto_sig, file = "./Sig_Male_autosomes_IVFvsSPONT_adjFunnorm.csv")
#male X 
SPONTvsIVF_fitM_X <- lmFit(Mval_male_X, SPONTvsIVF_modelM) 
SPONTvsIVF_fitM_X <- eBayes(SPONTvsIVF_fitM_X)
tt_SPONTvsIVF_M_X <- topTable(SPONTvsIVF_fitM_X, n = Inf, adjust = "fdr", coef = "GroupSpontaneous")
print(sum(tt_SPONTvsIVF_M_X$adj.P.Val < 0.05)) #1
tt_SPONTvsIVF_M_X$probes <- rownames(tt_SPONTvsIVF_M_X)
names(placmet_MaleAvgbetas_X)[which(names(placmet_MaleAvgbetas_X) == "ProbeSPONT")] <- "probes"
tt_SPONTvsIVF_M_X$probes <- as.factor(tt_SPONTvsIVF_M_X$probes)
placmet_MaleAvgbetas_X$probes <- as.factor(placmet_MaleAvgbetas_X$probes)
placmet_M_fulldata_X <- merge(merge(tt_SPONTvsIVF_M_X, placmet_MaleAvgbetas_X[, c("deltaB","probes")], by = "probes"),
                              price_anno[,c("Closest_TSS_gene_name", "probes")], by = "probes")
write.table(tt_SPONTvsIVF_M_X, sep = "\t", file = "./tt_SPONTvsIVF_M_X_study.tsv")
write.table(placmet_M_fulldata_X, sep = "\t", file = "./placmet_M_fulldata_X_study.tsv")
maleX_sig <- subset(placmet_M_fulldata_X[placmet_M_fulldata_X$adj.P.Val <0.05,])
write.csv(maleX_sig, file = "./Sig_Male_X_IVFvsSPONT_adjFunnorm.csv")
#male Y
SPONTvsIVF_fitM_Y <- lmFit(Mval_male_Y, SPONTvsIVF_modelM) 
SPONTvsIVF_fitM_Y <- eBayes(SPONTvsIVF_fitM_Y)
tt_SPONTvsIVF_M_Y <- topTable(SPONTvsIVF_fitM_Y, n = Inf, adjust = "fdr", coef = "GroupSpontaneous")
print(sum(tt_SPONTvsIVF_M_Y$adj.P.Val < 0.05)) #0 
tt_SPONTvsIVF_M_Y$probes <- rownames(tt_SPONTvsIVF_M_Y)
names(placmet_MaleAvgbetas_Y)[which(names(placmet_MaleAvgbetas_Y) == "ProbeSPONT")] <- "probes"
tt_SPONTvsIVF_M_Y$probes <- as.factor(tt_SPONTvsIVF_M_Y$probes)
placmet_MaleAvgbetas_Y$probes <- as.factor(placmet_MaleAvgbetas_Y$probes)
placmet_M_fulldata_Y <- merge(merge(tt_SPONTvsIVF_M_Y, placmet_MaleAvgbetas_Y[, c("deltaB","probes")], by = "probes"),
                              price_anno[,c("Closest_TSS_gene_name", "probes")], by = "probes")
write.table(tt_SPONTvsIVF_M_Y, sep = "\t", file = "./tt_SPONTvsIVF_M_Y.tsv")
write.table(placmet_M_fulldata_Y, sep = "\t", file = "./placmet_M_fulldata_Y.tsv")

#Female Autosomes 
SPONTvsIVF_modelF <- model.matrix(~ Group + Study, data = females)
SPONTvsIVF_fitF_auto <- lmFit(Mval_female_auto, SPONTvsIVF_modelF) 
SPONTvsIVF_fitF_auto <- eBayes(SPONTvsIVF_fitF_auto)
tt_SPONTvsIVF_F_all_auto <- topTable(SPONTvsIVF_fitF_auto, n = Inf, adjust = "fdr", coef = "GroupSpontaneous")
print(sum(tt_SPONTvsIVF_F_all_auto$adj.P.Val < 0.05)) #0
tt_SPONTvsIVF_F_all_auto$probes <- rownames(tt_SPONTvsIVF_F_all_auto)
names(placmet_FemaleAvgbetas_autosomes)[which(names(placmet_FemaleAvgbetas_autosomes) == "ProbeSPONT")] <- "probes"
tt_SPONTvsIVF_F_all_auto$probes <- as.factor(tt_SPONTvsIVF_F_all_auto$probes)
placmet_FemaleAvgbetas_autosomes$probes <- as.factor(placmet_FemaleAvgbetas_autosomes$probes)
placmet_F_fulldata_auto <- merge(merge(tt_SPONTvsIVF_F_all_auto, placmet_FemaleAvgbetas_autosomes[, c("deltaB","probes")], by = "probes"),
                                 price_anno[,c("Closest_TSS_gene_name", "probes")], by = "probes")
write.table(tt_SPONTvsIVF_F_all_auto, sep = "\t", file = "./tt_SPONTvsIVF_F_all_auto_study.tsv")
write.table(placmet_F_fulldata_auto, sep = "\t", file = "./placmet_F_fulldata_auto_study.tsv")
#Female X
SPONTvsIVF_fitF_X <- lmFit(Mval_female_X, SPONTvsIVF_modelF) 
SPONTvsIVF_fitF_X <- eBayes(SPONTvsIVF_fitF_X)
tt_SPONTvsIVF_F_X <- topTable(SPONTvsIVF_fitF_X, n = Inf, adjust = "fdr", coef = "GroupSpontaneous")
print(sum(tt_SPONTvsIVF_F_X$adj.P.Val < 0.05)) #0
tt_SPONTvsIVF_F_X$probes <- rownames(tt_SPONTvsIVF_F_X)
names(placmet_FemaleAvgbetas_X)[which(names(placmet_FemaleAvgbetas_X) == "ProbeSPONT")] <- "probes"
tt_SPONTvsIVF_F_X$probes <- as.factor(tt_SPONTvsIVF_F_X$probes)
placmet_FemaleAvgbetas_X$probes <- as.factor(placmet_FemaleAvgbetas_X$probes)
placmet_F_fulldata_X <- merge(merge(tt_SPONTvsIVF_F_X, placmet_FemaleAvgbetas_X[, c("deltaB","probes")], by = "probes"),
                              price_anno[,c("Closest_TSS_gene_name", "probes")], by = "probes")
write.table(tt_SPONTvsIVF_F_X, sep = "\t", file = "./tt_SPONTvsIVF_F_X_study.tsv")
write.table(placmet_F_fulldata_X, sep = "\t", file = "./placmet_F_fulldata_X_study.tsv")

#Plotting 
#Adding diff methylation information for colouring 
#Whole data 
placmet_wholepop_auto$diffmethylation <- "Not Biologically Significant"
placmet_wholepop_auto$diffmethylation[placmet_wholepop_auto$deltaB > 0.05 & placmet_wholepop_auto$adj.P.Val <0.05] <- "Hyper-methylated"
placmet_wholepop_auto$diffmethylation[placmet_wholepop_auto$deltaB < -0.05 & placmet_wholepop_auto$adj.P.Val <0.05] <- "Hypo-methylated"
#Male data 
placmet_M_fulldata_auto$diffmethylation <- "Not Biologically Significant"
placmet_M_fulldata_auto$diffmethylation[placmet_M_fulldata_auto$deltaB > 0.05 & placmet_M_fulldata_auto$adj.P.Val <0.05] <- "Hyper-methylated"
placmet_M_fulldata_auto$diffmethylation[placmet_M_fulldata_auto$deltaB < -0.05 & placmet_M_fulldata_auto$adj.P.Val <0.05] <- "Hypo-methylated"
placmet_M_fulldata_X$diffmethylation <- "Not Biologically Significant"
placmet_M_fulldata_X$diffmethylation[placmet_M_fulldata_X$deltaB > 0.05 & placmet_M_fulldata_X$adj.P.Val <0.05] <- "Hyper-methylated"
placmet_M_fulldata_X$diffmethylation[placmet_M_fulldata_X$deltaB < -0.05 & placmet_M_fulldata_X$adj.P.Val <0.05] <- "Hypo-methylated"
placmet_M_fulldata_Y$diffmethylation <- "Not Biologically Significant"
placmet_M_fulldata_Y$diffmethylation[placmet_M_fulldata_Y$deltaB > 0.05 & placmet_M_fulldata_Y$adj.P.Val <0.05] <- "Hyper-methylated"
placmet_M_fulldata_Y$diffmethylation[placmet_M_fulldata_Y$deltaB < -0.05 & placmet_M_fulldata_Y$adj.P.Val <0.05] <- "Hypo-methylated"
#Female Data 
placmet_F_fulldata_auto$diffmethylation <- "Not Biologically Significant"
placmet_F_fulldata_auto$diffmethylation[placmet_F_fulldata_auto$deltaB > 0.05 & placmet_F_fulldata_auto$adj.P.Val <0.05] <- "Hyper-methylated"
placmet_F_fulldata_auto$diffmethylation[placmet_F_fulldata_auto$deltaB < -0.05 & placmet_F_fulldata_auto$adj.P.Val <0.05] <- "Hypo-methylated"
placmet_F_fulldata_X$diffmethylation <- "Not Biologically Significant"
placmet_F_fulldata_X$diffmethylation[placmet_F_fulldata_X$deltaB > 0.05 & placmet_F_fulldata_X$adj.P.Val <0.05] <- "Hyper-methylated"
placmet_F_fulldata_X$diffmethylation[placmet_F_fulldata_X$deltaB < -0.05 & placmet_F_fulldata_X$adj.P.Val <0.05] <- "Hypo-methylated"

#Volcano Plots 
wholepop_auto <- ggplot(data = placmet_wholepop_auto, aes(x = deltaB, y = -log10(adj.P.Val), col = diffmethylation)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.3, size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  ylab("-log10(adjusted P.Value)") +
  xlab("Delta Beta") + 
  scale_y_continuous(breaks = seq(0, 3, by = 0.5), limits = c(0, 3)) +
  scale_x_continuous(breaks = seq(-0.35, 0.35, by = 0.1), limits = c(-0.35, 0.35)) +
  scale_color_manual(values = c("black", "#FFB518", "#0C7BDC"), 
                     labels = c("Not Biologically Significant","Hypo-Methylated", "Hyper-methylated"),
                     guide = "none")
placmet_M_fulldata_auto$siglabel <- ifelse(placmet_M_fulldata_auto$probes %in% maleauto_sig$probes, placmet_M_fulldata_auto$Closest_TSS_gene_name, NA)
male_auto <- ggplot(data = placmet_M_fulldata_auto, aes(x = deltaB, y = -log10(adj.P.Val), col = diffmethylation, label=siglabel)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.4, size = 3) + 
  geom_text_repel(max.overlaps = Inf,
                  segment.size = 0.5) +
  theme_bw() +
  ylab(" ") +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  xlab("Delta Beta") +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5), limits = c(0, 3)) +
  scale_x_continuous(breaks = seq(-0.35, 0.35, by = 0.1), limits = c(-0.35, 0.35)) +
  scale_color_manual(values = c("#FFB518", "#0C7BDC", "black"), 
                     labels = c("Not Biologically Significant", "Hyper-methylated", "Hypo-methylated"),
                     guide = "none")
female_auto <- ggplot(data = placmet_F_fulldata_auto, aes(x = deltaB, y = -log10(adj.P.Val), col = diffmethylation)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.4, size = 3) + 
  theme_bw() +
  ylab(" ") +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  xlab("Delta Beta") +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5), limits = c(0, 3)) +
  scale_x_continuous(breaks = seq(-0.35, 0.35, by = 0.1), limits = c(-0.35, 0.35)) +
  scale_color_manual(values = c("black", "#FFB518", "#0C7BDC"), 
                     labels = c("Not Biologically Significant", "Hyper-methylated", "Hypo-methylated"),
                     guide = "none")
png("./autosome_vol_adjFunnorm_panel.png", height = 9, width = 16, units = "in", res = 300)
grid.arrange(wholepop_auto, female_auto, male_auto, nrow = 1)
dev.off()


#plots of X chromosome 
placmet_M_fulldata_X$siglabel <- ifelse(placmet_M_fulldata_X$probes %in% maleX_sig$probes, placmet_M_fulldata_X$Closest_TSS_gene_name, NA)
male_X <- ggplot(data = placmet_M_fulldata_X, aes(x = deltaB, y = -log10(adj.P.Val), col = diffmethylation, label=siglabel)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.4, size = 3) +
  geom_text_repel(max.overlaps = Inf) +
  theme_bw() +
  ylab("-log10(adjusted P.Value)") +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  xlab("Delta Beta") +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5), limits = c(0, 3)) +
  scale_x_continuous(breaks = seq(-0.3, 0.3, by = 0.1), 
                     limits = c(-0.30,0.30), 
                     labels = label_number(accuracy = 0.1))  +
  scale_color_manual(values = c("black", "#FFB518", "#0C7BDC"), 
                     labels = c("Not Biologically Significant", "Hyper-methylated", "Hypo-methylated"),
                     guide = "none")

female_X <- ggplot(data = placmet_F_fulldata_X, aes(x = deltaB, y = -log10(adj.P.Val), col = diffmethylation)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.4, size = 3) +
  theme_bw() +
  ylab(" ") +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  xlab("Delta Beta") +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5), limits = c(0, 3)) +
  scale_x_continuous(breaks = seq(-0.3, 0.3, by = 0.1), 
                     limits = c(-0.30,0.30), 
                     labels = label_number(accuracy = 0.1)) +
  scale_color_manual(values = c("black", "#FFB518", "#0C7BDC"), 
                     labels = c("Not Biologically Significant", "Hyper-methylated", "Hypo-methylated"),
                     guide = "none")
png("./X_vol_adjFunnorm_panel.png", height = 9, width = 10, units = "in", res = 300)
grid.arrange(female_X, male_X, nrow = 1)
dev.off()

# Plot of Y chromosome
male_Y <- ggplot(data = placmet_M_fulldata_Y, aes(x = deltaB, y = -log10(adj.P.Val), col = diffmethylation)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.4, size = 3) + 
  theme_bw() +
  ylab("-log10(adjusted P.Value)") +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  xlab("Delta Beta") +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5), limits = c(0, 3)) +
  scale_x_continuous(breaks = seq(-0.2, 0.2, by = 0.1), limits = c(-0.2, 0.2)) +
  scale_color_manual(values = c("black", "#FFB518", "#0C7BDC"), 
                     labels = c("Not Biologically Significant", "Hyper-methylated", "Hypo-methylated"),
                     guide = "none")
png("./Y_vol_adjFunnorm.png", height = 9, width = 5, units = "in", res = 300)
male_Y
dev.off()