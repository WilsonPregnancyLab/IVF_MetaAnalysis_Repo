#R Version - R4.4.1
#Packages
library(minfi) #version 1.50.0
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) #version 0.6.1
library(IlluminaHumanMethylation450kmanifest) #version 0.4.0
library(magrittr) #version 2.0.3
library(wateRmelon). #version 2.11.2

#Load Base Directories (each directory should have all red and green IDATs + one sample sheet as a .csv)
baseDir120250 <- "/workspace/lab/wilsonslab/lemairem/Placenta_Based_Studies/Met_array/rerun_July2024_XYprobes/GSE120250_dir"
baseDir75248 <- "/workspace/lab/wilsonslab/lemairem/Placenta_Based_Studies/Met_array/rerun_July2024_XYprobes/GSE75248_dir"
baseDir208529 <- "/workspace/lab/wilsonslab/lemairem/Placenta_Based_Studies/Met_array/rerun_July2024_XYprobes/GSE208529_dir"
#Create targets 
targets120250 <- read.metharray.sheet(baseDir120250)
targets75248 <- read.metharray.sheet(baseDir75248)
targets208529 <- read.metharray.sheet(baseDir208529)

#Sex Prediction - Quality control 
GSE75248_sexpredict_RGset <- read.metharray.exp(targets = targets75248, verbose = TRUE)
GSE75248_MetSet_raw <- preprocessRaw(GSE75248_sexpredict_RGset)
GSE75248_GMetSet_raw <- mapToGenome(GSE75248_MetSet_raw)
GSE75248_sex <- getSex(GSE75248_GMetSet_raw)
CompareSex75248 <- merge(GSE75248_sex[,c("predictedSex")], GSE75248_metadata_sex[,c("reportedSex")], by = "row.names")
CompareSex75248$match <- CompareSex75248$predictedSex == CompareSex75248$reportedSex

GSE120250_sexpredict_RGset <- read.metharray.exp(targets = targets120250, verbose = TRUE)
GSE120250_MetSet_raw <- preprocessRaw(GSE120250_sexpredict_RGset)
GSE120250_GMetSet_raw <- mapToGenome(GSE120250_MetSet_raw)
GSE120250_sex <- getSex(GSE120250_GMetSet_raw)
CompareSex120250 <- merge(GSE120250_sex[,c("predictedSex")], GSE120250_metadata_sex[,c("reportedSex")], by = "row.names")
CompareSex120250$match <- CompareSex120250$predictedSex == CompareSex120250$reportedSex

GSE208529_sexpredict_RGset <- read.metharray.exp(targets = targets208529, verbose = TRUE)
GSE208529_MetSet_raw <- preprocessRaw(GSE208529_sexpredict_RGset)
GSE208529_GMetSet_raw <- mapToGenome(GSE208529_MetSet_raw)
GSE208529_sex <- getSex(GSE208529_GMetSet_raw)
CompareSex208529 <- merge(GSE208529_sex[,c("predictedSex")], GSE208529_metadata_sex[,c("reportedSex")], by = "row.names")
CompareSex208529$match <- CompareSex208529$predictedSex == CompareSex208529$reportedSex

#Edit targets to remove "EXCLUDE" samples from sample sheets 
edit_targets120250 <- subset(targets120250, Sample_Group != "EXCLUDE")
edit_targets75248 <- subset(targets75248, Sample_Group != "EXCLUDE")
edit_targets208529 <- subset(targets208529, Sample_Group != "EXCLUDE")
# Samples that were excluded
## Fetal sex did not match predicted (n=5) 
## in vivo ART was used (n = 21)
## Maternal facing placental samples (n = 195)
## Did not report conception type (n = 9)
## Outliers from 120250 (n = 9)

#Build RG sets 
RGset120250 <- read.metharray.exp(targets = edit_targets120250, verbose = TRUE) 
RGset75248 <- read.metharray.exp(targets = edit_targets75248, verbose = TRUE)
RGset208529 <- read.metharray.exp(targets = edit_targets208529, verbose = TRUE)
#Convert EPIC array (GSE208529) to a 450K array (removes EPIC only probes)
convertedRGset208529 <- convertArray(RGset208529, outType = "IlluminaHumanMethylation450k") 
#Combine the RG sets (can only do this two at a time)
combo_120240_75248 <- combineArrays(RGset120250, RGset75248, outType = c("IlluminaHumanMethylation450k"), verbose = TRUE)
placmet_combined_RGset <- combineArrays(combo_120240_75248, convertedRGset208529, outType = c("IlluminaHumanMethylation450k"), verbose = TRUE)
#Check number of probes (type I and type II probes added together = 485,512 for 450K)
placmet_combined_manifest <- getManifest(placmet_combined_RGset)
print(placmet_combined_manifest)


#Normalization (adjFunnorm)
placmet_adjFunnorm <- adjustedFunnorm(placmet_combined_RGset)


#Adapted filtering step to keep both autosomal and XY probes

#Filtering bad probes (bad detection p vals or missing betas)
detp <- minfi::detectionP(placmet_combined_RGset)
number_bad_P_before_Fstrat <- print(sum(rowSums(detp)>=(ncol(placmet_combined_RGset))*0.05)) #2498
write.table(detp, sep = "/t", file = "./detp_table.tsv")
## Annotations for X and Y probes
probeInfo <- as.data.frame(cbind(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations, 
                                 IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other, 
                                 IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Manifest)) 
probeInfo$probeID <- rownames(probeInfo)
chrXprobes <- subset(probeInfo, probeInfo$chr == "chrX")
chrYprobes <- subset(probeInfo, probeInfo$chr == "chrY")
metadata <- read.csv("/workspace/lab/wilsonslab/lemairem/Placenta_Based_Studies/Met_array/rerun_June2024_outlier_exclusion/PlacMet_MetaData_outliersremoved.csv")
males <- subset(metadata, metadata$Fetal_Sex == "Male")
females <- subset(metadata, metadata$Fetal_Sex == "Female")
## For females, set detp in Y chromosomes to 0 (these probes do not bind to anything in female samples)
detp[rownames(detp) %in% chrYprobes$probeID, females$Sample_Name] <- 0
head(detp[rownames(detp) %in% chrYprobes$probeID, females$Sample_Name]) #Make sure all 0
bad_detp <- detp > 0.01
#running bad probes
number_bad_detp <- print(sum(rowSums(bad_detp)>=(ncol(placmet_combined_RGset))*0.05)) #2189
# missing betas >5% of samples 
avgbeta <- getBeta(placmet_combined_RGset)
bad_beta <- is.na(avgbeta)
number_bad_beta <- print(sum(rowSums(bad_beta)>=(ncol(placmet_combined_RGset))*0.05)) #29
# remove probes with bad p-values or missing betas
badProbes <- rowSums(bad_detp)>=(ncol(placmet_combined_RGset))*0.05 | rowSums(bad_beta)>=(ncol(placmet_combined_RGset))*0.05 
number_badprobes <- print(sum(rowSums(badProbes)))
placmet_adjFunnorm_BPfilt <- placmet_adjFunnorm[!badProbes,]

# Remove SNP probes
price_anno <- read.csv("/workspace/lab/wilsonslab/lemairem/annotations/Price_anno_450K.tsv", header = TRUE, sep = "\t") #downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42409
price_anno_SNP <- subset(price_anno, price_anno$n_target.CpG.SNP>0)
placmet_adjFunnorm_SNPremoved <- placmet_adjFunnorm_BPfilt[!rownames(placmet_adjFunnorm_BPfilt) %in% price_anno_SNP$ID,]

# Remove Cross-Hybridizing probes
price_anno_Hyb <- subset(price_anno, price_anno$XY_Hits == "XY_YES" | price_anno$Autosomal_Hits == "A_YES")
placmet_adjFunnorm_HybRemoved <- placmet_adjFunnorm_SNPremoved[!rownames(placmet_adjFunnorm_SNPremoved) %in% price_anno_Hyb$ID,]

# Remove Non-variable placental probes 
nonvarplac_anno <- read.csv("/workspace/lab/wilsonslab/lemairem/annotations/Invariant_Placenta_CpGs.csv", sep = ",") #downloaded from https://github.com/redgar598/Tissue_Nonvariable_450K_CpGs
placmet_adjFunnorm_nonvarremoved <- placmet_adjFunnorm_HybRemoved[!rownames(placmet_adjFunnorm_HybRemoved) %in% nonvarplac_anno$CpG,]
placmet_adjFunnorm_allfiltered <- placmet_adjFunnorm_nonvarremoved

