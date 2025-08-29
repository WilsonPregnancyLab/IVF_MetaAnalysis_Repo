#R version 4.4.2
library(stringr) #version 1.5.1
library(dplyr) #version 1.1.4
library(ggplot2) #version 3.5.1
library(gridExtra) #version 2.3
library(scales) #version 1.3.0

#creating groups for sex stratification 
metadata <- read.csv("./PlacMet_MetaData_outliersremoved.csv")
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

placmet_filtfun_F <- placmet_adjFunnorm_filtbetas[, females$Sample_Name]
# 281 females - gives all the GSMs in either IVF or spont that are female
placmet_filtfun_M <- placmet_adjFunnorm_filtbetas[, males$Sample_Name]
# 294 males - gives all the GSMs in either IVF or spont that are male 

#creating groups for conception stratification 
Spontaneousmetadata_all <- subset(metadata, metadata$Group == "Spontaneous")
Spontaneousmetadata_F <- subset(females, females$Group == "Spontaneous")
Spontaneousmetadata_M <- subset(males, males$Group == "Spontaneous")
IVFmetadata_all <- subset(metadata, metadata$Group == "IVF")
IVFmetadata_F <- subset(females, females$Group == "IVF")
IVFmetadata_M <- subset(males, males$Group == "IVF")

#Fix Beta Tables to be (IVF - Spontaneous (final - initial))
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
placmet_AllAvgbetas_autosomes$deltaBA <- placmet_AllAvgbetas_autosomes$AvgBIVF - placmet_AllAvgbetas_autosomes$AvgBSPONT 
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
placmet_MaleAvgbetas_autosomes$deltaBM <- placmet_MaleAvgbetas_autosomes$AvgBIVF - placmet_MaleAvgbetas_autosomes$AvgBSPONT 
rownames(placmet_MaleAvgbetas_autosomes) <- placmet_MaleAvgbetas_autosomes$ProbeSPONT
write.table(placmet_MaleAvgbetas_autosomes, sep = "\t", file = "./placmet_MaleAvgbetas_autosomes.tsv")
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
placmet_FemaleAvgbetas_autosomes$deltaBF <- placmet_FemaleAvgbetas_autosomes$AvgBIVF - placmet_FemaleAvgbetas_autosomes$AvgBSPONT 
rownames(placmet_FemaleAvgbetas_autosomes) <- placmet_FemaleAvgbetas_autosomes$ProbeSPONT
write.table(placmet_FemaleAvgbetas_autosomes, sep = "\t", file = "./placmet_FemaleAvgbetas_autosomes.tsv")

#Making the Methylation Table with detlaB's, adjPvalues, and Gene names 
names(autosome_results)[names(autosome_results) == "CpG"] <- "Row.names"
price_anno <- read.csv("./Price_anno_450K.tsv", header = TRUE, sep = "\t") #downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42409
names(price_anno)[names(price_anno) == "SPOT_ID"] <- "Row.names"
#create autosome results table with all infor
autosome_results_allinfo <- autosome_results %>%
  inner_join(placmet_FemaleAvgbetas_autosomes[, c("deltaBF", "Row.names")], by = "Row.names") %>%
  inner_join(placmet_MaleAvgbetas_autosomes[, c("deltaBM", "Row.names")], by = "Row.names") %>%
  inner_join(placmet_AllAvgbetas_autosomes[, c("deltaBA", "Row.names")], by = "Row.names") %>%
  inner_join(price_anno[, c("Closest_TSS_gene_name", "Row.names")], by = "Row.names")
#create methylation information columns 
#all
autosome_results_allinfo$All_diffmethylation <- "Not Biologically Significant"
autosome_results_allinfo$All_diffmethylation[autosome_results_allinfo$deltaBA > 0.05 & autosome_results_allinfo$Whole_Adj_P <0.05] <- "Hyper-methylated"
autosome_results_allinfo$All_diffmethylation[autosome_results_allinfo$deltaBA < -0.05 & autosome_results_allinfo$Whole_Adj_P <0.05] <- "Hypo-methylated"
#Male
autosome_results_allinfo$M_diffmethylation <- "Not Biologically Significant"
autosome_results_allinfo$M_diffmethylation[autosome_results_allinfo$deltaBM > 0.05 & autosome_results_allinfo$Male_Adj_P <0.05] <- "Hyper-methylated"
autosome_results_allinfo$M_diffmethylation[autosome_results_allinfo$deltaBM < -0.05 & autosome_results_allinfo$Male_Adj_P <0.05] <- "Hypo-methylated"
#Female
autosome_results_allinfo$F_diffmethylation <- "Not Biologically Significant"
autosome_results_allinfo$F_diffmethylation[autosome_results_allinfo$deltaBF > 0.05 & autosome_results_allinfo$Female_Adj_P <0.05] <- "Hyper-methylated"
autosome_results_allinfo$F_diffmethylation[autosome_results_allinfo$deltaBF < -0.05 & autosome_results_allinfo$Female_Adj_P <0.05] <- "Hypo-methylated"

write.csv(autosome_results_allinfo, "./autosome_results_allinfo.csv")

#plotting Autosomes 
wholepop_auto <- ggplot(data = autosome_results_allinfo, aes(x = deltaBA, y = -log10(Whole_Adj_P), col = All_diffmethylation)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.3, size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  ylab("-log10(adjusted P.Value)") +
  xlab("Delta Beta") + 
  scale_y_continuous(breaks = seq(0, 5.5, by = 1), limits = c(0, 5.5)) +
  scale_x_continuous(breaks = seq(-0.4, 0.4, by = 0.1), limits = c(-0.4, 0.4)) +
  scale_color_manual(values = c("#FFB518", "#0C7BDC", "black"), 
                     labels = c("Not Biologically Significant","Hypo-Methylated", "Hyper-methylated"),
                     guide = "none")
malepop_auto <- ggplot(data = autosome_results_allinfo, aes(x = deltaBM, y = -log10(Male_Adj_P), col = M_diffmethylation)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.3, size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  ylab(" ") +
  xlab("Delta Beta") + 
  scale_y_continuous(breaks = seq(0, 5.5, by = 1), limits = c(0, 5.5)) +
  scale_x_continuous(breaks = seq(-0.4, 0.4, by = 0.1), limits = c(-0.4, 0.4)) +
  scale_color_manual(values = c("#FFB518", "#0C7BDC", "black"), 
                     labels = c("Not Biologically Significant","Hypo-Methylated", "Hyper-methylated"),
                     guide = "none")
femalepop_auto <- ggplot(data = autosome_results_allinfo, aes(x = deltaBF, y = -log10(Female_Adj_P), col = F_diffmethylation)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.3, size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  ylab(" ") +
  xlab("Delta Beta") + 
  scale_y_continuous(breaks = seq(0, 5.5, by = 1), limits = c(0, 5.5)) +
  scale_x_continuous(breaks = seq(-0.4, 0.4, by = 0.1), limits = c(-0.4, 0.4)) +
  scale_color_manual(values = c("black", "#FFB518", "#0C7BDC"), 
                     labels = c("Not Biologically Significant","Hypo-Methylated", "Hyper-methylated"),
                     guide = "none")

png("./autosome_all_panel.png", height = 9, width = 16, units = "in", res = 300)
grid.arrange(wholepop_auto, femalepop_auto, malepop_auto, nrow = 1)
dev.off()



#X Chromosome Beta Table
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
placmet_MaleAvgbetas_X$deltaBM <- placmet_MaleAvgbetas_X$AvgBIVF - placmet_MaleAvgbetas_X$AvgBSPONT #This is where the equation is flipped CORRECTLY 
rownames(placmet_MaleAvgbetas_X) <- placmet_MaleAvgbetas_X$ProbeSPONT
write.table(placmet_MaleAvgbetas_X, sep = "\t", file = "./placmet_MaleAvgbetas_X.tsv")
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
placmet_FemaleAvgbetas_X$deltaBF <- placmet_FemaleAvgbetas_X$AvgBIVF - placmet_FemaleAvgbetas_X$AvgBSPONT #This is where the equation is flipped CORRECTLY 
rownames(placmet_FemaleAvgbetas_X) <- placmet_FemaleAvgbetas_X$ProbeSPONT
write.table(placmet_FemaleAvgbetas_X, sep = "\t", file = "./placmet_FemaleAvgbetas_X.tsv")
#create X Chromosome results table with all information
names(XChromosome_results)[names(XChromosome_results) == "CpG"] <- "Row.names"
XChromosome_results_allinfo <- XChromosome_results %>%
  inner_join(placmet_FemaleAvgbetas_X[, c("deltaBF", "Row.names")], by = "Row.names") %>%
  inner_join(placmet_MaleAvgbetas_X[, c("deltaBM", "Row.names")], by = "Row.names") %>%
  inner_join(price_anno[, c("Closest_TSS_gene_name", "Row.names")], by = "Row.names")
#create methylation information columns 
#Male
XChromosome_results_allinfo$M_diffmethylation <- "Not Biologically Significant"
XChromosome_results_allinfo$M_diffmethylation[XChromosome_results_allinfo$deltaBM > 0.05 & XChromosome_results_allinfo$Males_X <0.05] <- "Hyper-methylated"
XChromosome_results_allinfo$M_diffmethylation[XChromosome_results_allinfo$deltaBM < -0.05 & XChromosome_results_allinfo$Males_X <0.05] <- "Hypo-methylated"
#Female
XChromosome_results_allinfo$F_diffmethylation <- "Not Biologically Significant"
XChromosome_results_allinfo$F_diffmethylation[XChromosome_results_allinfo$deltaBF > 0.05 & XChromosome_results_allinfo$Females_X <0.05] <- "Hyper-methylated"
XChromosome_results_allinfo$F_diffmethylation[XChromosome_results_allinfo$deltaBF < -0.05 & XChromosome_results_allinfo$Females_X <0.05] <- "Hypo-methylated"

write.csv(XChromosome_results_allinfo, "./XChromosome_results_allinfo")
#Plot
malepop_X <- ggplot(data = XChromosome_results_allinfo, aes(x = deltaBM, y = -log10(MaleX_Adj_P), col = M_diffmethylation)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.3, size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  ylab(" ") +
  xlab("Delta Beta") + 
  scale_y_continuous(breaks = seq(0, 2, by = 0.5), limits = c(0, 2)) +
  scale_x_continuous(breaks = seq(-0.30, 0.30, by = 0.1), limits = c(-0.30, 0.30), labels = label_number(accuracy = 0.1)) +
  scale_color_manual(values = c("black","black", "black"), 
                     labels = c("Not Biologically Significant","Hypo-Methylated", "Hyper-methylated"),
                     guide = "none")
femalepop_X <- ggplot(data = XChromosome_results_allinfo, aes(x = deltaBF, y = -log10(FemaleX_Adj_P), col = F_diffmethylation)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.3, size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  ylab("-log10(adjusted P.Value") +
  xlab("Delta Beta") + 
  scale_y_continuous(breaks = seq(0, 2, by = 0.5), limits = c(0, 2)) +
  scale_x_continuous(breaks = seq(-0.30, 0.30, by = 0.1), limits = c(-0.30, 0.30), labels = label_number(accuracy = 0.1)) +
  scale_color_manual(values = c("black", "black", "black"), 
                     labels = c("Not Biologically Significant","Hypo-Methylated", "Hyper-methylated"),
                     guide = "none")

png("./X_chromosome_volcano_panel.png", height = 9, width = 10, units = "in", res = 300)
grid.arrange(femalepop_X, malepop_X, nrow = 1)
dev.off()


#Y Chromosome BetaTable 
placmet_IVF_M_Y <- as.data.frame(placmet_filtfun_M[rownames(placmet_filtfun_M) %in% chrYprobes$probeID, IVFmetadata_M$Sample_Name]) #dim: 227, 56
placmet_SPONT_M_Y <- as.data.frame(placmet_filtfun_M[rownames(placmet_filtfun_M) %in% chrYprobes$probeID, Spontaneousmetadata_M$Sample_Name]) #dim: 227, 238
# calculate male Y average betas
placmet_SPONT_M_Y$AvgBSPONT <- rowMeans(placmet_SPONT_M_Y)
placmet_SPONT_M_Y$ProbeSPONT <- rownames(placmet_SPONT_M_Y)
placmet_IVF_M_Y$AvgBIVF <- rowMeans(placmet_IVF_M_Y)
placmet_IVF_M_Y$ProbeIVF <- rownames(placmet_IVF_M_Y)
# Merge male Y table
placmet_MaleAvgbetas_Y<- merge(placmet_SPONT_M_Y[,c("AvgBSPONT","ProbeSPONT")], placmet_IVF_M_Y[,c("AvgBIVF","ProbeIVF")], by = "row.names")
placmet_MaleAvgbetas_Y$deltaBM <- placmet_MaleAvgbetas_Y$AvgBIVF - placmet_MaleAvgbetas_Y$AvgBSPONT #This is where the equation is flipped CORRECTLY 
rownames(placmet_MaleAvgbetas_Y) <- placmet_MaleAvgbetas_Y$ProbeSPONT
write.table(placmet_MaleAvgbetas_Y, sep = "\t", file = "./placmet_MaleAvgbetas_Y.tsv")  
#Plotting 
names(YChromosome_results)[names(YChromosome_results) == "CpG"] <- "Row.names"
#create Y Chromosome results table with all inforation
YChromosome_results_allinfo <- YChromosome_results %>%
  inner_join(placmet_MaleAvgbetas_Y[, c("deltaBM", "Row.names")], by = "Row.names") %>%
  inner_join(price_anno[, c("Closest_TSS_gene_name", "Row.names")], by = "Row.names")
#create methylation information columns 
#Male
YChromosome_results_allinfo$M_diffmethylation <- "Not Biologically Significant"
YChromosome_results_allinfo$M_diffmethylation[YChromosome_results_allinfo$deltaBM > 0.05 & YChromosome_results_allinfo$Males_Y <0.05] <- "Hyper-methylated"
YChromosome_results_allinfo$M_diffmethylation[YChromosome_results_allinfo$deltaBM < -0.05 & YChromosome_results_allinfo$Males_Y <0.05] <- "Hypo-methylated"
write.csv(YChromosome_results_allinfo, "./YChromosome_results_allinfo.csv")

malepop_Y <- ggplot(data = YChromosome_results_allinfo, aes(x = deltaBM, y = -log10(MaleY_Adj_P), col = M_diffmethylation)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.3, size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  ylab(" ") +
  xlab("Delta Beta") + 
  scale_y_continuous(breaks = seq(0, 2, by = 0.5), limits = c(0, 2)) +
  scale_x_continuous(breaks = seq(-0.30, 0.30, by = 0.1), limits = c(-0.30, 0.30), labels = label_number(accuracy = 0.1)) +
  scale_color_manual(values = c("black","black", "black"), 
                     labels = c("Not Biologically Significant","Hypo-Methylated", "Hyper-methylated"),
                     guide = "none")

png("./Y_chromosome_volcano.png", height = 9, width = 5, units = "in", res = 300)
malepop_Y
dev.off() 
