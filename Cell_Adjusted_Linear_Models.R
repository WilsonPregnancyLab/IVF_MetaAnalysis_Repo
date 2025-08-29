
#Data Table for Cell Decon Only Linear Models with Cell Type 
celldecon_probes <- read.csv("./Cell_Decon_Probelist.csv") #From Yuan et.al. cell deconvolution methodology paper 
placmet_adjFunnorm_mvalues_celldecon <- placmet_adjFunnorm_mvalues[rownames(placmet_adjFunnorm_mvalues) %in% celldecon_probes$cpg, ]
celldeconDNAm <- t(placmet_adjFunnorm_mvalues_celldecon)
celldeconDNAm2 <- as.data.frame(apply(celldeconDNAm, 2, as.numeric))
celldeconDNAm2$Sample_Name <- rownames(celldeconDNAm)
celldeconmergedd <- merge(pheno, estF_conceptandsex)
celldeconmergedd <- merge(celldeconmergedd, celldeconDNAm2)
head(celldeconmergedd[, 1:16], 5)
#Data Table for All Probe Linear Models with Cell Type (allDNAm2 comes from DNAmethylation_Analyses.R Script)
allmergedd_withcelltype <- merge(pheno, estF_conceptandsex)
allmergedd_withcelltype <- merge(allmergedd_withcelltype, allDNAm2)
head(allmergedd_withcelltype[, 1:16], 5)

# Function to run the model and extract adjusted p-values
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
# Define phenotype columns (adjust based on your dataset)
phenotype_cols <- 15:574
# Run models for the whole sample
celldecon_whole_adj_p <- run_model(celldeconmergedd, ~ Group + Fetal_Sex + Study, phenotype_cols)
# Run models for males only
celldecon_male_adj_p <- run_model(celldeconmergedd[celldeconmergedd$Fetal_Sex == "Male", ], ~ Group + Study, phenotype_cols)
# Run models for females only
celldecon_female_adj_p <- run_model(celldeconmergedd[celldeconmergedd$Fetal_Sex == "Female", ], ~ Group + Study, phenotype_cols)
# Run models for the whole sample
celldecon_whole_adj_p_celladjust <- run_model(celldeconmergedd, ~ Group + Fetal_Sex + Study + Trophoblasts + Stromal + Hofbauer + Endothelial + nRBC, phenotype_cols)
# Run models for males only
celldecon_male_adj_p_celladjust <- run_model(celldeconmergedd[celldeconmergedd$Fetal_Sex == "Male", ], ~ Group + Study + Trophoblasts + Stromal + Hofbauer + Endothelial + nRBC, phenotype_cols)
# Run models for females only
celldecon_female_adj_p_celladjust <- run_model(celldeconmergedd[celldeconmergedd$Fetal_Sex == "Female", ], ~ Group + Study + Trophoblasts + Stromal + Hofbauer + Endothelial + nRBC, phenotype_cols)
# Combine results into a data frame
celldecon_probe_results <- data.frame(
  CpG = colnames(celldeconmergedd[, phenotype_cols]),
  Whole_Raw_P = celldecon_whole_adj_p$Raw_P,
  Whole_Adj_P = celldecon_whole_adj_p$Adj_P,
  Whole_Raw_P_celladjust = celldecon_whole_adj_p_celladjust$Raw_P,
  Whole_Adj_P_celladjust = celldecon_whole_adj_p_celladjust$Adj_P,
  Male_Raw_P = celldecon_male_adj_p$Raw_P,
  Male_Adj_P = celldecon_male_adj_p$Adj_P,
  Male_Raw_P_celladjust = celldecon_male_adj_p_celladjust$Raw_P,
  Male_Adj_P_celladjust = celldecon_male_adj_p_celladjust$Adj_P,
  Female_Raw_P = celldecon_female_adj_p$Raw_P,
  Female_Adj_P = celldecon_female_adj_p$Adj_P,
  Female_Raw_P_celladjust = celldecon_female_adj_p_celladjust$Raw_P,
  Female_Adj_P_celladjust = celldecon_female_adj_p_celladjust$Adj_P
)
write.csv(celldecon_probe_results, "./celldecon_probe_LM_CellLM_results.csv")

# Define phenotype columns (adjust based on your dataset)
phenotype_cols <- 15:296559
# Run models for the whole sample
whole_adj_p_celladj <- run_model(allmergedd_withcelltype, ~ Group + Fetal_Sex + Study + Trophoblasts + Stromal + Hofbauer + Endothelial + nRBC, phenotype_cols)
# Run models for males only
male_adj_p_celladj <- run_model(allmergedd_withcelltype[allmergedd_withcelltype$Fetal_Sex == "Male", ], ~ Group + Study + Trophoblasts + Stromal + Hofbauer + Endothelial + nRBC, phenotype_cols)
# Run models for females only
female_adj_p_celladj <- run_model(allmergedd_withcelltype[allmergedd_withcelltype$Fetal_Sex == "Female", ], ~ Group + Study + Trophoblasts + Stromal + Hofbauer + Endothelial + nRBC, phenotype_cols)
# Combine results into a data frame
autosome_celladjusted_results <- data.frame(
  CpG = colnames(allmergedd_withcelltype[, phenotype_cols]),
  Whole_Raw_P = whole_adj_p_celladj$Raw_P,
  Whole_Adj_P = whole_adj_p_celladj$Adj_P,
  Male_Raw_P = male_adj_p_celladj$Raw_P,
  Male_Adj_P = male_adj_p_celladj$Adj_P,
  Female_Raw_P = female_adj_p_celladj$Raw_P,
  Female_Adj_P = female_adj_p_celladj$Adj_P
)
sum(autosome_celladjusted_results$Whole_Adj_P < 0.05, na.rm = TRUE) #138 CpGs vs 127 in original model
sum(autosome_celladjusted_results$Male_Adj_P < 0.05, na.rm = TRUE) #18 CpGs vs 9 in original model 
sum(autosome_celladjusted_results$Female_Adj_P < 0.05, na.rm = TRUE) #2 CpGs vs 0 in original model 

write.csv(autosome_celladjusted_results, "./autosome_celladjusted_results.csv")

#Making the Methylation Table with detlaB's, adjPvalues, and Gene names 
names(autosome_celladjusted_results)[names(autosome_celladjusted_results) == "CpG"] <- "Row.names"
price_anno <- read.csv("/workspace/lab/wilsonslab/lemairem/annotations/Price_anno_450K.tsv", header = TRUE, sep = "\t") #downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42409
names(price_anno)[names(price_anno) == "SPOT_ID"] <- "Row.names"
#create autosome results table with all infor
autosome_celladjusted_results_allinfo <- autosome_celladjusted_results %>%
  inner_join(placmet_FemaleAvgbetas_autosomes[, c("deltaBF", "Row.names")], by = "Row.names") %>%
  inner_join(placmet_MaleAvgbetas_autosomes[, c("deltaBM", "Row.names")], by = "Row.names") %>%
  inner_join(placmet_AllAvgbetas_autosomes[, c("deltaBA", "Row.names")], by = "Row.names") %>%
  inner_join(price_anno[, c("Closest_TSS_gene_name", "Row.names")], by = "Row.names")
#create methylation information columns 
#all
autosome_celladjusted_results_allinfo$All_diffmethylation <- "Not Biologically Significant"
autosome_celladjusted_results_allinfo$All_diffmethylation[autosome_celladjusted_results_allinfo$deltaBA > 0.05 & autosome_celladjusted_results_allinfo$Whole_Adj_P <0.05] <- "Hyper-methylated"
autosome_celladjusted_results_allinfo$All_diffmethylation[autosome_celladjusted_results_allinfo$deltaBA < -0.05 & autosome_celladjusted_results_allinfo$Whole_Adj_P <0.05] <- "Hypo-methylated"
#Male
autosome_celladjusted_results_allinfo$M_diffmethylation <- "Not Biologically Significant"
autosome_celladjusted_results_allinfo$M_diffmethylation[autosome_celladjusted_results_allinfo$deltaBM > 0.05 & autosome_celladjusted_results_allinfo$Male_Adj_P <0.05] <- "Hyper-methylated"
autosome_celladjusted_results_allinfo$M_diffmethylation[autosome_celladjusted_results_allinfo$deltaBM < -0.05 & autosome_celladjusted_results_allinfo$Male_Adj_P <0.05] <- "Hypo-methylated"
#Female
autosome_celladjusted_results_allinfo$F_diffmethylation <- "Not Biologically Significant"
autosome_celladjusted_results_allinfo$F_diffmethylation[autosome_celladjusted_results_allinfo$deltaBF > 0.05 & autosome_celladjusted_results_allinfo$Female_Adj_P <0.05] <- "Hyper-methylated"
autosome_celladjusted_results_allinfo$F_diffmethylation[autosome_celladjusted_results_allinfo$deltaBF < -0.05 & autosome_celladjusted_results_allinfo$Female_Adj_P <0.05] <- "Hypo-methylated"

write.csv(autosome_celladjusted_results_allinfo, "./autosome_celladjusted_results_allinfo.csv")

#Plotting 
wholepop_auto <- ggplot(data = autosome_celladjusted_results_allinfo, aes(x = deltaBA, y = -log10(Whole_Adj_P), col = All_diffmethylation)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.3, size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  ylab("-log10(adjusted P.Value)") +
  xlab("Delta Beta") + 
  scale_y_continuous(breaks = seq(0, 6.5, by = 1), limits = c(0, 6.5)) +
  scale_x_continuous(breaks = seq(-0.4, 0.4, by = 0.1), limits = c(-0.4, 0.4)) +
  scale_color_manual(values = c("#FFB518", "#0C7BDC", "black"), 
                     labels = c("Not Biologically Significant","Hypo-Methylated", "Hyper-methylated"),
                     guide = "none")
malepop_auto <- ggplot(data = autosome_celladjusted_results_allinfo, aes(x = deltaBM, y = -log10(Male_Adj_P), col = M_diffmethylation)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.3, size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  ylab(" ") +
  xlab("Delta Beta") + 
  scale_y_continuous(breaks = seq(0, 6.5, by = 1), limits = c(0, 6.5)) +
  scale_x_continuous(breaks = seq(-0.4, 0.4, by = 0.1), limits = c(-0.4, 0.4)) +
  scale_color_manual(values = c("#FFB518", "#0C7BDC", "black"), 
                     labels = c("Not Biologically Significant","Hypo-Methylated", "Hyper-methylated"),
                     guide = "none")
femalepop_auto <- ggplot(data = autosome_celladjusted_results_allinfo, aes(x = deltaBF, y = -log10(Female_Adj_P), col = F_diffmethylation)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.3, size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  ylab(" ") +
  xlab("Delta Beta") + 
  scale_y_continuous(breaks = seq(0, 6.5, by = 1), limits = c(0, 6.5)) +
  scale_x_continuous(breaks = seq(-0.4, 0.4, by = 0.1), limits = c(-0.4, 0.4)) +
  scale_color_manual(values = c("black","#FFB518", "#0C7BDC"), 
                     labels = c("Not Biologically Significant","Hypo-Methylated", "Hyper-methylated"),
                     guide = "none")
png("./autosome_celladjusted_panel.png", height = 9, width = 16, units = "in", res = 300)
grid.arrange(wholepop_auto, femalepop_auto, malepop_auto, nrow = 1)
dev.off()

#Overlap between Non-adjusted and cell adjusted models 
sig_whole_autosome_celladjusted_results <- autosome_celladjusted_results[autosome_celladjusted_results$Whole_Adj_P < 0.05,] #138
sig_male_autosome_celladjusted_results <- autosome_celladjusted_results[autosome_celladjusted_results$Male_Adj_P < 0.05,] #18
sig_female_autosome_celladjusted_results <- autosome_celladjusted_results[autosome_celladjusted_results$Female_Adj_P < 0.05, ] #2 

adjusted_overlap_whole <- full_join(sig_whole_autosome_celladjusted_results, sig_whole_autosome_results, by = "Row.names") 
sum(sig_whole_autosome_celladjusted_results$Row.names %in% sig_whole_autosome_results$Row.names) #90 
write.csv(adjusted_overlap_whole, "./adjusted_overlap_whole.csv")
adjusted_overlap_male <- full_join(sig_male_autosome_celladjusted_results, sig_male_autosome_results, by = "Row.names")
sum(sig_male_autosome_celladjusted_results$Row.names %in% sig_male_autosome_results$Row.names) #8 
write.csv(adjusted_overlap_male, "./adjusted_overlap_male.csv")
adjusted_overlap_female <- full_join(sig_female_autosome_celladjusted_results, sig_female_autosome_results, by = "Row.names")
sum(sig_female_autosome_celladjusted_results$Row.names %in% sig_female_autosome_results$Row.names) #0
write.csv(adjusted_overlap_female, "./adjusted_overlap_female.csv")
