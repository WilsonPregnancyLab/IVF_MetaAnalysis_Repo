#R version 4.4.2
library(limma) #version 3.60.4

#Checking validity of these models when you add in different changes 
#Uses the run_model function made in DNAMethylation_Analysis.R

# Run models for the whole sample excluding study with great imbalance
whole_adj_pEX <- run_model(allmergedd[allmergedd$Study!="GSE75248",], ~ Group + Fetal_Sex + Study, phenotype_cols)
# Run models for males only
male_adj_pEX <- run_model(allmergedd[allmergedd$Study!="GSE75248" & allmergedd$Fetal_Sex == "Male", ], ~ Group + Study, phenotype_cols)
# Run models for females only
female_adj_pEX <- run_model(allmergedd[allmergedd$Study!="GSE75248" & allmergedd$Fetal_Sex == "Female", ], ~ Group + Study, phenotype_cols)
# Combine results into a data frame
autosome_resultsEX <- data.frame(
  CpG = colnames(allmergedd[, phenotype_cols]),
  Whole_Raw_P = whole_adj_pEX$Raw_P,
  Whole_Adj_P = whole_adj_pEX$Adj_P,
  Male_Raw_P = male_adj_pEX$Raw_P,
  Male_Adj_P = male_adj_pEX$Adj_P,
  Female_Raw_P = female_adj_pEX$Raw_P,
  Female_Adj_P = female_adj_pEX$Adj_P
)
sum(autosome_resultsEX$Whole_Adj_P < 0.05, na.rm = TRUE) #108
sum(autosome_resultsEX$Male_Adj_P < 0.05, na.rm = TRUE) #24
sum(autosome_resultsEX$Female_Adj_P < 0.05, na.rm = TRUE) #0 
write.csv(autosome_resultsEX, "linearmodel_exclude_GSE75248_autosome_results.csv", row.names = FALSE)

#Repeated analyses on significant associations from the initial autosomal analysis 
#Table of only significant probes for this linear model 
significant_autosome_results <- autosome_results[autosome_results$Whole_Adj_P < 0.05, ]
#Significant probes only
placmet_adjFunnorm_mvalues_significant <- as.data.frame(
  placmet_adjFunnorm_mvalues[rownames(placmet_adjFunnorm_mvalues) %in% 
                               c(significant_autosome_results$Row.names), ]) #127 probes 
sigDNAm <- t(placmet_adjFunnorm_mvalues_significant)
sigDNAm2 <- as.data.frame(apply(sigDNAm, 2, as.numeric))
sigDNAm2$Sample_Name <- rownames(sigDNAm)
sigmergedd <- merge(pheno, sigDNAm2)
head(sigmergedd[, 1:10], 5) 
#reset phenotype cols
phenotype_cols <- 5:131
# Run autosome models for the whole sample
whole_adj_p_sig <- run_model(sigmergedd, ~ Group + Fetal_Sex + Study, phenotype_cols)
# Run autosome models for males only
male_adj_p_sig <- run_model(sigmergedd[sigmergedd$Fetal_Sex == "Male", ], ~ Group + Study, phenotype_cols)
# Run autosome models for females only
female_adj_p_sig <- run_model(sigmergedd[sigmergedd$Fetal_Sex == "Female", ], ~ Group + Study, phenotype_cols)
# Combine results into a data frame
sig_results <- data.frame(
  CpG = colnames(sigmergedd[, phenotype_cols]),
  Whole_Raw_P = whole_adj_p_sig$Raw_P,
  Whole_Adj_P = whole_adj_p_sig$Adj_P,
  Male_Raw_P = male_adj_p_sig$Raw_P,
  Male_Adj_P = male_adj_p_sig$Adj_P,
  Female_Raw_P = female_adj_p_sig$Raw_P,
  Female_Adj_P = female_adj_p_sig$Adj_P
)
#checking how many significant hits found 
sum(sig_results$Whole_Adj_P < 0.05, na.rm = TRUE) #127
sum(sig_results$Male_Adj_P < 0.05, na.rm = TRUE) #126
sum(sig_results$Female_Adj_P < 0.05, na.rm = TRUE) #118 
write.csv(sig_results, "./sig_results.csv")

phenotype_cols <- 5:296549  
whole_adj_pr <- run_model(allmergedd, ~ Group + Fetal_Sex + Study, phenotype_cols, robust=T)
# Run models for males only
male_adj_pr <- run_model(allmergedd[allmergedd$Fetal_Sex == "Male", ], ~ Group + Study, phenotype_cols, robust=T)
# Run models for females only
female_adj_pr <- run_model(allmergedd[allmergedd$Fetal_Sex == "Female", ], ~ Group + Study, phenotype_cols, robust=T)
# Combine results into a data frame
all_resultsR <- data.frame(
  CpG = colnames(allmergedd[, phenotype_cols]),
  Whole_Raw_P = whole_adj_pr$Raw_P,
  Whole_Adj_P = whole_adj_pr$Adj_P,
  Male_Raw_P = male_adj_pr$Raw_P,
  Male_Adj_P = male_adj_pr$Adj_P,
  Female_Raw_P = female_adj_pr$Raw_P,
  Female_Adj_P = female_adj_pr$Adj_P
)
sum(all_resultsR$Whole_Adj_P < 0.05, na.rm = TRUE) #127
sum(all_resultsR$Male_Adj_P < 0.05, na.rm = TRUE) #126
sum(all_resultsR$Female_Adj_P < 0.05, na.rm = TRUE) #118
write.csv(sig_resultsR, "linearmodel_robust_sig_results.csv", row.names = FALSE)

  
whole_LM_validity <- autosome_results_allinfo %>%
  select(
    Row.names, 
    Closest_TSS_gene_name, 
    deltaBA,
    Adj_P_LM = Whole_Adj_P
  ) %>%
  left_join(
    autosome_resultsEX %>% 
      select(
        Adj_P_EX = Whole_Adj_P,
        CpG
      ),
    by = c("Row.names" = "CpG")) %>%
  left_join(
    all_resultsR %>%
      select(
        Adj_P_R = Whole_Adj_P,
        CpG
      ),
    by = c("Row.names" = "CpG")) 

top15 <- whole_LM_validity[order(whole_LM_validity$Adj_P_LM), ][1:15, ]
write.csv(top15, "./top15_sig_validity.csv")
