sa <- autosome_results_allinfo %>%
  select(
    Row.names,
    Closest_TSS_gene_name,
    Whole_LM_Adj_P = Whole_Adj_P,
    Whole_LM_Delta = deltaBA,
    Male_LM_Adj_P = Male_Adj_P,
    Male_LM_Delta = deltaBM,
    Female_LM_Adj_P = Female_Adj_P,
    Female_LM_Delta = deltaBF
  )
sa_celladj <- autosome_celladjusted_results_allinfo %>%
  select(
    Row.names, 
    Whole_celladjustLM_Adj_P = Whole_Adj_P,
    Whole_celladjustLM_Delta = deltaBA,
    Male_celladjustLM_Adj_P = Male_Adj_P,
    Male_celladjustLM_Delta = deltaBM,
    Female_celladjustLM_Adj_P = Female_Adj_P,
    Female_celladjustLM_Delta = deltaBF
  )
sa_overlap_significant <- sa %>%
  full_join (sa_celladj, by = c("Row.names")) %>%
  filter(if_any(ends_with("Adj_P"), ~ .x < 0.05))

sig_autosomes <- autosome_results_allinfo$Whole_Adj_P < 0.05
hyper_autosomes <- autosome_results_allinfo$deltaBA > 0.05 & autosome_results_allinfo$Whole_Adj_P < 0.05
hypo_autosomes <- autosome_results_allinfo$deltaBA < -0.05 & autosome_results_allinfo$Whole_Adj_P < 0.05
sum(sig_autosomes) #127 sig
sum(hyper_autosomes) #8 Hyper
sum(hypo_autosomes) #39 Hypo
sig_autosomes_celladjust <- autosome_celladjusted_results_allinfo$Whole_Adj_P < 0.05
hyper_autosomes_celladjust <- autosome_celladjusted_results_allinfo$deltaBA > 0.05 & autosome_celladjusted_results_allinfo$Whole_Adj_P < 0.05
hypo_autosomes_celladjust <- autosome_celladjusted_results_allinfo$deltaBA < -0.05 & autosome_celladjusted_results_allinfo$Whole_Adj_P < 0.05
sum(sig_autosomes_celladjust) #138 sig
sum(hyper_autosomes_celladjust) #3 Hyper
sum(hypo_autosomes_celladjust) #31 Hypo 
sig_autosome_overlap <- sig_autosomes & sig_autosomes_celladjust #90 sig overlap 
sig_hyper_autosome_overlap <- hyper_autosomes & hyper_autosomes_celladjust #2 overlap 
sig_hypo_autosome_overlap <- hypo_autosomes & hypo_autosomes_celladjust #28 overlap 
sig_autosomesM <- autosome_results_allinfo$Male_Adj_P < 0.05
hyper_autosomesM <- autosome_results_allinfo$deltaBM > 0.05 & autosome_results_allinfo$Male_Adj_P < 0.05
hypo_autosomesM <- autosome_results_allinfo$deltaBM < -0.05 & autosome_results_allinfo$Male_Adj_P < 0.05
sum(sig_autosomesM) #9 sig
sum(hyper_autosomesM) #1 Hyper
sum(hypo_autosomesM) #2 Hypo
sig_autosomes_celladjustM <- autosome_celladjusted_results_allinfo$Male_Adj_P < 0.05
hyper_autosomes_celladjustM <- autosome_celladjusted_results_allinfo$deltaBM > 0.05 & autosome_celladjusted_results_allinfo$Male_Adj_P < 0.05
hypo_autosomes_celladjustM <- autosome_celladjusted_results_allinfo$deltaBM < -0.05 & autosome_celladjusted_results_allinfo$Male_Adj_P < 0.05
sum(sig_autosomes_celladjustM) #18 sig
sum(hyper_autosomes_celladjustM) #2 Hyper
sum(hypo_autosomes_celladjustM) #3 Hypo 
sum(sig_autosomesM & sig_autosomes_celladjustM) #8 sig overlap 
sum(hyper_autosomesM & hyper_autosomes_celladjustM) #1 overlap 
sum(hypo_autosomesM & hypo_autosomes_celladjustM) #2 overlap
sig_autosomesF <- autosome_results_allinfo$Female_Adj_P < 0.05
hyper_autosomesF <- autosome_results_allinfo$deltaBF > 0.05 & autosome_results_allinfo$Female_Adj_P < 0.05
hypo_autosomesF <- autosome_results_allinfo$deltaBF < -0.05 & autosome_results_allinfo$Female_Adj_P < 0.05
sum(sig_autosomesF) #0 sig
sum(hyper_autosomesF) #0 Hyper
sum(hypo_autosomesF) #0 Hypo
sig_autosomes_celladjustF <- autosome_celladjusted_results_allinfo$Female_Adj_P < 0.05
hyper_autosomes_celladjustF <- autosome_celladjusted_results_allinfo$deltaBF > 0.05 & autosome_celladjusted_results_allinfo$Female_Adj_P < 0.05
hypo_autosomes_celladjustF <- autosome_celladjusted_results_allinfo$deltaBF < -0.05 & autosome_celladjusted_results_allinfo$Female_Adj_P < 0.05
sum(sig_autosomes_celladjustF) #2 sig
sum(hyper_autosomes_celladjustF) #0 Hyper
sum(hypo_autosomes_celladjustF) #0 Hypo 
sum(sig_autosomesF & sig_autosomes_celladjustF) #0 sig overlap 
sum(hyper_autosomesF & hyper_autosomes_celladjustF) #0 overlap 
sum(hypo_autosomesF & hypo_autosomes_celladjustF) #0 Overlap

#Cell decon probes 
celldeconprobes <- celldecon_probe_results %>%
  filter(if_any(matches("Adj_P$|Adj_P_celladjust$"), ~ .x < 0.05))

sa <- autosome_results_allinfo %>%
  select(
    Row.names,
    Closest_TSS_gene_name,
    Whole_LM_Adj_P = Whole_Adj_P,
    Whole_LM_Delta = deltaBA,
    Male_LM_Adj_P = Male_Adj_P,
    Male_LM_Delta = deltaBM,
    Female_LM_Adj_P = Female_Adj_P,
    Female_LM_Delta = deltaBF
  )

supplementary1_Whole <- autosome_results_allinfo %>%
  filter(Whole_Adj_P < 0.05) %>%
  select(
    Row.names,
    Closest_TSS_gene_name,
    Whole_Adj_P = Whole_Adj_P,
    Whole_Delta_beta = deltaBA,
    Beta_Threshold = All_diffmethylation
  )
write.csv(supplementary1_Whole, "./supplementary1_Whole.csv")

supplementary1_Male <- autosome_results_allinfo %>%
  filter(Male_Adj_P < 0.05) %>%
  select(
    Row.names,
    Closest_TSS_gene_name,
    Male_Adj_P = Male_Adj_P,
    Male_Delta_beta = deltaBM,
    Beta_Threshold = M_diffmethylation
  )
write.csv(supplementary1_Male, "./supplementary1_Male.csv")

supplementary2_Whole <- sa_overlap_significant %>% 
  filter(Whole_LM_Adj_P <0.05 | Whole_celladjustLM_Adj_P < 0.05) %>%
  select(
    Row.names,
    Closest_TSS_gene_name,
    Whole_Original_Model_Adj_P = Whole_LM_Adj_P,
    Whole_CellAdj_Model_Adj_P = Whole_celladjustLM_Adj_P,
    Whole_Delta_beta = Whole_LM_Delta
  )
supplementary2_Whole$Overlap <- "Not Overlap" 
supplementary2_Whole$Overlap[supplementary2_Whole$Whole_Original_Model_Adj_P < 0.05 & supplementary2_Whole$Whole_CellAdj_Model_Adj_P > 0.05] <- "Original Model Only"
supplementary2_Whole$Overlap[supplementary2_Whole$Whole_Original_Model_Adj_P > 0.05 & supplementary2_Whole$Whole_CellAdj_Model_Adj_P < 0.05] <- "Cell Adjusted Model Only" 
supplementary2_Whole$Overlap[supplementary2_Whole$Whole_Original_Model_Adj_P < 0.05 & supplementary2_Whole$Whole_CellAdj_Model_Adj_P < 0.05] <- "Overlap in Both Models" 
supplementary2_Whole$All_diffmethylation <- "Not Biologically Significant"
supplementary2_Whole$All_diffmethylation[supplementary2_Whole$Whole_Delta_beta > 0.05 & (supplementary2_Whole$Whole_Original_Model_Adj_P <0.05 | supplementary2_Whole$Whole_CellAdj_Model_Adj_P <0.05)] <- "Hyper-methylated"
supplementary2_Whole$All_diffmethylation[supplementary2_Whole$Whole_Delta_beta < -0.05 & (supplementary2_Whole$Whole_Original_Model_Adj_P <0.05 | supplementary2_Whole$Whole_CellAdj_Model_Adj_P <0.05)] <- "Hypo-methylated"

write.csv(supplementary2_Whole, "./supplementary2_Whole.csv")

supplementary2_Male <- sa_overlap_significant %>% 
  filter(Male_LM_Adj_P <0.05 | Male_celladjustLM_Adj_P < 0.05) %>%
  select(
    Row.names,
    Closest_TSS_gene_name,
    Male_Original_Model_Adj_P = Male_LM_Adj_P,
    Male_CellAdj_Model_Adj_P = Male_celladjustLM_Adj_P,
    Male_Delta_beta = Male_LM_Delta
  )
supplementary2_Male$Overlap <- "Not Overlap" 
supplementary2_Male$Overlap[supplementary2_Male$Male_Original_Model_Adj_P < 0.05 & supplementary2_Male$Male_CellAdj_Model_Adj_P > 0.05] <- "Original Model Only"
supplementary2_Male$Overlap[supplementary2_Male$Male_Original_Model_Adj_P > 0.05 & supplementary2_Male$Male_CellAdj_Model_Adj_P < 0.05] <- "Cell Adjusted Model Only" 
supplementary2_Male$Overlap[supplementary2_Male$Male_Original_Model_Adj_P < 0.05 & supplementary2_Male$Male_CellAdj_Model_Adj_P < 0.05] <- "Overlap in Both Models" 
supplementary2_Male$All_diffmethylation <- "Not Biologically Significant"
supplementary2_Male$All_diffmethylation[supplementary2_Male$Male_Delta_beta > 0.05 & (supplementary2_Male$Male_Original_Model_Adj_P <0.05 | supplementary2_Male$Male_CellAdj_Model_Adj_P <0.05)] <- "Hyper-methylated"
supplementary2_Male$All_diffmethylation[supplementary2_Male$Male_Delta_beta < -0.05 & (supplementary2_Male$Male_Original_Model_Adj_P <0.05 | supplementary2_Male$Male_CellAdj_Model_Adj_P <0.05)] <- "Hypo-methylated"

write.csv(supplementary2_Male, "./supplementary2_Male.csv")

supplementary2_Female <- sa_overlap_significant %>%
  filter(Female_LM_Adj_P <0.05 | Female_celladjustLM_Adj_P < 0.05) %>%
  select(
    Row.names,
    Closest_TSS_gene_name,
    Female_Original_Model_Adj_P = Female_LM_Adj_P,
    Female_CellAdj_Model_Adj_P = Female_celladjustLM_Adj_P,
    Female_Delta_beta = Female_LM_Delta
  )
supplementary2_Female$Overlap <- "Not Overlap" 
supplementary2_Female$Overlap[supplementary2_Female$Female_Original_Model_Adj_P < 0.05 & supplementary2_Female$Female_CellAdj_Model_Adj_P > 0.05] <- "Original Model Only"
supplementary2_Female$Overlap[supplementary2_Female$Female_Original_Model_Adj_P > 0.05 & supplementary2_Female$Female_CellAdj_Model_Adj_P < 0.05] <- "Cell Adjusted Model Only" 
supplementary2_Female$Overlap[supplementary2_Female$Female_Original_Model_Adj_P < 0.05 & supplementary2_Female$Female_CellAdj_Model_Adj_P < 0.05] <- "Overlap in Both Models" 
supplementary2_Female$All_diffmethylation <- "Not Biologically Significant"
supplementary2_Female$All_diffmethylation[supplementary2_Female$Female_Delta_beta > 0.05 & (supplementary2_Female$Female_Original_Model_Adj_P <0.05 | supplementary2_Female$Female_CellAdj_Model_Adj_P <0.05)] <- "Hyper-methylated"
supplementary2_Female$All_diffmethylation[supplementary2_Female$Female_Delta_beta < -0.05 & (supplementary2_Female$Female_Original_Model_Adj_P <0.05 | supplementary2_Female$Female_CellAdj_Model_Adj_P <0.05)] <- "Hypo-methylated"
write.csv(supplementary2_Female, "./supplementary2_Female.csv")


#Overlap with Auvinen 
#Whole 
whole_overlap_Auvinen <- supplementary2_Whole %>%
  left_join(
    Auvinen2024_sigprobes %>% select(ARTProbe, ARTdeltaB),     # select only needed columns
    by = c("Row.names" = "ARTProbe")) %>% # match Row.names to ARTProbe
  left_join(
    Auvinen2024_sigprobes %>% select(IVFProbe, IVFdeltaB),
    by = c("Row.names" = "IVFProbe")) 
write.csv(whole_overlap_Auvinen, "./whole.overlap_Auvinen.csv")
