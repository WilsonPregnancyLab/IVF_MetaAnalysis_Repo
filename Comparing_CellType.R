library(ggplot2) #Version 3.5.1
library(dplyr) #Version 1.1.4
library(tidyr) #Version 1.3.1
library(lme4) #Version 1.1.35.4
library(lmerTest) #Version 3.1.3
library(tidyr) #Version 1.3.1

estF_conceptandsex <- estF %>%
  mutate(concep_type = case_when( #Creating the concep_type variable column
    group == "IVFF" | group == "IVFM" ~ "IVF",
    group == "SPONTF" | group == "SPONTM" ~ "SPONT",
    .default = NA)) %>%
  mutate(sex = case_when( #Creating the sex variable column
    group == "IVFF" | group == "SPONTF" ~ "F",
    group == "IVFM" | group == "SPONTM" ~ "M",
    .default = NA))
estF_conceptandsex$concep_type <- as.factor(estF_conceptandsex$concep_type)
estF_conceptandsex$sex <- as.factor(estF_conceptandsex$sex)
estF_conceptandsex$Sample_Name <- rownames(estF_conceptandsex)
estF_dat <- merge(estF_conceptandsex, pheno, by = "Sample_Name")
estF_dat$sex <- as.factor(estF_dat$sex)
cell_types <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC", "Syncytiotrophoblast")

#Linear Model Comparing Cell Proportions Between Groups 
#Males: Comparing IVF sample cell proportions to Spontaneous sample cell proportions
male_cell_results <- data.frame(
  Cell_Type = character(),
  Beta_ConcepType = numeric(),
  P_ConcepType = numeric(),
  stringsAsFactors = FALSE
)
male_data <- subset(estF_dat, sex == "M")
for(cell_type in cell_types) {
  cat("Analyzing cell types:", cell_type, "\n")
  
  #Fit linear model (comparing conception type)
  formula <- as.formula(paste(cell_type, "~ concep_type + (1 | Study)"))
  lmm <- lmer(formula, data = male_data)
  
  #model summary 
  model_summary <- summary(lmm)
  
  #Extract fixed effects 
  fixed_effects <-  model_summary$coefficients 
  
  #extract results for conception type 
  beta_concep_type <- fixed_effects["concep_typeSPONT", "Estimate"]
  p_concep_type <- fixed_effects["concep_typeSPONT", "Pr(>|t|)"]
  
  # Append results to the data frame
  male_cell_results <- rbind(male_cell_results, data.frame(
    Cell_Type = cell_type,
    Beta_ConcepType = beta_concep_type,
    P_ConcepType = p_concep_type
  ))
}
write.csv(male_cell_results, "male_stratified_celltype_results.csv", row.names = FALSE)
#Females: Comparing IVF sample cell proportions to Spontaneous sample cell proportions
female_cell_results <- data.frame(
  Cell_Type = character(),
  Beta_ConcepType = numeric(),
  P_ConcepType = numeric(),
  stringsAsFactors = FALSE
)
female_data <- subset(estF_dat, sex == "F")
for(cell_type in cell_types) {
  cat("Analyzing cell types:", cell_type, "\n")
  
  #Fit linear model (comparing conception type)
  formula <- as.formula(paste(cell_type, "~ concep_type + (1 | Study)"))
  lmm <- lmer(formula, data = female_data)
  
  #model summary 
  model_summary <- summary(lmm)
  
  #Extract fixed effects 
  fixed_effects <-  model_summary$coefficients 
  
  #extract results for conception type 
  beta_concep_type <- fixed_effects["concep_typeSPONT", "Estimate"]
  p_concep_type <- fixed_effects["concep_typeSPONT", "Pr(>|t|)"]
  
  # Append results to the data frame
  female_cell_results <- rbind(female_cell_results, data.frame(
    Cell_Type = cell_type,
    Beta_ConcepType = beta_concep_type,
    P_ConcepType = p_concep_type
  ))
}
write.csv(female_cell_results, "female_stratified_celltype_results.csv", row.names = FALSE)
#IVF only Samples: Comparing Male and Female Samples
IVF_cell_results <- data.frame(
  Cell_Type = character(),
  Beta_Sex = numeric(),
  P_Sex = numeric(),
  stringsAsFactors = FALSE
)
IVF_data <- subset(estF_dat, concep_type == "IVF")
for(cell_type in cell_types) {
  cat("Analyzing cell types:", cell_type, "\n")
  
  #Fit linear model (comparing conception type)
  formula <- as.formula(paste(cell_type, "~ sex + (1 | Study)"))
  lmm <- lmer(formula, data = IVF_data)
  
  #model summary 
  model_summary <- summary(lmm)
  
  #Extract fixed effects 
  fixed_effects <-  model_summary$coefficients 
  
  #extract results for conception type 
  beta_sex <- fixed_effects["sexM", "Estimate"]
  p_sex <- fixed_effects["sexM", "Pr(>|t|)"]
  
  # Append results to the data frame
  IVF_cell_results <- rbind(IVF_cell_results, data.frame(
    Cell_Type = cell_type,
    Beta_Sex = beta_sex,
    P_Sex = p_sex
  ))
}
write.csv(IVF_cell_results, "IVF_stratified_celltype_results.csv", row.names = FALSE)
#Spontaneous only Samples: Comparing Male and Female Samples
SPONT_cell_results <- data.frame(
  Cell_Type = character(),
  Beta_Sex = numeric(),
  P_Sex = numeric(),
  stringsAsFactors = FALSE
)
SPONT_data <- subset(estF_dat, concep_type == "SPONT")
for(cell_type in cell_types) {
  cat("Analyzing cell types:", cell_type, "\n")
  
  #Fit linear model (comparing conception type)
  formula <- as.formula(paste(cell_type, "~ sex + (1 | Study)"))
  lmm <- lmer(formula, data = SPONT_data)
  
  #model summary 
  model_summary <- summary(lmm)
  
  #Extract fixed effects 
  fixed_effects <-  model_summary$coefficients 
  
  #extract results for conception type 
  beta_sex <- fixed_effects["sexM", "Estimate"]
  p_sex <- fixed_effects["sexM", "Pr(>|t|)"]
  
  # Append results to the data frame
  SPONT_cell_results <- rbind(SPONT_cell_results, data.frame(
    Cell_Type = cell_type,
    Beta_Sex = beta_sex,
    P_Sex = p_sex
  ))
}
write.csv(SPONT_cell_results, "SPONT_stratified_celltype_results.csv", row.names = FALSE)




#Mixed Effects Model - No stratification of sex or conception type done prior
# Initialize a results data frame
cell_unstratified_results <- data.frame(
  Cell_Type = character(),
  Beta_Sex = numeric(),
  P_Sex = numeric(),
  Beta_ConcepType = numeric(),
  P_ConcepType = numeric(),
  stringsAsFactors = FALSE
)
# Loop through each cell type
for (cell_type in cell_types) {
  cat("Analyzing cell type:", cell_type, "\n")
  
  # Fit the linear mixed-effects model
  formula <- as.formula(paste(cell_type, "~ sex + concep_type  + (1 | Study)"))
  lmm <- lmer(formula, data = estF_dat)
  
  # Summarize the model
  model_summary <- summary(lmm)
  
  # Extract fixed effects
  fixed_effects <- model_summary$coefficients
  
  # Extract results for sex and conception type
  beta_sex <- fixed_effects["sexM", "Estimate"]
  p_sex <- fixed_effects["sexM", "Pr(>|t|)"]
  
  beta_concep_type <- fixed_effects["concep_typeSPONT", "Estimate"]
  p_concep_type <- fixed_effects["concep_typeSPONT", "Pr(>|t|)"]
  
  # Append results to the data frame
  cell_unstratified_results <- rbind(cell_unstratified_results, data.frame(
    Cell_Type = cell_type,
    Beta_Sex = beta_sex,
    P_Sex = p_sex,
    Beta_ConcepType = beta_concep_type,
    P_ConcepType = p_concep_type
  ))
}
write.csv(cell_unstratified_results, "LMM_CellType_Results.csv", row.names = FALSE)


#Plotting 
#Reshape data into long format for ggplot2
#Plot separated by study 
long_data <- estF_dat %>%
  pivot_longer(
    cols = c(Trophoblasts, Stromal, Hofbauer, Endothelial, nRBC, Syncytiotrophoblast),
    names_to = "Cell_Type",
    values_to = "Proportion"
  )
png("celltype.png", width=4000, height=1500, res=200)
ggplot(long_data, aes(x = Cell_Type, y = Proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.1, size = 1) +
  facet_wrap(~ Study, scales = "free_y") +
  labs(
    title = "Boxplot of Cell Type Proportions by Group and Sex",
    x = "Group",
    y = "Proportion",
    fill = "Sex"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off()
#plot with all samples in one figure
png("celltype_oneplot.png", width = 1500, height = 1500, res = 200)
ggplot(long_data, aes(x = Cell_Type, y = Proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  labs(
    title = "Boxplot of Cell Type Proportions by Group",
    x = "Cell Type",
    y = "Proportion",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
dev.off()
