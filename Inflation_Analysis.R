library(qqman) #Version 0.1.9

#Globlal inflation test 
#whole all autosomes
whole_autosome_RawP <- as.numeric(autosome_results$Whole_Raw_P)
chisq_whole_autosomeprobes <- qchisq(1 - whole_autosome_RawP, df = 1)
lambda_whole_autosomeprobes <- median(chisq_whole_autosomeprobes) / qchisq(0.5, df = 1)
cat("Inflation lambda:", lambda_whole_autosomeprobes, "\n") #Inflation lambda: 1.34531
#male all autosomes 
male_autosome_RawP <- as.numeric(autosome_results$Male_Raw_P)
chisq_male_autosomeprobes <- qchisq(1 - male_autosome_RawP, df = 1)
lambda_male_autosomeprobes <- median(chisq_male_autosomeprobes) / qchisq(0.5, df = 1)
cat("Inflation lambda:", lambda_male_autosomeprobes, "\n") #Inflation lambda: 0.9397216
#female all autosomes
female_autosome_RawP <- as.numeric(autosome_results$Female_Raw_P)
chisq_female_autosomeprobes <- qchisq(1 - female_autosome_RawP, df = 1)
lambda_female_autosomeprobes <- median(chisq_female_autosomeprobes) / qchisq(0.5, df = 1)
cat("Inflation lambda:", lambda_female_autosomeprobes, "\n") #Inflation lambda: 1.44585
#whole all autosomes cell adjusted model 
whole_autosome_celladjusted_RawP <- as.numeric(autosome_celladjusted_results$Whole_Raw_P)
chisq_whole_autosome_celladjustedprobes <- qchisq(1 - whole_autosome_celladjusted_RawP, df = 1)
lambda_whole_autosome_celladjustedprobes <- median(chisq_whole_autosome_celladjustedprobes) / qchisq(0.5, df = 1)
cat("Inflation lambda:", lambda_whole_autosome_celladjustedprobes, "\n") #Inflation lambda: 0.9880773
#male all autosomes cell adjusted model
male_autosome_celladjusted_RawP <- as.numeric(autosome_celladjusted_results$Male_Raw_P)
chisq_male_autosome_celladjustedprobes <- qchisq(1 - male_autosome_celladjusted_RawP, df = 1)
lambda_male_autosome_celladjustedprobes <- median(chisq_male_autosome_celladjustedprobes) / qchisq(0.5, df = 1)
cat("Inflation lambda:", lambda_male_autosome_celladjustedprobes, "\n") #Inflation lambda: 0.9739169
#female all autosomes cell adjusted model 
female_autosome_celladjusted_RawP <- as.numeric(autosome_celladjusted_results$Female_Raw_P)
chisq_female_autosome_celladjustedprobes <- qchisq(1 - female_autosome_celladjusted_RawP, df = 1)
lambda_female_autosome_celladjustedprobes <- median(chisq_female_autosome_celladjustedprobes) / qchisq(0.5, df = 1)
cat("Inflation lambda:", lambda_female_autosome_celladjustedprobes, "\n") #Inflation lambda: 1.25226

#Plotting 
#All Autosomes 
png("./autosome_qq_panel.png", height = 5, width = 15, units = "in", res = 300)
par(mfrow = c(1, 3))
qq(whole_autosome_RawP, main = "Whole Autosome", ylim = c(0, 12))
qq(female_autosome_RawP, main = "Female Autosome", ylim = c(0, 12))
qq(male_autosome_RawP, main = "Male Autosome", ylim = c(0, 12))
dev.off()
#All Autosomes cell adjusted model 
png("./autosome_celladjusted_qq_panel.png", height = 5, width = 15, units = "in", res = 300)
par(mfrow = c(1, 3))
qq(whole_autosome_celladjusted_RawP, main = "Whole Autosome Cell Adjusted", ylim = c(0, 12))
qq(female_autosome_celladjusted_RawP, main = "Female Autosome Cell Adjusted", ylim = c(0, 12))
qq(male_autosome_celladjusted_RawP, main = "Male Autosome Cell Adjusted", ylim = c(0, 12))
dev.off()


#Cell decon probes 
#whole cell decon probes
celldecon_whole_autosome_RawP <- as.numeric(celldecon_probe_results$Whole_Raw_P)
chisq_whole_celldecon <- qchisq(1 - celldecon_whole_autosome_RawP, df = 1)
lambda_whole_celldecon <- median(chisq_whole_celldecon) / qchisq(0.5, df = 1)
cat("Inflation lambda:", lambda_whole_celldecon, "\n") #Inflation lambda: 1.98726
#male cell decon probes 
celldecon_male_autosome_RawP <- as.numeric(celldecon_probe_results$Male_Raw_P)
chisq_male_celldecon <- qchisq(1 - celldecon_male_autosome_RawP, df = 1)
lambda_male_celldecon <- median(chisq_male_celldecon) / qchisq(0.5, df = 1)
cat("Inflation lambda:", lambda_male_celldecon, "\n") #Inflation lambda: 0.9939228
#female cell decon probes
celldecon_female_autosome_RawP <- as.numeric(celldecon_probe_results$Female_Raw_P)
chisq_female_celldecon <- qchisq(1 - celldecon_female_autosome_RawP, df = 1)
lambda_female_celldecon <- median(chisq_female_celldecon) / qchisq(0.5, df = 1)
cat("Inflation lambda:", lambda_female_celldecon, "\n") #Inflation lambda: 2.016082
#whole cell decon probes cell adjusted model 
celldecon_whole_autosomes_celladjusted_RawP <- as.numeric(celldecon_probe_results$Whole_Raw_P_celladjust)
chisq_whole_celldecon_celladjustedprobes <- qchisq(1 - celldecon_whole_autosomes_celladjusted_RawP, df = 1)
lambda_whole_celldecon_celladjustedprobes <- median(chisq_whole_celldecon_celladjustedprobes) / qchisq(0.5, df = 1)
cat("Inflation lambda:", lambda_whole_celldecon_celladjustedprobes, "\n") #Inflation lambda: 0.77914 
#male cell decon probes cell adjusted model
celldecon_male_autosomes_celladjusted_RawP <- as.numeric(celldecon_probe_results$Male_Raw_P_celladjust)
chisq_male_celldecon_celladjustedprobes <- qchisq(1 - celldecon_male_autosomes_celladjusted_RawP, df = 1)
lambda_male_celldecon_celladjustedprobes <- median(chisq_male_celldecon_celladjustedprobes) / qchisq(0.5, df = 1)
cat("Inflation lambda:", lambda_male_celldecon_celladjustedprobes, "\n")  #Inflation lambda: 1.044324
#female cell decon probes cell adjusted model 
celldecon_female_autosomes_celladjusted_RawP <- as.numeric(celldecon_probe_results$Female_Raw_P_celladjust)
chisq_female_celldecon_celladjustedprobes <- qchisq(1 - celldecon_female_autosomes_celladjusted_RawP, df = 1)
lambda_female_celldecon_celladjustedprobes <- median(chisq_female_celldecon_celladjustedprobes) / qchisq(0.5, df = 1)
cat("Inflation lambda:", lambda_female_celldecon_celladjustedprobes, "\n")  #Inflation lambda: 1.335826

#Plotting 
#Cell decon probes  
png("./celldecon_qq_panel.png", height = 5, width = 15, units = "in", res = 300)
par(mfrow = c(1, 3))
qq(celldecon_whole_autosome_RawP, main = "Whole Cell Decon Probes", ylim = c(0, 8))
qq(celldecon_female_autosome_RawP, main = "Female Cell Decon Probes", ylim = c(0, 8))
qq(celldecon_male_autosome_RawP, main = "Male Cell Decon Probes", ylim = c(0, 8))
dev.off()
#Cell decon probes cell adjusted model 
png("./celldecon_celladjusted_qq_panel.png", height = 5, width = 15, units = "in", res = 300)
par(mfrow = c(1, 3))
qq(celldecon_whole_autosomes_celladjusted_RawP, main = "Whole Cell Decon Probes Cell Adjusted", ylim = c(0, 8))
qq(celldecon_female_autosomes_celladjusted_RawP, main = "Female Cell Decon Probes Cell Adjusted", ylim = c(0, 8))
qq(celldecon_male_autosomes_celladjusted_RawP, main = "Male Cell Decon Probes Cell Adjusted", ylim = c(0, 8))
dev.off()
