library(ggplot2)
library(ggsignif)
library(ggpubr)
library(tidyverse)
library(planet)
library(minfi)
library(EpiDISH)
data("plBetas")
data("plCellCpGsThird")

#Cell Deconvolution - to predict placental cell composition in each sample

epidish_RPC_IVFF <- epidish(
  beta.m = placmet_IVF_F_betasonly[rownames(plCellCpGsThird), ],
  ref.m = plCellCpGsThird, 
  method = "RPC")

epidish_RPC_IVFM <- epidish(
  beta.m = placmet_IVF_M_betasonly[rownames(plCellCpGsThird), ],
  ref.m = plCellCpGsThird, 
  method = "RPC")

epidish_RPC_SPONTF <- epidish(
  beta.m = placmet_SPONT_F_betasonly[rownames(plCellCpGsThird), ],
  ref.m = plCellCpGsThird, 
  method = "RPC")

epidish_RPC_SPONTM <- epidish(
  beta.m = placmet_SPONT_M_betasonly[rownames(plCellCpGsThird), ],
  ref.m = plCellCpGsThird, 
  method = "RPC")



estF_IVFF <- epidish_RPC_IVFF$estF %>% 
  as.data.frame() %>% mutate(algorithm = "RPC")

estF_IVFM <- epidish_RPC_IVFM$estF %>% 
  as.data.frame() %>% mutate(algorithm = "RPC")

estF_SPONTF <- epidish_RPC_SPONTF$estF %>% 
  as.data.frame() %>% mutate(algorithm = "RPC")

estF_SPONTM <- epidish_RPC_SPONTM$estF %>% 
  as.data.frame() %>% mutate(algorithm = "RPC")



estF_IVFF <- estF_IVFF[!rownames(estF_IVFF) %in% c('AvgBIVF'),]
estF_IVFM <- estF_IVFM[!rownames(estF_IVFM) %in% c('AvgBIVF'),]
estF_SPONTF <- estF_SPONTF[!rownames(estF_SPONTF) %in% c('AvgBSPONT'),]
estF_SPONTM <- estF_SPONTM[!rownames(estF_SPONTM) %in% c('AvgBSPONT'),]

estF_IVFF$group <- 'IVFF'
estF_IVFM$group <- 'IVFM'
estF_SPONTF$group <- 'SPONTF'
estF_SPONTM$group <- 'SPONTM'

estF <- rbind(estF_IVFF, estF_IVFM, estF_SPONTF, estF_SPONTM)


# ANOVA for each cell type 
Troph.aov <- aov(Trophoblasts ~ group, data = estF)
summary(Troph.aov)
strom.aov <- aov(Stromal ~ group, data = estF)
summary(strom.aov)
hof.aov <- aov(Hofbauer ~ group, data = estF)
summary(hof.aov)
endo.aov <- aov(Endothelial ~ group, data = estF)
summary(endo.aov)
nRBC.aov <- aov(nRBC ~ group, data = estF)
summary(nRBC.aov)
syn.aov <- aov(Syncytiotrophoblast ~ group, data = estF)
summary(syn.aov)

# Bonferroni post-hoc test
troph.bonf <- pairwise.t.test(estF$Trophoblasts, estF$group, p.adjust.method = 'bonferroni')
strom.bonf <- pairwise.t.test(estF$Stromal, estF$group, p.adjust.method = 'bonferroni')
hof.bonf <- pairwise.t.test(estF$Hofbauer, estF$group, p.adjust.method = 'bonferroni')
endo.bonf <- pairwise.t.test(estF$Endothelial, estF$group, p.adjust.method = 'bonferroni')
nRBC.bonf <- pairwise.t.test(estF$nRBC, estF$group, p.adjust.method = 'bonferroni')
syn.bonf <- pairwise.t.test(estF$Syncytiotrophoblast, estF$group, p.adjust.method = 'bonferroni')

# Plotting 
cell_table <- data.frame(type = rep(c('troph', 'strom', 'hof', 'endo', 'nRBC', 'syn'), each = 4),
                         u = c(mean(estF_IVFF$Trophoblasts), mean(estF_IVFM$Trophoblasts), mean(estF_SPONTF$Trophoblasts), mean(estF_SPONTM$Trophoblasts),
                               mean(estF_IVFF$Stromal), mean(estF_IVFM$Stromal), mean(estF_SPONTF$Stromal), mean(estF_SPONTM$Stromal),
                               mean(estF_IVFF$Hofbauer), mean(estF_IVFM$Hofbauer), mean(estF_SPONTF$Hofbauer), mean(estF_SPONTM$Hofbauer),
                               mean(estF_IVFF$Endothelial), mean(estF_IVFM$Endothelial), mean(estF_SPONTF$Endothelial), mean(estF_SPONTM$Endothelial),
                               mean(estF_IVFF$nRBC), mean(estF_IVFM$nRBC), mean(estF_SPONTF$nRBC), mean(estF_SPONTM$nRBC),
                               mean(estF_IVFF$Syncytiotrophoblast), mean(estF_IVFM$Syncytiotrophoblast), mean(estF_SPONTF$Syncytiotrophoblast), mean(estF_SPONTM$Syncytiotrophoblast)),
                         s = c(sd(estF_IVFF$Trophoblasts), sd(estF_IVFM$Trophoblasts), sd(estF_SPONTF$Trophoblasts), sd(estF_SPONTM$Trophoblasts),
                               sd(estF_IVFF$Stromal), sd(estF_IVFM$Stromal), sd(estF_SPONTF$Stromal), sd(estF_SPONTM$Stromal),
                               sd(estF_IVFF$Hofbauer), sd(estF_IVFM$Hofbauer), sd(estF_SPONTF$Hofbauer), sd(estF_SPONTM$Hofbauer),
                               sd(estF_IVFF$Endothelial), sd(estF_IVFM$Endothelial), sd(estF_SPONTF$Endothelial), sd(estF_SPONTM$Endothelial),
                               sd(estF_IVFF$nRBC), sd(estF_IVFM$nRBC), sd(estF_SPONTF$nRBC), sd(estF_SPONTM$nRBC),
                               sd(estF_IVFF$Syncytiotrophoblast), sd(estF_IVFM$Syncytiotrophoblast), sd(estF_SPONTF$Syncytiotrophoblast), sd(estF_SPONTM$Syncytiotrophoblast)),
                         group = rep(c('IVFF','IVFM','SPONTF','SPONTM'), 6)
)

png(filename = "./cell_decon_anova_adjFunnorm.png", height = 7.5, width = 10, units = "in", res = 750)
ggplot(cell_table, aes(fill = group, y = u, x = type)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  scale_fill_manual(values = c("#E67000", "#0060B5","#F7D56E", "#8EC9FF")) +
  geom_errorbar(aes(ymin = u-s, ymax = u+s), width = .2, position = position_dodge(0.9)) +
  theme_classic()
dev.off()
