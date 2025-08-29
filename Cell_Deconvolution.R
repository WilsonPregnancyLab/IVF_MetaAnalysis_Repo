library(tidyverse)
library(planet) #Version 1.12.0
library(minfi) #version 1.50.0
library(EpiDISH) #Version 2.20.0
data("plBetas")
data("plCellCpGsThird")

#Cell Deconvolution - to predict placental cell composition in each sample
placmet_IVF_M_betasonly <- as.data.frame(placmet_filtfun_M[!rownames(placmet_filtfun_M) %in% c(chrXprobes$probeID,chrYprobes$probeID), IVFmetadata_M$Sample_Name]) 
placmet_SPONT_M_betasonly <- as.data.frame(placmet_filtfun_M[!rownames(placmet_filtfun_M) %in% c(chrXprobes$probeID,chrYprobes$probeID), Spontaneousmetadata_M$Sample_Name]) 
placmet_IVF_F_betasonly <- as.data.frame(placmet_filtfun_F[!rownames(placmet_filtfun_F) %in% c(chrXprobes$probeID,chrYprobes$probeID), IVFmetadata_F$Sample_Name]) 
placmet_SPONT_F_betasonly <- as.data.frame(placmet_filtfun_F[!rownames(placmet_filtfun_F) %in% c(chrXprobes$probeID,chrYprobes$probeID), Spontaneousmetadata_F$Sample_Name])

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

#Estimates of cell type proportion breakdown for each sample and information if they are from IVF Females, IVF Males, Spontaneous Females or Spontaneous Males
estF <- rbind(estF_IVFF, estF_IVFM, estF_SPONTF, estF_SPONTM)
