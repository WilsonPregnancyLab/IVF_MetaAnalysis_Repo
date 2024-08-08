#R Version - R4.4.1
#Packages
library(limma) #version 3.60.4 
library(lumi) #version 2.56.0 
library(stringr) #version 1.5.1
library(ggplot2) #version 3.5.1
library(gridExtra) #version 2.3
library(ggrepel) # version 0.9.5
library(rtracklayer) # version 1.64.0
library(GenomicRanges) #version 1.56.1
library(ggpubr) #version 0.6.0
library(Gviz) #version 1.48.0
library(png) #version 0.1.8
library("TxDb.Hsapiens.UCSC.hg19.knownGene") #version 3.2.2
library(EnsDb.Hsapiens.v75) #version 2.99.0
library(Homo.sapiens) #version 1.3.1

placmet_IVF_F_betasonly <- placmet_IVF_F_autosomes[, !(names(placmet_IVF_F_autosomes) == "ProbeIVF")]
placmet_IVF_F_betasonly <- placmet_IVF_F_betasonly[, !(names(placmet_IVF_F_betasonly) == "AvgBIVF")]

placmet_SPONT_F_betasonly <- placmet_SPONT_F_autosomes[, !(names(placmet_SPONT_F_autosomes) == "ProbeSPONT")]
placmet_SPONT_F_betasonly <- placmet_SPONT_F_betasonly[, !(names(placmet_SPONT_F_betasonly) == "AvgBSPONT")]

placmet_IVF_M_betasonly <- placmet_IVF_M_autosomes[, !(names(placmet_IVF_M_autosomes) == "ProbeIVF")]
placmet_IVF_M_betasonly <- placmet_IVF_M_betasonly[, !(names(placmet_IVF_M_betasonly) == "AvgBIVF")]

placmet_SPONT_M_betasonly <- placmet_SPONT_M_autosomes[, !(names(placmet_SPONT_M_autosomes) == "ProbeSPONT")]
placmet_SPONT_M_betasonly <- placmet_SPONT_M_betasonly[, !(names(placmet_SPONT_M_betasonly) == "AvgBSPONT")]

#Calculate Coefficients of Variance
IVF_F_CV <- placmet_IVF_F_betasonly
IVF_F_CV$CV_IVF_F <- apply(IVF_F_CV, 1, function(row) {
  cv <- (sd(row, na.rm=TRUE) / mean(row, na.rm = TRUE)) * 100
  return (cv) })

IVF_M_CV <- placmet_IVF_M_betasonly
IVF_M_CV$CV_IVF_M <- apply(IVF_M_CV, 1, function(row) {
  cv <- (sd(row, na.rm=TRUE) / mean(row, na.rm = TRUE)) * 100
  return (cv) })

SPONT_F_CV <- placmet_SPONT_F_betasonly
SPONT_F_CV$CV_SPONT_F <- apply(SPONT_F_CV, 1, function(row) {
  cv <- (sd(row, na.rm=TRUE) / mean(row, na.rm = TRUE)) * 100
  return (cv) })

SPONT_M_CV <- placmet_SPONT_M_betasonly
SPONT_M_CV$CV_SPONT_M <- apply(SPONT_M_CV, 1, function(row) {
  cv <- (sd(row, na.rm=TRUE) / mean(row, na.rm = TRUE)) * 100
  return (cv) })

#Coefficients of variance density plot
png(filename = "./variability_autosomes_adjFunnorm.png", height = 9, width = 9, units = "in", res = 300)
plot(density(SPONT_F_CV$CV_SPONT_F), col = "#F7D56E", lwd = 3,xlim = c(0,60), ylim = c(0,0.05), 
     xlab = "Coefficient of Variance (N = 305160 probes)", ylab = "Density", main = " ", cex.lab = 1.2) +
  lines(density(SPONT_M_CV$CV_SPONT_M), col = "#8EC9FF",lwd = 3) +
  lines(density(IVF_F_CV$CV_IVF_F), col = "#E67000", lwd = 3) +
  lines(density(IVF_M_CV$CV_IVF_M), col = "#0060B5", lwd = 3) 
legend("topright", legend = c("Spontaneous Female", "Spontaneous Male", "IVF Female", "IVF Male"),
       col = c("#F7D56E", "#8EC9FF", "#E67000", "#0060B5"),
       lwd = 3, bty = "n", cex = 1.2)
dev.off()

#comparing distributions 
ks.test(SPONT_F_CV$CV_SPONT_F, SPONT_M_CV$CV_SPONT_M)
ks.test(IVF_F_CV$CV_IVF_F, IVF_M_CV$CV_IVF_M)
ks.test(IVF_F_CV$CV_IVF_F, SPONT_F_CV$CV_SPONT_F)
ks.test(IVF_M_CV$CV_IVF_M, SPONT_M_CV$CV_SPONT_M)

#Creating table for genome browser 
IVF_F_CV$probes <- rownames(IVF_F_CV)
SPONT_F_CV$probes <- rownames(SPONT_F_CV)
IVF_M_CV$probes <- rownames(IVF_M_CV)
SPONT_M_CV$probes <- rownames(SPONT_M_CV)

IVF_table <- merge(merge(merge(merge(merge(IVF_F_CV[,c("CV_IVF_F", "probes")], SPONT_F_CV[, c("CV_SPONT_F", "probes")], by = "probes"),
                                      IVF_M_CV[, c("CV_IVF_M", "probes")], by = "probes"),
                                SPONT_M_CV[,c("CV_SPONT_M", "probes")], by = "probes"),
                          placmet_MaleAvgbetas_autosomes[,c("deltaB", "probes")], by = "probes"),
                    placmet_FemaleAvgbetas_autosomes[,c("deltaB", "probes")], by = "probes")
names(IVF_table)[names(IVF_table) == 'deltaB.x'] <- 'deltaB_F'
names(IVF_table)[names(IVF_table) == 'deltaB.y'] <- 'deltaB_M'
IVF_table$delta_CoV_M <- IVF_table$CV_SPONT_M - IVF_table$CV_IVF_M
IVF_table$delta_CoV_F <- IVF_table$CV_SPONT_F - IVF_table$CV_IVF_F

#Genome Browser - showing coefficients of variance, deltaB, and deltaCoV
IVF_table$X <- NULL
hm450.annotation <- read.table('HM450.hg19.manifest.tsv.gz', sep = '\t', header = T) #downloaded from https://zwdzwd.github.io/InfiniumAnnotation
cpgIslands <- read.table('cpgIslandExthg19.txt.gz', sep = '\t', header = F) #cpgIslandExt.txt.gz downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/
TSS.annotation <- read.table('Price_anno_450K.tsv', sep = '\t', header = T) #downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42409

cpgIslands_rel <- cpgIslands[,c('V2','V3','V4')]
names(cpgIslands_rel) <- c('chr','start','end')

IVF_table <- merge(IVF_table, hm450.annotation[,c('CpG_chrm','CpG_beg', 'Probe_ID')], by.x = 'probes', by.y = 'Probe_ID')
names(IVF_table) <- c('Probe','CV_IVF_F','CV_IVF_M', 'CV_SPONT_F', 'CV_SPONT_M', 'deltaB_F', 'deltaB_M', 'delta_CoV_F', 'delta_CoV_M', 'chr', 'C_pos')


IVF_F <- IVF_table[,c('chr', 'C_pos', 'CV_IVF_F')]
IVF_M <- IVF_table[,c('chr', 'C_pos', 'CV_IVF_M')]
SPONT_F <- IVF_table[,c('chr', 'C_pos', 'CV_SPONT_F')]
SPONT_M <- IVF_table[,c('chr', 'C_pos', 'CV_SPONT_M')]

names(IVF_F) <- c('chr','C_pos','CV')
names(IVF_M) <- c('chr','C_pos','CV')
names(SPONT_F) <- c('chr','C_pos','CV')
names(SPONT_M) <- c('chr','C_pos','CV')


dB_F <- IVF_table[,c('chr', 'C_pos', 'deltaB_F')]
dB_M <- IVF_table[,c('chr', 'C_pos', 'deltaB_M')]

names(dB_F) <- c('chr','C_pos','dB')
names(dB_M) <- c('chr','C_pos','dB')


dCoV_F <- IVF_table[,c('chr', 'C_pos', 'delta_CoV_F')]
dCoV_M <- IVF_table[,c('chr', 'C_pos', 'delta_CoV_M')]

names(dCoV_F) <- c('chr','C_pos','dCV')
names(dCoV_M) <- c('chr','C_pos','dCV')


lg_df <- data.frame(x = 1:7, y = 0, color = factor(c('IVF Female CoV','IVF Male CoV','Spontaneous Female CoV','Spontaneous Male CoV', '', 'Female delta B/CoV','Male delta B/CoV'), levels = c('IVF Female CoV','IVF Male CoV','Spontaneous Female CoV','Spontaneous Male CoV', '', 'Female delta B/CoV','Male delta B/CoV')))
legend <- ggplot(lg_df, aes(x = x, y = y, color = color)) +
  geom_line() +
  scale_color_manual(values = c('#E67000','#0060B5','#F7D56E','#8EC9FF', '#FFFFFF', '#9400CD','#228B22')) +
  labs(color='Group') +
  theme_classic() +
  theme(text = element_text(size = 12),
        legend.text=element_text(size=12)
  )
lg <- as_ggplot(cowplot::get_legend(legend))
png(file = 'lg_temp.png', width = 3, height = 3, unit = 'in', res = 300)
lg
dev.off()

get_gb_plot_panel <- function(chrom_num = 0, start = 0, stop = Inf, probe, pm, plot_height = 10, justSide = 'left'){
  
  if(!missing(probe)){
    chrom <- IVF_table[IVF_table$Probe == probe,]$chr
    chrom_num <- substr(chrom, 4, nchar(chrom))
    pos <- IVF_table[IVF_table$Probe == probe,]$C_pos
    start <- pos - pm
    stop <- pos + pm
    closest_TSS <- TSS.annotation[TSS.annotation$ID == probe,]$Closest_TSS_gene_name
  }
  
  
  else{
    chrom = paste0('chr', chrom_num)
  }
  
  IVF_F_sub <- subset(IVF_F, chr == chrom)
  IVF_M_sub <- subset(IVF_M, chr == chrom)
  SPONT_F_sub <- subset(SPONT_F, chr == chrom)
  SPONT_M_sub <- subset(SPONT_M, chr == chrom)
  
  dB_F_sub <- subset(dB_F, chr == chrom)
  dB_M_sub <- subset(dB_M, chr == chrom)
  
  dCoV_F_sub <- subset(dCoV_F, chr == chrom)
  dCoV_M_sub <- subset(dCoV_M, chr == chrom)
  
  cpgIslands_sub <- subset(cpgIslands_rel, chr == chrom)
  
  gene <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, chromosome = chrom, geneSymbol = TRUE, start = start, end = stop, transcriptAnnotation = "symbol", name = "Gene")
  gene_out <<- gene
  IVF_F_Grange <- makeGRangesFromDataFrame(IVF_F_sub, start.field = 'C_pos', end.field = 'C_pos')
  mcols(IVF_F_Grange)$value <- IVF_F_sub$CV
  
  IVF_M_Grange <- makeGRangesFromDataFrame(IVF_M_sub, start.field = 'C_pos', end.field = 'C_pos')
  mcols(IVF_M_Grange)$value <- IVF_M_sub$CV
  
  SPONT_F_Grange <- makeGRangesFromDataFrame(SPONT_F_sub, start.field = 'C_pos', end.field = 'C_pos')
  mcols(SPONT_F_Grange)$value <- SPONT_F_sub$CV
  
  SPONT_M_Grange <- makeGRangesFromDataFrame(SPONT_M_sub, start.field = 'C_pos', end.field = 'C_pos')
  mcols(SPONT_M_Grange)$value <- SPONT_M_sub$CV
  
  
  dB_F_Grange <- makeGRangesFromDataFrame(dB_F_sub, start.field = 'C_pos', end.field = 'C_pos')
  mcols(dB_F_Grange)$value <- dB_F_sub$dB
  
  dB_M_Grange <- makeGRangesFromDataFrame(dB_M_sub, start.field = 'C_pos', end.field = 'C_pos')
  mcols(dB_M_Grange)$value <- dB_M_sub$dB
  
  
  dCoV_F_Grange <- makeGRangesFromDataFrame(dCoV_F_sub, start.field = 'C_pos', end.field = 'C_pos')
  mcols(dCoV_F_Grange)$value <- dCoV_F_sub$dCV
  
  dCoV_M_Grange <- makeGRangesFromDataFrame(dCoV_M_sub, start.field = 'C_pos', end.field = 'C_pos')
  mcols(dCoV_M_Grange)$value <- dCoV_M_sub$dCV
  
  cpgIslands_gr <- makeGRangesFromDataFrame(cpgIslands_sub, start.field = 'start', end.field = 'end')
  
  ids <- ranges(gene)
  
  if(length(ids) > 0){
    mappings <- select(Homo.sapiens, keys(TxDb.Hsapiens.UCSC.hg19.knownGene, 'TXNAME'), 'SYMBOL', 'TXNAME')
    rownames(mappings) <- mappings$TXNAME
    ids$symbol <- mappings[ids$symbol,]$SYMBOL
  }
  
  ranges(gene) <- ids
  gene <- gene[!is.na(ranges(gene)$symbol),]
  
  
  if(!missing(probe)){
    marker <- makeGRangesListFromDataFrame(
      data.frame(
        'chr' = chrom,
        'start' = pos - pm/200,
        'stop' = pos + pm/200
      )
    )
    mtrack <- AnnotationTrack(marker, showAxis = F, name = 'Probe', col = 'black', fill = 'black')
  }
  
  displayPars(gene) <- list(fontsize.group = 20, fontcolor.group = 'black', fontfamily.group = 'arial', fontface.group = 'plain', background.title = 'brown', col.line = '#555555', col = '#555555', lwd.group = 0.1, just.group = justSide)
  gtrack <- GenomeAxisTrack(name = paste0('Chr. ', chrom_num), col = 'black', fontcolor = 'black', showTitle = T)
  CVtrack <- OverlayTrack(trackList = list(DataTrack(IVF_F_Grange, col = '#E67000', name = 'Coefficient of Variance (CoV)', type = 'a'), DataTrack(IVF_M_Grange, col = '#0060B5', type = 'a'), DataTrack(SPONT_F_Grange, col = '#F7D56E', type = 'a'), DataTrack(SPONT_M_Grange, col = '#8EC9FF', type = 'a')), legend = T)
  dBtrack <- OverlayTrack(trackList = list(DataTrack(dB_F_Grange, col = '#9400CD', name = 'delta B', type = 'a'), DataTrack(dB_M_Grange, col = '#228B22', type = 'a')), legend = T)
  dCVtrack <- OverlayTrack(trackList = list(DataTrack(dCoV_F_Grange, col = '#9400CD', name = 'delta CoV', type = 'a'), DataTrack(dCoV_M_Grange, col = '#228B22', type = 'a')), legend = T)
  
  cpgIslandTrack <- AnnotationTrack(cpgIslands_gr, name = 'CpG', col = '#444444', fill = '#444444')
  
  png(file = paste0(probe, '_tmp.png'), width = 7, height = plot_height, unit = 'in', res = 300)
  
  if(!missing(probe)){            
    #track height calculations
    CVtrack_h = 0.25*(10 - 10/plot_height)/10
    dBtrack_h = 0.15*(10 - 10/plot_height)/10
    dCVtrack_h = 0.15*(10 - 10/plot_height)/10
    mtrack_h = 0.05*(10 - 10/plot_height)/10
    gtrack_h = 0.1*(10 - 10/plot_height)/10
    cpgIslandTrack_h = 0.05*(10 - 10/plot_height)/10
    gene_h = 0.25*(10 - 10/plot_height)/10
    plotTracks(list(CVtrack, dCVtrack, dBtrack, mtrack, gtrack,cpgIslandTrack, gene), sizes = c(CVtrack_h, dCVtrack_h, dBtrack_h, mtrack_h, gtrack_h, cpgIslandTrack_h, gene_h), title.width = 1, fontsize = 12, cex.axis = 1, background.title = 'brown', lwd = 2, main = paste0('Probe ID: ',probe, "\nClosest TSS: ", closest_TSS), fontface.main = 'plain', cex.main = 1.5, from = start, to = stop)
  }
  else{
    #track height calculations
    CVtrack_h = 0.25*(10 - 10/plot_height)/10
    dBtrack_h = 0.15*(10 - 10/plot_height)/10
    dCVtrack_h = 0.15*(10 - 10/plot_height)/10
    gtrack_h = 0.1*(10 - 10/plot_height)/10
    gene_h = 0.3*(10 - 10/plot_height)/10
    cpgIslandTrack_h = 0.05*(10 - 10/plot_height)/10
    plotTracks(list(CVtrack, dCVtrack, dBtrack, gtrack, cpgIslandTrack, gene), sizes = c(CVtrack_h, dCVtrack_h, dBtrack_h, gtrack_h,cpgIslandTrack_H, gene_h), title.width = 1, fontsize = 12, cex.axis = 1, background.title = 'brown', lwd = 2, fontface.main = 'plain', cex.main = 1.5, from = start, to = stop)
  }
  
  dev.off()
  
}

dev.off()
get_gb_plot_panel('1', 0, Inf, probe = 'cg10048212', pm = 60000, plot_height = 14, justSide = 'right')
dev.off()
get_gb_plot_panel('1', 0, Inf, probe = 'cg25822783', pm = 50000, plot_height = 14)
dev.off()
get_gb_plot_panel('1', 0, Inf, probe = 'cg26186239', pm = 30000, plot_height = 14, justSide = 'right')
dev.off()


cg10048212_panel <- ggplot() + annotation_custom(grid::rasterGrob(readPNG('cg10048212_tmp.png'),
                                                                  width=unit(1,"npc"),
                                                                  height=unit(1,"npc")),
                                                 -Inf, Inf, -Inf, Inf) + theme_void()

cg25822783_panel <- ggplot() + annotation_custom(grid::rasterGrob(readPNG('cg25822783_tmp.png'),
                                                                  width=unit(1,"npc"),
                                                                  height=unit(1,"npc")),
                                                 -Inf, Inf, -Inf, Inf) + theme_void()

cg26186239_panel <- ggplot() + annotation_custom(grid::rasterGrob(readPNG('cg26186239_tmp.png'),
                                                                  width=unit(1,"npc"),
                                                                  height=unit(1,"npc")),
                                                 -Inf, Inf, -Inf, Inf) + theme_void()

legend_panel <- ggplot() + annotation_custom(grid::rasterGrob(readPNG('lg_temp.png'),
                                                              width=unit(1,"npc"),
                                                              height=unit(1,"npc")),
                                             -Inf, Inf, -Inf, Inf) + theme_void() + coord_fixed()

plots <- cowplot::plot_grid(cg10048212_panel,cg25822783_panel, cg26186239_panel, legend_panel, ncol=4, align = "h", rel_widths = c(7, 7, 7, 3))


png(file = 'variability_panels.png', width = 24, height = 14, unit = 'in', res = 300)
plots
dev.off()
