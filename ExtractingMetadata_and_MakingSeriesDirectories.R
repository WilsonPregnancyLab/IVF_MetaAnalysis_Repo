#R Version - R4.4.1
#Packages
library(GEOquery) #version 2.72.0

#Extracting Metadata from GEO
#GSE208529 
gse208529 <- getGEO('GSE208529', GSEMatrix = TRUE)
## Preview metadata fields
titles_gse208529 <- pData(phenoData(gse208529[[1]]))[1:3,]
## Selected metadata
metadata208529 <- pData(phenoData(gse208529[[1]]))[,c("title","geo_accession","platform_id","gender:ch1","group:ch1","sample type:ch1","ethnicity:ch1")]
write.csv(metadata208529, "./metadata208529.csv")

#GSE120250
gse120250 <- getGEO('GSE120250', GSEMatrix = TRUE) 
titles_gse120250 <- pData(phenoData(gse120250[[1]]))[1:3,]
metadata120250<-pData(phenoData(gse120250[[1]]))[,c("title","geo_accession","platform_id","gender:ch1","art treatment:ch1","source_name_ch1","outlier:ch1")]
write.csv(metadata120250, "./metadata120250.csv")

#GSE75248
gse75248 <- getGEO('GSE75248', GSEMatrix = TRUE) 
titles_gse75248 <- pData(phenoData(gse75248[[1]]))[1:3,]
metadata75248<-pData(phenoData(gse75248[[1]]))[,c("title","geo_accession","platform_id","characteristics_ch1.1","characteristics_ch1.2","characteristics_ch1.3","characteristics_ch1.4","characteristics_ch1.5","characteristics_ch1.6","characteristics_ch1.7","characteristics_ch1.8","characteristics_ch1.9","characteristics_ch1.10","characteristics_ch1.11","characteristics_ch1.11","characteristics_ch1.12","characteristics_ch1.13","characteristics_ch1.14","characteristics_ch1.15","characteristics_ch1.16","characteristics_ch1.17","characteristics_ch1.18","characteristics_ch1.19","characteristics_ch1.20","characteristics_ch1.21","characteristics_ch1.22","characteristics_ch1.23","characteristics_ch1.24","characteristics_ch1.25","source_name_ch1")]
write.csv(metadata75248, "./metadata75248.csv")


#Downloading IDAT Files from GEO 
##Go to GEO, select and download .tar of all raw IDAT files
##Transfer .tar file from local computer to server if needed 
##Extract individual files in a folder that will become the base directory (all IDATS + StudySampleSheet.csv) - one for each GSE study 
