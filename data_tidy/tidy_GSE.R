setwd("Renv")

library(GEOquery)
library(pbapply)
library(data.table) #fread
library(tidyverse)

gse_ano <- getGEO("GSE180474", GSEMatrix= TRUE, destdir=".",
                  AnnotGPL = F,
                  getGPL = F)
phenotypes.df <- pData(gse_ano[[1]])
gse <- read.table(file = "/mnt/local-disk/data/guoxiaolong/data/GSE213478/GSE213478_methylation_DNAm_noob_final_BMIQ_all_tissues_987.txt")
gse_matrix <- gse$V2
gse <- read_csv("/mnt/local-disk/data/guoxiaolong/data/eGTEX/GSE180474_MatrixProcessed.csv")

test <- pblapply(gse_matrix, function(x) {strsplit(x[[1]], ",")})
test[[1]] <- NULL
test1 <- pbsapply(test, function(x) {as.numeric(x[[1]])})
Allbeta.m <- t(test1[-1,]) #beta value
rm(test,test1)
# CpGs
gse_cpg <- gse$V1
test <- pblapply(gse_cpg, function(x) {strsplit(x[[1]], ",")})
test[[1]] <- NULL
rownames(Allbeta.m) <- unlist(test)
colnames(Allbeta.m) <- rownames(phenotypes.df)


#age summary
boxplot(age ~ tissue, data = pheno.df, 
        main = "Age Distribution by Org Type", 
        xlab = "Organization Type", ylab = "Age")
