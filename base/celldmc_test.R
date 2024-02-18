setwd("Renv")

library(EpiDISH)
library(dplyr)
library(qvalue)
library(ggplot2)
library(pbapply)
library(broom)
library(tidyverse)
library(magick)
library(gridExtra)
library(ggplotify)
library(gplots)
library(grid)
library(png)

load("~/data/MRC-NSHD/mergedBUC.RData") #new gold-standard smk-DMCs
load("~/data/MRC-NSHD/PhenoTypesBUC.RData") #buccal swab samples
PhenoTypesMRG.lv$sample_name <- paste(PhenoTypesMRG.lv$SentrixPos, PhenoTypesMRG.lv$SentrixID, sep = "_")
colnames(mergedBUC.m) <- PhenoTypesMRG.lv$sample_name
smk <- factor(PhenoTypesMRG.lv$SmokStat53)
#Epi,Lym and Mye
buccalFrac.m <- hepidish(beta.m = mergedBUC.m, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m, h.CT.idx = 3, method = 'RPC')
boxplot(buccalFrac.m)
buccalFrac.m <- subset(buccalFrac.m , select = -c(Fib, CD8T, Eosino))
sum_cols <- apply(buccalFrac.m, 1 , sum)
buccalFrac.m <- as.data.frame(apply(buccalFrac.m, 2,function(x) {x/sum_cols}))
buccalFrac.m$Lym <- buccalFrac.m[,"NK"] + buccalFrac.m[,"B"] + buccalFrac.m[,"CD4T"]
buccalFrac.m$Mye <- buccalFrac.m[,"Mono"] + buccalFrac.m[,"Neutro"]
generalFrac.m <- as.matrix(subset(buccalFrac.m, select = c(Epi, Lym, Mye)))
generalFrac.m <- subset(buccalFrac.m, select = c(Epi, Lym, Mye))

#Epi and IC
buccalw2.m <- subset(epidish(beta.m = mergedBUC.m, ref.m = centEpiFibIC.m, method = 'RPC')$estF, select = -c(Fib))
sum_cols <- apply(buccalw2.m, 1 , sum)
buccalw2.m <- as.data.frame(apply(buccalw2.m, 2,function(x) {x/sum_cols}))
#cellDMC
celldmc.o <- CellDMC(mergedBUC.m, smk ,as.matrix(generalFrac.m), #3
                     adjPMethod = "fdr",
                     adjPThresh = 0.05,
                     cov.mod = NULL,
                     sort = FALSE,
                     mc.cores = 40)
dmcts<- celldmc.o$dmct
coefs <- celldmc.o$coe

celldmcW2.o <- CellDMC(mergedBUC.m, smk ,as.matrix(buccalw2.m), #2
                     adjPMethod = "fdr",
                     adjPThresh = 0.05,
                     cov.mod = NULL,
                     sort = FALSE,
                     mc.cores = 40)
dmctw2<- celldmcW2.o$dmct
coefw2<- celldmcW2.o$coe


#t-stat
all_t <- as.data.frame(lapply(coefs, function(x) x$t))
rownames(all_t) <- rownames(coefs$Epi)




#FDR threshold
fdr_threshold_epi <- min(abs(coefs$Epi$t[sapply(coefs$Epi$adjP, function(x) {x<=0.05})]))





#scatter plot
points_df <- all_t[rownames(gsSMKjh.m), ]
#points_df<- points_df[which(abs(points_df$Epi) >= fdr_threshold_epi), ]
scatterplot <- function(x){
  print(x)
pdf(paste0(paste0("Epi_vs",x),".pdf"), height = 5, width = 5)
smoothScatter(all_t[, "Epi"], all_t[, x], xlab = "t -stat. (CellDMC:Epi)", 
              ylab = paste("t -stat. (CellDMC:", x, ")", sep = ""), 
              xlim = c(-7,7), ylim = c(-7,7), cex = 1.5, 
              main = "Smoking DMCTs predicted by CellDMC")
abline(v = 0, lty = 2, col = "black")
abline(h = 0, lty = 2, col = "black")
abline(v = fdr_threshold_epi, lty = 2, col = "blue")
abline(v = -fdr_threshold_epi, lty = 2, col = "blue")
abline(h = fdr_threshold_oth, lty = 2, col = "blue")
abline(h = -fdr_threshold_oth, lty = 2, col = "blue")
#points_df<- points_df[which(abs(points_df$Epi) >= fdr_threshold_epi | abs(points_df[[x]]) >= fdr_threshold_oth), ]
#points(points_df$Epi, points_df[[x]], col = "red", pch = 15, cex = 0.3)
#legend("bottomleft", cex = 0.5,legend = "gold standard smk-DMCs", pch = 15, col = "red")
dev.off()
}
n <- colnames(dmcts)[3:4]
p <- list()
for (i in 1:length(n)){
  fdr_threshold_oth <- min(abs(coefs[[i+1]]$t[sapply(coefs[[i+1]]$adjP, function(m) {m<=0.05})]))
scatterplot(n[i])
}
# merge pdfs
pdf_dir <- "/mnt/local-disk/data/guoxiaolong/Renv"
pdf_files <- list.files(pdf_dir, pattern = "*.pdf", full.names = TRUE)
png_files <- list()
# convert each PDF file to a PNG file and add it to the list
for (i in 1:length(pdf_files)) {
  pdf_path <- pdf_files[i]
  png_path <- gsub(".pdf", ".png", pdf_path)
  system(paste("convert -density 300", shQuote(pdf_path), shQuote(png_path)))
  png_files[[length(png_files) + 1]] <- png_path
}
# create a list to store the grid objects
grid_list <- list()
# import each PNG file as a grid object and add it to the list
for (i in 1:length(png_files)) {
  grid_list[[i]] <- rasterGrob(readPNG(png_files[[i]]))
}
# arrange the grid objects in a grid
grid_arrange = function(grid_list, nrow, ncol) {
  grid_matrix = matrix(list(), nrow=nrow, ncol=ncol)
  for (i in 1:length(grid_list)) {
    row = ((i-1) %/% ncol) + 1
    col = ((i-1) %% ncol) + 1
    grid_matrix[[row, col]] = grid_list[[i]]
  }
  do.call("grid.arrange", c(grid_matrix, ncol=ncol, nrow=nrow))
}
pdf("t-stat.pdf", width = 11.69, height = 8.27, paper = "a4r")
grid_arrange(grid_list, nrow=1, ncol=2)
dev.off()

# heatmap
dmcts<-as.data.frame(dmcts)
DMCs_df <- dmcts[rownames(gsSMKjh.m),]
DMCs_smk <- t(as.matrix(DMCs_df[DMCs_df$DMC == 1 & !is.na(DMCs_df$DMC),][-1]))
DMCs <- t(as.matrix(dmcts[dmcts$DMC == 1,][-1]))
heatmap_f <- function(x){
dev.new()
my_colors <- c("#FFC0CB", "#E5E5E1", "#B0E0E6")
names(my_colors) <- c("Hypo", "No DM", "Hyper")
par(cex.main = 0.75)
heatmap.2(
  x, 
  Rowv = FALSE, 
  Colv = FALSE, 
  trace = "none",
  dendrogram = "none",
  col = my_colors, 
  cexRow = 0.6, 
  key = FALSE, 
  keysize = 1.5,
  margins = c(6,4),
  xlab = "CpG",
  ylab = "Cell type",
  key.ylab = NULL,
  main = "CellDMC predictions of gold-standard smk-DMCs in blood"
)
legend(
  "bottomleft",
  legend = c("Hypo", "No DM", "Hyper"),
  fill = my_colors,
  bty = "o",
  cex = 0.6,
  pt.cex = 0.5
)

}
heatmap_f(DMCs_smk)
dev.off()



































