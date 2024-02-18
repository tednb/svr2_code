library(parallel)
library(data.table)
library(progress)
setwd("~/data/RRBS/mouse/GSE134238")
files <- list.files(pattern = "\\.cov$")
# genes' promoter position
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
# 24528 genes
promoter<-promoters(genes(txdb), upstream = 500,downstream = 0)
start <- promoter@ranges@start
gene_id <- promoter$gene_id
# build 200 range list
ranges.m <- cbind(start, start + 199)
rownames(ranges.m) <- gene_id
sample_name<-sapply(files,function(x) {regmatches(x, regexpr("GSM\\d{7}", x))})
mf.m <- matrix(data = NA,ncol = length(files),nrow = length(gene_id))
rownames(mf.m) <- gene_id
colnames(mf.m) <- sample_name

calc_mean_fraction <- function(range,raw.m) {
  judge <- raw.m$mm10 >= range[1] & raw.m$mm10 <= range[2]
  mean(raw.m$fraction[judge],na.rm = TRUE)
}
pb <- progress_bar$new(total = length(files))
for (i in seq_along(files)) {
  pb$tick()
  raw.m <- fread(files[[i]])
  raw.m$V3 <- NULL
  setnames(raw.m, c("Chr","mm10","fraction","m","um"))
  idx <- which((raw.m$m + raw.m$um) > 5)
  raw.m <- raw.m[idx,]
  # Calculate mean fraction in parallel
  result <- mclapply(1:nrow(ranges.m), function(x) calc_mean_fraction(ranges.m[x,], raw.m), mc.cores = 100)
  mf.m[,i] <- unlist(result)
}
mf.m[is.nan(mf.m)] <- NA
colSums(!is.na(mf.m))
colSums(!is.na(mf.m)) / nrow(ranges.m)
