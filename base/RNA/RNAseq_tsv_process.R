
setwd("~/data/MouseDNAmAtlas/ENCODE/adult")
library(data.table)
files<-list.dirs(recursive = FALSE)
fpkm_ad.m <- matrix(NA,nrow = 81881, ncol = length(files))

for (i in 1:length(files)){
  print(i)
  # 设置i为工作目录
  setwd(files[i])
  file_names <- list.files(pattern = "\\.tsv$", full.names = TRUE)
  biore <- matrix(NA,nrow = 81881,ncol = 2)
  for (j in 1:length(file_names)){
    data <- fread(file_names[j])
    fpkm <- data$FPKM
    gene <- data$gene_id
    biore[,j] <- fpkm
  }
  fpkm_ad.m[,i] <- rowMeans(biore)
  setwd("..")
}
tissue <- unlist(strsplit(basename(files), "\\."))
colnames(fpkm_ad.m) <- tissue
rownames(fpkm_ad.m) <- gene




setwd("~/data/MouseDNAmAtlas/ENCODE/adult")
tpm_ad.m <- matrix(NA,nrow = 81881, ncol = length(files))

for (i in 1:length(files)){
  print(i)
  # 设置i为工作目录
  setwd(files[i])
  file_names <- list.files(pattern = "\\.tsv$", full.names = TRUE)
  biore <- matrix(NA,nrow = 81881,ncol = 2)
  for (j in 1:length(file_names)){
    data <- fread(file_names[j])
    tpm <- data$TPM
    gene <- data$gene_id
    biore[,j] <- tpm
  }
  tpm_ad.m[,i] <- rowMeans(biore)
  setwd("..")
}
tissue <- unlist(strsplit(basename(files), "\\."))
colnames(tpm_ad.m) <- tissue
rownames(tpm_ad.m) <- gene

setwd("~/data/MouseDNAmAtlas/ENCODE/Embryonic/RNA/E_10")
library(data.table)
files<-list.dirs(recursive = FALSE)
fpkm_10.m <- matrix(NA,nrow = 81881, ncol = length(files))

for (i in 1:length(files)){
  print(i)
  # 设置i为工作目录
  setwd(files[i])
  file_names <- list.files(pattern = "\\.tsv$", full.names = TRUE)
  biore <- matrix(NA,nrow = 81881,ncol = 2)
  for (j in 1:length(file_names)){
    data <- fread(file_names[j])
    fpkm <- data$FPKM
    gene <- data$gene_id
    biore[,j] <- fpkm
  }
  fpkm_10.m[,i] <- rowMeans(biore)
  setwd("..")
}
tissue <- unlist(strsplit(basename(files), "\\."))
colnames(fpkm_10.m) <- tissue
rownames(fpkm_10.m) <- gene




setwd("~/data/MouseDNAmAtlas/ENCODE/Embryonic/RNA/E_10")
tpm_10.m <- matrix(NA,nrow = 81881, ncol = length(files))

for (i in 1:length(files)){
  print(i)
  # 设置i为工作目录
  setwd(files[i])
  file_names <- list.files(pattern = "\\.tsv$", full.names = TRUE)
  biore <- matrix(NA,nrow = 81881,ncol = 2)
  for (j in 1:length(file_names)){
    data <- fread(file_names[j])
    tpm <- data$TPM
    gene <- data$gene_id
    biore[,j] <- tpm
  }
  tpm_10.m[,i] <- rowMeans(biore)
  setwd("..")
}
tissue <- unlist(strsplit(basename(files), "\\."))
colnames(tpm_10.m) <- tissue
rownames(tpm_10.m) <- gene

setwd("~/data/MouseDNAmAtlas/ENCODE/Embryonic/RNA/E_16")
library(data.table)
files<-list.dirs(recursive = FALSE)
fpkm_16.m <- matrix(NA,nrow = 81881, ncol = length(files))

for (i in 1:length(files)){
  print(i)
  # 设置i为工作目录
  setwd(files[i])
  file_names <- list.files(pattern = "\\.tsv$", full.names = TRUE)
  biore <- matrix(NA,nrow = 81881,ncol = 2)
  for (j in 1:length(file_names)){
    data <- fread(file_names[j])
    fpkm <- data$FPKM
    gene <- data$gene_id
    biore[,j] <- fpkm
  }
  fpkm_16.m[,i] <- rowMeans(biore)
  setwd("..")
}
tissue <- unlist(strsplit(basename(files), "\\."))
colnames(fpkm_16.m) <- tissue
rownames(fpkm_16.m) <- gene




setwd("~/data/MouseDNAmAtlas/ENCODE/Embryonic/RNA/E_16")
tpm_16.m <- matrix(NA,nrow = 81881, ncol = length(files))

for (i in 1:length(files)){
  print(i)
  # 设置i为工作目录
  setwd(files[i])
  file_names <- list.files(pattern = "\\.tsv$", full.names = TRUE)
  biore <- matrix(NA,nrow = 81881,ncol = 2)
  for (j in 1:length(file_names)){
    data <- fread(file_names[j])
    tpm <- data$TPM
    gene <- data$gene_id
    biore[,j] <- tpm
  }
  tpm_16.m[,i] <- rowMeans(biore)
  setwd("..")
}
tissue <- unlist(strsplit(basename(files), "\\."))
colnames(tpm_16.m) <- tissue
rownames(tpm_16.m) <- gene

colnames(fpkm_ad.m) <-c("ENCSR000CGT","ENCSR000CGX","ENCSR137GMB","ENCSR000CGZ","ENCSR000CHA","ENCSR784VGW","ENCSR000CHB","ENCSR000CHC","ENCSR000CHL","ENCSR000CHI","ENCSR000CGW","ENCSR000CHK")
colnames(tpm_ad.m) <-c("ENCSR000CGT","ENCSR000CGX","ENCSR137GMB","ENCSR000CGZ","ENCSR000CHA","ENCSR784VGW","ENCSR000CHB","ENCSR000CHC","ENCSR000CHL","ENCSR000CHI","ENCSR000CGW","ENCSR000CHK")

fpkm.m <- cbind(fpkm_10.m,fpkm_16.m)
fpkm.m <- cbind(fpkm.m,fpkm_ad.m)

tpm.m <- cbind(tpm_10.m,tpm_16.m)
tpm.m <- cbind(tpm.m,tpm_ad.m)

idx1<- which(rowSums(fpkm.m) == 0)
idx2<- which(rowSums(tpm.m) == 0)

fpkm.m <- fpkm.m[-idx1,]
tpm.m <- tpm.m[-idx2,]

idx<-match(setdiff(rownames(tpm.m),rownames(fpkm.m)),rownames(tpm.m))

tpm.m <- tpm.m[-idx,]


setwd("~/data/MouseDNAmAtlas/GSE42836")
files<-list.dirs(recursive = FALSE)
tpm_geo.m <- matrix(NA,nrow = 81881, ncol = length(files))

for (i in 1:length(files)){
  print(i)
  # 设置i为工作目录
  setwd(files[i])
  file_names <- list.files(pattern = "\\.tsv$", full.names = TRUE)
  biore <- matrix(NA,nrow = 81881,ncol = 2)
  for (j in 1:length(file_names)){
    data <- fread(file_names[j])
    tpm <- data$TPM
    gene <- data$gene_id
    biore[,j] <- tpm
  }
  tpm_geo.m[,i] <- rowMeans(biore)
  setwd("..")
}
tissue <- unlist(strsplit(basename(files), "\\."))
colnames(tpm_geo.m) <- tissue
rownames(tpm_geo.m) <- gene

idx<- which(rowSums(tpm_geo.m) == 0)
tpm_geo.m <- tpm_geo.m[-idx,]
tpm_geo.m <- log2(tpm_geo.m+1)
