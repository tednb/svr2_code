library(data.table)
library(pbapply)
library(parallel)
library(progress)
# 196 age dtr
pdf("age.pdf",width = 5,height = 5)
ggplot(pheno, aes(x = age, fill = tissue)) +
  geom_bar(position = "stack") +
  labs(x = "Age", y = "Count", fill = "Tissue", title = "RRBS data (n = 196)") +
  scale_y_continuous(breaks = seq(0, 70, by = 10)) +
  theme_bw()
dev.off()
########################### start
load("~/data/scage/REF/sample_info.Rd")
folder_path <- "~/data/scage/REF/data/All_21103/"
file_names <- list.files(folder_path, pattern = "\\.txt$", full.names = TRUE)
files <- list.files(folder_path, pattern = "\\.txt$", full.names = F)
pattern <- "GSM(\\d+)"
match <- regexpr(pattern, files)
sample_name<-regmatches(files, match)
cpgs <- list()

for (i in seq_along(file_names)) {
  print(i)
  raw.m <- fread(file_names[[i]])
  setnames(raw.m, c("Chr","cg","mm10","m","allm","fraction","interval"))
  cpgs[[i]] <- raw.m$`mm10`
}
cpg <- Reduce(union, cpgs)
m_all <- matrix(data = NA,ncol = length(file_names),nrow = length(cpg))
colnames(m_all) <- sample_name
rownames(m_all) <- cpg
for (i in seq_along(file_names)) {
  print(i)
  raw.m <- fread(file_names[[i]])
  setnames(raw.m, c("Chr","cg","mm10","m","allm","fraction","interval"))
  cpg_s <- raw.m$`mm10`
  idx_na <- which(raw.m$`allm` < 10)
  raw.m$`fraction`[idx_na] <- NA
  idx <- match(cpg_s,cpg)
  m_all[idx,i] <- raw.m$`fraction`
}
idx <- match(rownames(pheno), colnames(m_all))
m_all <- m_all[,idx]
# 90%
judge <- function(x) {
  sum(!is.na(x))/length(x) >= 0.9
}
idx_ok <- unlist(mclapply(1:nrow(m_all), function(i) judge(m_all[i, ]),mc.cores = 100))
m_all <- m_all[idx_ok,] # 983555 * 196

library(impute)
m_all<- impute.knn(m_all,k=5)$data

save(m_all,file = "merge_10COV.Rd")
source("~/code/Dmathy/code/Damthy.R")
qc.o <- new("qc.o")

qc.o@m <- m_all
qc.o@s <- pheno

qc.o@pca.o <- PCA(qc.o,2)
#load("~/data/scage/REF/sample_info.Rd")

# mapping
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm10)

#######################INPUT GENES OF INTEREST####
load("~/result/REFERENCES/ref_LIVER/mouse_ref.Rd")
genes = rownames(mice_ref)


#######################GET TSS#############################
#get transcriptional start sites for all genes of interest

ensembl <- useMart('ensembl',dataset = "mmusculus_gene_ensembl",host = "https://asia.ensembl.org/")
tss <- getBM(attributes=c('mgi_symbol','start_position', 'end_position','strand'), filters = c("mgi_symbol"), values=genes, mart=ensembl)
save(tss,file = "marker_codin.Rd")
########################## MAP ##############################
load("~/Renv/ref/marker_codin.Rd")
mf.m <- matrix(data = NA,ncol = ncol(m_all),nrow = nrow(tss))
colnames(mf.m) <- colnames(m_all)
rownames(mf.m) <- tss$mgi_symbol


judge <- function(i,cpg){
  if(tss[i,4] == 1){
    idx <- which(cpg<=tss[i,2] & cpg >= tss[i,2]-200)
  }else{
    idx <- which(cpg>=tss[i,3] & cpg <= tss[i,3]+200)
  }
  if(length(idx) > 1){
    mf.m[i,] <- colMeans(m_all[idx,])
  }else if(length(idx) == 1){
    mf.m[i,] <- m_all[idx,]
  }else{
    mf.m[i,] <-NA
  }
}

cpg <- as.numeric(rownames(m_all))
mf.m <- do.call(rbind,mclapply(1:nrow(tss),function(i) {judge(i,cpg)},mc.cores = 100))
rownames(mf.m) <- tss$mgi_symbol
idx <- which(is.na(rowSums(mf.m)))
data.m <- mf.m[-idx,]
save(data.m,file = "map_cov10_200_gse60012.Rd")
qc.o <- new("qc.o")

qc.o@m <- data.m
qc.o@s <- pheno

qc.o@pca.o <- PCA(qc.o,2)
# check
load("h_markers.Rd")
hep_df <- homologene(hep_DMC,inTax = 9606, outTax = 10090)
other_df <- homologene(rest_DMC,inTax = 9606, outTax = 10090)

length(intersect(rownames(data.m),hep_df$`10090`))
length(intersect(rownames(data.m),other_df$`10090`))

# data

data_liver.m <- data.m[,pheno$tissue == "Liver"]
data_lung.m <- data.m[,pheno$tissue == "Lung"]
idx_other <- setdiff(1:nrow(pheno),which(pheno$tissue == "Liver"))
data_other.m <- data.m[,idx_other]
# ref

idx <- match(rownames(data.m),rownames(mice_ref))
ref.m <- mice_ref[idx,-3]
identical(rownames(ref.m),rownames(data.m))
#################################Convolution#########
RPC(ref.m,data_liver.m) #3.806
RPC(ref.m,data_other.m) #3.884
#QP(ref.m,data_liver.m)

hepf.lv <- list()
for ( i in levels(factor(pheno$tissue))) {
  idx <- which(pheno$tissue == i)
  m <- data.m[,idx]
  hepf.lv[[i]] <- QP(ref.m,m)
}
pdf("hep_f_qp_500.pdf",width = 10,height = 5)
par(mar=c(8,4,4,4))
boxplot(hepf.lv,ylab = "Fraction",las = 3,main = "Hepatocyte")
abline(h=0)
dev.off()


QP(ref.m,data_liver.m) # 4.565
QP(ref.m,data_other.m) # 4.732


###
qc.o <- new("qc.o")

qc.o@m <- data.m
qc.o@s <- pheno

qc.o@pca.o <- PCA(qc.o,2)
