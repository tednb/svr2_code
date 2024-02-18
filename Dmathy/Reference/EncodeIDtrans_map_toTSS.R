source("~/code/Dmathy/code/hts_map_to_tss_mm10.R")
source("~/code/Dmathy/code/ID_trans.R")
pheno <- read.csv("~/data/MouseDNAmAtlas/ENCODE/pheno.csv")
pheno$X <- NULL
# 1. id transform and location ---------------------------------------------
load("~/data/MouseDNAmAtlas/ENCODE/tpm.Rd")
result <- EnsemblID_to_symbol(tpm.m,'~/chain/gencode.vM21.basic.annotation.gtf.gz')
expr_tpm <- result$expr_tpm_averaged
gene_info <- result$gene_location_info
expr_tpm$ENCSR611PTP <- NULL
save(expr_tpm,gene_info,file = "tpm+location.Rdata")

# 3. map to tss200 --------------------------------------------------------
load("~/data/MouseDNAmAtlas/ENCODE/tpm+location.Rdata")
load("~/data/MouseDNAmAtlas/ENCODE/Encode_beta.Rdata")
load("~/data/MouseDNAmAtlas/ENCODE/Encode_depth.Rdata")
mat <- depth_mask(Encode.beta,Encode.depth,thre = 10)
#saveRDS(mat,file = "~/data/MouseDNAmAtlas/ENCODE/beta_10X.rds")
mat <- readRDS("~/data/MouseDNAmAtlas/ENCODE/beta_10X.rds")
mf_10depth <- map_to_tss(mat,gene_info,region = 200)
save(mf_10depth,file = "~/data/MouseDNAmAtlas/ENCODE/mf_tss200_10X.Rdata")

# cluster -----------------------------------------------------------------
load("~/data/MouseDNAmAtlas/ENCODE/mf_tss500_10X.Rdata")
library(broom)
library(isva)
wgbs <- gsub("_rep1|_rep2","",colnames(mf_10depth)[1:67])
tname <- pheno$tissue[match(wgbs,pheno$WGBS)]
t_name <- paste(tname,sapply(strsplit(colnames(mf_10depth)[1:67],"_"),function(x) x[2]),sep = "_")
t_name <- c(t_name,colnames(mf_10depth)[68:91])
m <- mf_10depth - rowMeans(mf_10depth)
n <- EstDimRMT(m)$dim
svd.o <- svd(m)
v_top <- svd.o$v[,1:n] #n=4
rownames(v_top) <- t_name
dist_matrix <- as.dist(1 - cor(t(v_top))) #相关性距离
hc <- hclust(dist_matrix, method = "complete")
stage <- c(rep("E10",27),rep("E16",40),rep("Adult",24))
idx_e10 <- which(stage %in% c("E10"))
idx_e16 <- which(stage %in% c("E16"))
idx_adu <- which(stage %in% c("Adult"))
colors <- c("E10" = "blue", "E16" = "red", "Adult" = "green")
# Apply colors to the stages vector
stage_colors <- colors[stage]
# Convert hclust object to dendrogram
dend <- as.dendrogram(hc)
labels_order <- order.dendrogram(dend)
ordered_colors <- stage_colors[labels_order]
# Color the labels of the dendrogram
labels_colors(dend) <- ordered_colors
# 假设dend是您已经创建好并想要调整显示的dendrogram对象
pdf("cluster.pdf", width = 14, height = 10)
# 获取当前y轴的范围
y_range <- range(rev(hc$height))
# 绘制树状图时扩大y轴显示范围
plot(dend, ylim=c(-1, y_range[2] * 1.5), main = "Methylation level on TSS500 (n = 46)", 
     xlab = paste0("Correlation distance,", n," PCs"))
legend("bottomleft", legend = names(colors), fill = colors, horiz = TRUE, cex = 0.8)
dev.off()



### distribution
pheno<-read.csv("~/data/MouseDNAmAtlas/ENCODE/pheno.csv")
pheno$X <- NULL
pheno <- pheno[,-c(47:52)]
load("~/data/MouseDNAmAtlas/ENCODE/mf_tss200_10dep.Rdata")
pdf("wgbs_gene_dtr_10dep.pdf",width = 6,height = 3)
for(i in 1:ncol(mf.m)){
  plot(density(mf.m[,i]),main = paste0(pheno$tissue[i]," ",pheno$stage[i]))
}
dev.off()

pdf("wgbs_dtr.pdf",width = 6,height = 3)
for(i in 1:ncol(mf_10depth)){
  plot(density(mf_10depth[,i]),main = colnames(mf_10depth)[i])
}
dev.off()

pdf("wgbs_dtr_10dep.pdf",width = 6,height = 3)
for(i in 1:ncol(m)){
  plot(density(m[,i]),main = paste0(pheno$tissue[i]," ",pheno$stage[i]))
}
dev.off()
