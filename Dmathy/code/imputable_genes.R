setwd("~/data/MouseDNAmAtlas/ENCODE")
load("tpm+location.Rdata")
load("mf_tss500_10X.Rdata")
pheno <- read.csv("pheno.csv")
pheno$X <- NULL

# density plot for each sample --------------------------------------------
library(ggplot2)
mf.df <- as.data.frame(mf_10depth)
# 对于数据框中的每一列
pdf("wgbs_density.pdf",width = 4,height = 4)
for (col_name in colnames(mf.df)) {
  # 打印当前列名
  print(col_name)
  
  # 绘制当前列的密度图
  
  p <- ggplot(mf.df, aes_string(x=col_name)) + 
    geom_density() +
    ggtitle(paste("Density plot for", col_name)) +
    theme_minimal()
  
  print(p)
}
dev.off()

# average_replicates ------------------------------------------------------

colnames <- names(mf.df)
# 使用正则表达式提取前缀，即 "_rep" 之前的部分
prefixes <- sapply(strsplit(colnames, "_rep"), `[`, 1)
# 创建一个新的数据框来存储平均值
df_avg <- data.frame(row.names = rownames(mf.df))
# 对于每个唯一的前缀，计算平均值或者直接复制单个重复值
for (prefix in unique(prefixes)) {
  # 找到具有相同前缀的列
  pattern <- paste0("^", prefix, "_rep\\d+$")
  cols_to_avg <- grep(pattern, colnames, value = TRUE)
  
  # 如果找到多于一个列，计算平均值；如果只有一个列，直接复制该列
  if (length(cols_to_avg) > 1) {
    # 计算平均值
    df_avg[[prefix]] <- rowMeans(mf.df[, cols_to_avg], na.rm = TRUE)
  } else {
    # 如果只有一个重复值，复制该列并在结果数据框中使用不带后缀的列名
    single_rep_colname <- sub("_rep\\d+$", "", cols_to_avg)
    df_avg[[single_rep_colname]] <- mf.df[, cols_to_avg]
  }
}
colnames(df_avg)[35:46] <- pheno$WGBS[match(colnames(expr_tpm)[35:46],pheno$RNA)]
idx <- match(pheno$WGBS,colnames(df_avg))
mf.m <- df_avg[,idx]
idx <- match(pheno$RNA,colnames(expr_tpm))
expr_tpm <- expr_tpm[,idx]
identical(match(colnames(mf.m),pheno$WGBS),match(colnames(expr_tpm),pheno$RNA))
#save(expr_tpm,mf.m,file = "matched.Rdata")
# range > 1
load("matched.Rdata")
tpm.m <- log2(expr_tpm+1)
idx <- match(rownames(mf.m),rownames(tpm.m))
tpm.m <- tpm.m[idx,]
identical(rownames(tpm.m),rownames(mf.m))
tpm.m <- as.matrix(tpm.m)
mf.m <- as.matrix(mf.m)# 17592
idx_1<- which(apply(tpm.m,1,function(x) {(max(x)-min(x)) <= 1}))
idx_2<- which(apply(mf.m,1,function(x) {(max(x)-min(x)) <= 0.1}))
idx <- c(idx_1,idx_2)
tpm.m <- tpm.m[-idx,]
mf.m <- mf.m[-idx,]
shared_genes <- intersect(rownames(mf.m),rownames(tpm.m))
beta.m <- as.matrix(mf.m[shared_genes,])
tpm.m <- as.matrix(tpm.m[shared_genes,])
save(beta.m,tpm.m,file = "End_matched.Rdata")
# pcc
load("")
library(parallel)
pcc.lv <- unlist(mclapply(1:nrow(beta.m),function(i) {cor(beta.m[i,],tpm.m[i,])}, mc.cores = 10))
p.lv <- unlist(mclapply(1:nrow(beta.m),function(i) {cor.test(beta.m[i,],tpm.m[i,])$p.value}, mc.cores = 10))
p.lv <- p.adjust(p.lv,method = "fdr")
idx <- which(p.lv < 0.05 & pcc.lv < 0)
pcc.sig <- pcc.lv[idx]
p.sig <- p.lv[idx]
imputable_genes <- rownames(beta.m)[idx]
impg.df <- data.frame(imgene=imputable_genes,pcc = pcc.sig, p = p.sig) # 2668
impg.df <- impg.df[order(impg.df$p), ]
#save(impg.df,file = "imputable_genes.Rd")

pdf("imputable_genes.pdf",width = 5,height = 5)
par(cex.main = 0.7)

for (i in 1:nrow(impg.df)) {
  gene_name <- impg.df$imgene[i]
  
  # 绘制基本图形
  plot(tpm.m[gene_name,], beta.m[gene_name,], xlab = "Log2(TPM + 1)", ylab = "beta",
       main = paste0("r = ", round(impg.df$pcc[i], 2), " P = ", format(impg.df$p[i], digits = 2), " (FDR < 0.05)"))
  
  # 标注 E10.5 和 E11.5 阶段的点
  points(tpm.m[gene_name, idx_e10], beta.m[gene_name, idx_e10], col = 'red', pch = 19)
  
  # 标注 E15.5 和 E16.5 阶段的点
  points(tpm.m[gene_name, idx_e16], beta.m[gene_name, idx_e16], col = 'blue', pch = 19)
  
  # 标注 Adult 阶段的点
  points(tpm.m[gene_name, idx_adu], beta.m[gene_name, idx_adu], col = 'green', pch = 19)
  
  # 添加图例到右下角
  legend("topright", legend = c("E10.5 & E11.5", "E15.5 & E16.5", "Adult"), col = c('red', 'blue', 'green'), pch = 19, cex = 0.7)
}

dev.off()


# positive ----------------------------------------------------------------

idx_pos <- which(p.lv < 0.05 & pcc.lv > 0) # 859
pcc.sig <- pcc.lv[idx_pos]
p.sig <- p.lv[idx_pos]
imputable_genes <- rownames(beta.m)[idx_pos]

impg.df <- data.frame(imgene=imputable_genes,pcc = pcc.sig, p = p.sig)
impg.df <- impg.df[order(impg.df$p), ]

idx_e10 <- which(pheno$stage %in% c("E10.5","E11.5"))
idx_e16 <- which(pheno$stage %in% c("E15.5","E16.5"))
idx_adu <- which(pheno$stage %in% c("Adult"))

pdf("positive_genes.pdf", width = 5, height = 5)
par(cex.main = 0.7)

for (i in 1:nrow(impg.df)) {
  gene_name <- impg.df$imgene[i]
  
  # 绘制基本图形
  plot(tpm.m[gene_name,], beta.m[gene_name,], xlab = "Log2(TPM + 1)", ylab = "beta",
       main = paste0("r = ", round(impg.df$pcc[i], 2), " P = ", format(impg.df$p[i], digits = 2), " (FDR < 0.05)"))
  
  # 标注 E10.5 和 E11.5 阶段的点
  points(tpm.m[gene_name, idx_e10], beta.m[gene_name, idx_e10], col = 'red', pch = 19)
  
  # 标注 E15.5 和 E16.5 阶段的点
  points(tpm.m[gene_name, idx_e16], beta.m[gene_name, idx_e16], col = 'blue', pch = 19)
  
  # 标注 Adult 阶段的点
  points(tpm.m[gene_name, idx_adu], beta.m[gene_name, idx_adu], col = 'green', pch = 19)
  
  # 添加图例到右下角
  legend("bottomright", legend = c("E10.5 & E11.5", "E15.5 & E16.5", "Adult"), col = c('red', 'blue', 'green'), pch = 19, cex = 0.7)
}

dev.off()














load("~/Renv/lung_ref/ref_lung.Rd")

impg.marker <- impg.df[match(rownames(lung_ref),impg.df$imgene),]
pdf("imputable_markers_lung.pdf",width = 5,height = 5)
par(cex.main = 0.7)
for (i in 1:nrow(impg.marker)){
  plot(tpm.m[rownames(lung_ref)[i],],beta.m[rownames(lung_ref)[i],],xlab = "Log2(TPM + 1)", ylab = "beta",
       main = paste0("r = ",round(impg.marker$pcc[i],2)," P = ",format(impg.marker$p[i],digits = 2)," (FDR < 0.05)"))
  
}
dev.off()
