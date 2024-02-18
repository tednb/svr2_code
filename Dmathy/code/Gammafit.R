load("~/data/MouseDNAmAtlas/ENCODE/tpm.Rd")
library(dplyr)
library(tibble)
library(limma)
tpm.m <- as.data.frame(tpm.m) %>% rownames_to_column(., var = "Ensembl_ID")
tpm.m$ENCSR611PTP <- NULL
gtf_vM21 <- rtracklayer::import('~/chain/gencode.vM21.basic.annotation.gtf.gz') %>% as.data.frame()
gtf_mrna_vM21 <- dplyr::select(gtf_vM21,c("gene_id","gene_type", "gene_name"))%>%
  subset(., gene_type == "protein_coding") %>% unique()
expr_tpm <- tpm.m %>%
  inner_join(gtf_mrna_vM21, by = c("Ensembl_ID" = "gene_id")) %>%
  select(gene_name, starts_with("ENCS") )
expr_tpm = as.matrix(avereps(expr_tpm[,-1],ID = expr_tpm$gene_name) )
expr_tpm <- log2(expr_tpm + 1)
### GammaFit.R
library(mixtools);
### assume that your normalized log2(TPM+1) expression matrix for all protein coding genes and 46 samples are in expr_tpm
tmp.v <- as.vector(expr_tpm)
epsilon <- 1e-6
gfit.o <-  gammamixEM(tmp.v+epsilon,lambda = NULL, alpha = NULL, beta = NULL, k = 2,verb=TRUE, maxit = 50000);
save(gfit.o,file = "EM_re.Rd")
# print(gfit.o[2:4]);
lambda.v <- gfit.o$lambda;
parsFit.m <- gfit.o$gamma.pars;

pdf("gamma_fit.pdf",width = 6,height = 4)
plot(density(tmp.v), xlim=c(min(tmp.v)-2, max(tmp.v)+2), main = "Mixture of Gamma distribution",xlab = "Log(TPM+1)")

colors <- c("green", "magenta","black") # 定义颜色向量以便复用和引用
legends <- c("non-expression", "expression","True distribution") # 定义每条线对应的图例文本

for(k in 1:2){
  d.o <- density(rgamma(length(tmp.v), shape=parsFit.m[1,k], scale=parsFit.m[2,k]))
  points(d.o$x, lambda.v[k]*d.o$y, type="l", col=colors[k])
}

# 添加图例
legend("topright", # 图例位置
       legend = legends, # 图例文本
       col = colors, # 图例颜色
       lty = 1, # 线型（这里是实线）
       cex = 0.8,
       title = "Fitting distribution") # 图例缩放比例

dev.off()
# 
# pXgE.m <- pgamma(expr_tpm,shape=parsFit.m[1,2],scale=parsFit.m[2,2],lower.tail=T)
# pXgNE.m <- pgamma(expr_tpm,shape=parsFit.m[1,1],scale=parsFit.m[2,1],lower.tail=F)
# print(cbind(expr_tpm[1:10,1],pXgE.m[1:10,1],pXgNE.m[1:10,1]));
# 
# pEgX.m <- pXgE.m*lambda.v[2];
# pNEgX.m <- pXgNE.m*lambda.v[1];
# 
# ### estimate probability of expression given the observed value
# ### find Bayesian posterior
# 
# sum.m <- pEgX.m + pNEgX.m;
# pEgX.m <- pEgX.m/sum.m;
# pNEgX.m <- pNEgX.m/sum.m;
# print(cbind(expr_tpm[1:10,1],pEgX.m[1:10,1],pNEgX.m[1:10,1]));


### with pEgX.m you can now proceed by imputing DNAm values for imputable marker genes in the cell-types where it is not expressed by using the median DNAm value for the samples in your ENCODE dataset where that gene is not expressed using the threshold pEgX.m[g,] < 0.2
load("Data/dbENCODE.Rd")
post.v <- as.vector(pEgXeode.m)
plot(density(post.v), xlim=c(min(post.v)-0.2, max(post.v)+0.2), main = "Posterior probability of gene expression",xlab = "Posterior")
abline(v=0.2,col = "red")
