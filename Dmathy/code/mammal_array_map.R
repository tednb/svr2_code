library(data.table)
source("~/code/Dmathy/code/Damthy.R")
load("~/data/infinium/mammal_array/GSE199979_liver/raw_GSE199979.Rd")
load("~/data/infinium/mammal_array/HBP/raw_GSE224361.Rd")
GPL <- fread("~/data/infinium/GPL_files/GPL28271/GPL28271-57075.txt")

mapping<-function(m,GPL){
  library(parallel)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  idx<- match(rownames(m),GPL$ID)
  mm10_odin <- GPL$Mouse.GRCm38.100_CGstart[idx]
  mod_rela <- GPL$Mouse.GRCm38.100_main_Categories[idx]
  gene_rela <- GPL$Mouse.GRCm38.100_SYMBOL[idx]
  mm10_info <- data.frame(mm10_odin,mod_rela,gene_rela)
  rownames(mm10_info) <- rownames(m)
  m_all <- m
  mm10_info <- na.omit(mm10_info) # coodinate
  m_all <- m_all[rownames(mm10_info),] # preparation
  mm10 <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
  mgene <- mm10$gene_id
  start <- mm10@ranges@start
  end <- end(mm10@ranges)
  strand <- as.character(strand(mm10))
  mm10<-data.frame(start,end,strand)
  rownames(mm10) <- mgene
  mf.m <- matrix(data = NA,ncol = ncol(m_all),nrow = nrow(mm10))
  colnames(mf.m) <- colnames(m_all)
  rownames(mf.m) <- mgene
  
  mm10_odin <- mm10_info$mm10_odin
  
  
  judge <- function(i,cpg){
    if(mm10[i,3] == "+"){
      idx <- which(cpg<=mm10[i,1] & cpg >= mm10[i,1]-200)
    }else{
      idx <- which(cpg>=mm10[i,2] & cpg <= mm10[i,2]+200)
    }
    if(length(idx) > 1){
      mf.m[i,] <- colMeans(m_all[idx,])
    }else if(length(idx) == 1){
      mf.m[i,] <- m_all[idx,]
    }else{
      mf.m[i,] <-NA
    }
    
  }
  mf.m <- do.call(rbind,mclapply(1:nrow(mm10),function(i) {judge(i,mm10_odin)},mc.cores = 100))
  colnames(mf.m) <- colnames(m_all)
  rownames(mf.m) <- mgene
  
  idx <- which(is.na(rowSums(mf.m)))
  mf.m <- mf.m[-idx,]
  
  return(mf.m)
}

mf.m<-mapping(raw.o@raw.m,GPL)



############################################ over #############################
save(mf.m,file = "map_gene_GSE199979.Rd")
############################################ build specific mouse DNAm reference
load("~/data/infinium/mammal_array/GSE199979_liver/map_gene_GSE199979.Rd")
# limma
library(limma)
library(minfi)
load("~/data/WGBS/Loyfer_et_al/ref.Rd")
cellType <- lv$CT
cellType <- factor(cellType)
# use the above to create a design matrix
hypo <- function(fit){
  marker <- list()
  DMPs <- topTable(fit, num = Inf, coef = 1, sort.by = "p")
  DMPs_hypo <- DMPs[which(DMPs$t < 0 & DMPs$adj.P.Val <= 0.05),]
  markers <- rownames(DMPs_hypo)
  print(length(markers))
  return(markers)
}
hyper <- function(fit){
  marker <- list()
  DMPs <- topTable(fit, num = Inf, coef = 1, sort.by = "p")
  DMPs_hyper <- DMPs[which(DMPs$t > 0 & DMPs$adj.P.Val <= 0.05),]
  markers <- rownames(DMPs_hyper)
  print(length(markers))
  return(markers)
}
hypo_c <- list()
# compare one cell type with another three cell types together

med <- as.character(cellType)
other<-setdiff(levels(cellType), "Hepatocyte")
med[med %in% other] <- "other"
med <- factor(med)
design <- model.matrix(~0+med)
colnames(design) <- levels(med)
# fit the linear model
fit <- lmFit(ref.m, design)
# create a contrast matrix for specific comparisons
contrast_formula <- paste("Hepatocyte ", "- other", sep="")
contMatrix <- makeContrasts(contrast_formula, levels = design)
# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
# look at the numbers of DM CpGs at FDR < 0.05
print(summary(decideTests(fit2)))
# filter hypoM
hypo_c[["hepatocyte"]] <- hypo(fit2)
hypo_c[["other"]] <- hyper(fit2)
# match array
library(homologene)
h_m <- homologene(rownames(mf.m),inTax = 10090, outTax = 9606)
m_gen <- h_m$`9606_ID`
hep_marker <- intersect(hypo_c[["hepatocyte"]],m_gen)
rest_marker <- intersect(hypo_c[["other"]],m_gen)
# effect size
by_diff_max <- function(mark,ref.m,CT,cell,THRE){
  m <- ref.m[mark,]
  idx_cell <- which(CT == cell)
  idx_other <- setdiff(1:length(CT),idx_cell)
  g_diff.c <- apply(m,1,function(x) {min(x[idx_other])-max(x[idx_cell])}) # rank
  idx_diff <- which(g_diff.c >= THRE)
  g_diff.c <- g_diff.c[idx_diff]
  mark <- mark[idx_diff]
  names(g_diff.c) <- 1:length(g_diff.c)
  result<-sort(g_diff.c, decreasing = TRUE)
  idx <- as.numeric(names(result))
  print(length(idx))
  return(mark[idx])
}
hep_DMC<-by_diff_max(hep_marker,ref.m,med,"Hepatocyte",THRE = 0.05)
rest_DMC<-by_diff_max(rest_marker,ref.m,med,"other",THRE = 0.05)
hep_DMC <- hep_DMC[1:5] #232
rest_DMC <- rest_DMC[1:5] #232
pdf("hep.pdf",width = 9.5,height = 7)
par(mfrow=c(2,2))
sapply(hep_DMC [c(1,2,length(hep_DMC )-1,length(hep_DMC))], function(cpg){
  plotCpg(ref.m, cpg=cpg, pheno=med, ylab = "methylation fraction")
})
dev.off()
pdf("rest.pdf",width = 9.5,height = 7)
par(mfrow=c(2,2))
sapply(rest_DMC [c(1,2,length(rest_DMC)-1,length(rest_DMC))], function(cpg){
  plotCpg(ref.m, cpg=cpg, pheno=med, ylab = "methylation fraction")
})
dev.off()
###### build
mars <- as.character(unique(c(hep_DMC,rest_DMC)))
mars <- intersect(rownames(ref.m),mars)
ref <- ref.m[mars,]
ref.mx <- matrix(data = NA,nrow = length(mars),ncol = 2)
rownames(ref.mx) <- mars
colnames(ref.mx) <- levels(med)
for ( i in levels(med)){
  idx <- which(med == i)
  if (length(idx) > 1){
    med.c<-apply(ref[,idx],1,median)}
  else {
    med.c<-ref[,idx]
  }
  ref.mx[,i] <- med.c
}
weight <- rep(NA,times = nrow(ref.mx))
ref_liver <- cbind(ref.mx,weight)

h_gene <- rownames(ref_liver)
h_m <- homologene(h_gene,inTax = 9606, outTax = 10090)
idx <- match(h_gene,h_m$`9606_ID`)
h_ID <- h_m$`9606_ID`[idx]
m_name <- h_m$`10090_ID`[idx]


idx <- which(!is.na(m_name))
mice_ref <- ref_liver[idx,]
rownames(mice_ref) <- na.omit(m_name)



##############################################
save(mice_ref,file = "mref_form40k.Rd")
############################################### deconvolution
### RPC
load("~/Renv/ref_for_mmal/mref_form40k.Rd")
RPC <- function(mice_ref,mf.m){
  library(MASS)
  common.v <- intersect(rownames(mice_ref),rownames(mf.m))
  map.idx <- match(common.v,rownames(mf.m))
  rep.idx <- match(common.v,rownames(mice_ref))
  data.m <- mf.m[map.idx,]
  ref.m <- mice_ref[rep.idx,-ncol(mice_ref)]
  est.m <- matrix(data = NA,nrow=ncol(data.m),ncol=ncol(ref.m)) 
  colnames(est.m) <- colnames(ref.m)
  rownames(est.m) <- colnames(data.m)
  sum_squared_residuals_liver <- vector("numeric", ncol(data.m))
  for (s in 1:ncol(data.m)) {
    rlm.o <- rlm(data.m[, s] ~ ref.m, maxit = 500)
    coef.v <- summary(rlm.o)$coef[2:(ncol(ref.m)+1),1];
    coef.v[which(coef.v<0)] <- 0;
    total <- sum(coef.v);
    coef.v <- coef.v/total;
    est.m[s,] <- coef.v;
  }
  # pdf("liver_fraction.pdf",width = 6, height = 4)
  # boxplot(est.m, main = "Estimation of liver fractions", xlab = "Cell types", ylab = "Fraction")
  # dev.off()
  est.m
}
idx_liver <- raw.o@raw.s$tissue == "Liver"
mf_liver.m <- mf.m[,idx_liver]
est<-RPC(mice_ref,mf_liver.m) # 4.69

hepf.lv <- list()
for ( i in levels(factor(raw.o@raw.s$tissue))) {
  idx <- which(raw.o@raw.s$tissue == i)
  m <- mf.m[,idx]
  est.m <- RPC(mice_ref,m)
  hepf.lv[[i]] <- est.m[,1]
}
pdf("hep_f_rpc.pdf",width = 10,height = 5)
par(mar=c(8,4,4,4))
boxplot(hepf.lv,ylab = "Fraction",las = 3,main = "Hepatocyte")
abline(h=0.55)
dev.off()

### QP

QP <- function(mice_ref,mf.m,choice = "<="){
  common.v <- intersect(rownames(mice_ref),rownames(mf.m))
  map.idx <- match(common.v,rownames(mf.m))
  rep.idx <- match(common.v,rownames(mice_ref))
  data.m <- mf.m[map.idx,]
  ref.m <- mice_ref[rep.idx,-ncol(mice_ref)]
  library(quadprog)
  nCT <- ncol(ref.m)
  D <- 2 * apply(ref.m, 2, function(x) colSums(x * ref.m)) #
  
  if (choice == "=") {
    coe.v <- c(1, 1)
  } else coe.v <- c(-1, 0)
  
  A.m <- matrix(0, nrow = nCT, ncol = nCT)
  diag(A.m) <- rep(1, nCT)
  A.m <- cbind(rep(coe.v[1], nCT), A.m)
  b0.v <- c(coe.v[1], rep(0, nCT))
  
  est.m <- matrix(nrow=ncol(data.m),ncol=ncol(ref.m))
  colnames(est.m) <- colnames(ref.m)
  rownames(est.m) <- colnames(data.m)
  
  residual_sum_of_squares_liver <- 0
  
  for (s in seq_len(ncol(data.m))) {
    data<-matrix(data.m[,s], nrow = 1)
    d.v <- as.vector(2 * data %*% ref.m)
    qp.o <- solve.QP(D, d.v, A.m, b0.v, meq = coe.v[2])
    
    est.m[s, ] <- qp.o$sol
    #print(qp.o$sol)
    residuals <- ref.m %*% qp.o$sol - data.m[,s]
    
    
    residual_sum_of_squares_liver[s] <- sum(residuals^2)/nrow(data.m)
    #print(sum(residuals^2))
  }
  #return(residual_sum_of_squares_liver)
  return(est.m)
}

idx_liver <- raw.o@raw.s$tissue == "Liver"
mf_liver.m <- mf.m[,idx_liver]
GOF_liver<-QP(mice_ref,mf_liver.m)
idx_other <- setdiff(1:nrow(raw.o@raw.s),idx_liver)
mf_other.m <- mf.m[,idx_other]
GOF_other<-QP(mice_ref,mf_other.m)

GOF<-list(GOF_liver,GOF_other)
names(GOF) <- c("liver","other")
pdf("goodness.pdf",width = 5,height = 5)
# draw boxplot
boxplot(GOF, ylab="SSR", xlab="Tissue", main="goodness of fit", col=c("lightblue","lightgreen"))
# 使用Wilcoxon rank-sum test比较两列数据的差异，将p值添加到图中
p <- wilcox.test(GOF_liver,GOF_other)$p.value
text(1.5, 0.08, paste("p=", round(p, 3)), cex=1)
dev.off()
########################################## function transform

### CIBERSORT
CBS <- function(beta.m, ref.m,nu.v = c(0.25, 0.5, 0.75)) {
  library(e1071)
  
  map.idx <- match(rownames(ref.m), rownames(beta.m))
  rep.idx <- which(is.na(map.idx) == FALSE)
  
  data2.m <- beta.m[map.idx[rep.idx], , drop=FALSE]
  ref2.m <- ref.m[rep.idx, -ncol(ref.m)]
  
  est.lm <- list()
  nui <- 1
  for (nu in nu.v) {
    est.m <- matrix(nrow = ncol(data2.m), ncol = ncol(ref2.m))
    colnames(est.m) <- colnames(ref2.m)
    rownames(est.m) <- colnames(data2.m)
    for (s in seq_len(ncol(data2.m))) {
      svm.o <- svm(x = ref2.m, y = data2.m[, s], scale = TRUE, type = "nu-regression", 
                   kernel = "linear", nu = nu)
      coef.v <- t(svm.o$coefs) %*% svm.o$SV
      coef.v[which(coef.v < 0)] <- 0
      total <- sum(coef.v)
      coef.v <- coef.v/total
      est.m[s, ] <- coef.v
    }
    est.lm[[nui]] <- est.m
    nui <- nui + 1
  }
  
  #### select best nu
  rmse.m <- matrix(NA, nrow = ncol(beta.m), ncol = length(nu.v))
  for (nui in seq_along(nu.v)) {
    reconst.m <- ref2.m %*% t(est.lm[[nui]])
    s <- seq_len(ncol(beta.m))
    rmse.m[s, nui] <- sqrt(colMeans((data2.m[, s, drop = FALSE] - reconst.m[, s,  drop = FALSE])^2))
    message(nui)
  }
  colnames(rmse.m) <- nu.v
  nu.idx <- apply(rmse.m, 1, which.min)
  estF.m <- est.m
  for (s in seq_len(nrow(estF.m))) {
    estF.m[s, ] <- est.lm[[nu.idx[s]]][s, ]
  }
  return(list(estF = estF.m, nu = nu.v[nu.idx], ref = ref2.m, dataREF = data2.m))
}
















