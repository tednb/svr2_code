load("~/data/infinium/MM285/Rdata/normal_solid/map_DB.Rd")
load("~/data/infinium/MM285/Rdata/normal_solid/hey.saminfo.Rd")
load("~/Renv/DMC2/mouse_ref.Rd")


idx_liver <- which(anno$tissue == "Liver")
idx_other <- setdiff(1:nrow(anno),idx_liver)
map_gene_liver.m <- map_gene.m[,idx_liver]
map_gene_other.m <- map_gene.m[,idx_other]

common.v <- intersect(rownames(mice_ref),rownames(map_gene_liver.m))
map.idx <- match(common.v,rownames(map_gene_liver.m))
rep.idx <- match(common.v,rownames(mice_ref))
data_liver.m <- map_gene_liver.m[map.idx,] # 61 * 98
data_other.m <- map_gene_other.m[map.idx,]
data_all.m <- map_gene.m[map.idx,]
ref.m <- mice_ref[rep.idx,-ncol(mice_ref)] # 61 * 4

# RPC
library(MASS)
RPC <- function(ref.m,data.m){
  est.m <- matrix(data = NA,nrow=ncol(data.m),ncol=ncol(ref.m)) 
  colnames(est.m) <- colnames(ref.m)
  rownames(est.m) <- colnames(data.m)
sum_squared_residuals_liver <- vector("numeric", ncol(data.m))
for (s in 1:ncol(data.m)) {
  rlm.o <- rlm(data.m[, s] ~ ref.m, maxit = 500)
  coef.v <- summary(rlm.o)$coef[2:(ncol(ref.m)+1),1];
  # coef.v[which(coef.v<0)] <- 0;
  # total <- sum(coef.v);
  # coef.v <- coef.v/total;
  est.m[s,] <- coef.v;

  residuals <- rlm.o$residuals
  sum_squared_residuals_liver[s] <- sum(residuals^2)  # 平均平方残差
}
# pdf("MMA_other.pdf",width = 6, height = 4)
boxplot(est.m, main = "Estimation of liver fractions", xlab = "Cell types", ylab = "Fraction")
# dev.off()
#mean(sum_squared_residuals_liver)
est.m[,1]
}

RPC(ref.m,data_liver.m) # 4.69
RPC(ref.m,data_other.m) # 4.60

hepf.lv <- list()
for ( i in levels(anno$tissue)) {
  idx <- which(anno$tissue == i)
  m <- map_gene.m[,idx]
  data.m <- m[map.idx,]
  hepf.lv[[i]] <- RPC(ref.m,data.m)
}
pdf("hep_f_rpc_o.pdf",width = 10,height = 5)
par(mar=c(8,4,4,4))
boxplot(hepf.lv,ylab = "Fraction",las = 3,main = "Hepatocyte")
dev.off()


#QP/CP

QP <- function(ref.m,data.m,choice = "<="){
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

  
  residual_sum_of_squares_liver[s] <- sum(residuals^2)
  #print(sum(residuals^2))
}
# pdf("MMA_other.pdf",width = 6, height = 4)
boxplot(est.m, main = "Estimation of liver fractions", xlab = "Cell types", ylab = "Fraction")
# dev.off()
return(residual_sum_of_squares_liver)
#return(est.m[,1])

}

sr_liver<-QP(ref.m,data_liver.m) # 4.714
QP(ref.m,data_liver.m,"=") #4.819

sr_other<-QP(ref.m,data_other.m) # 4.763
QP(ref.m,data_other.m,"=") #4.818

# validation
ssr.lv <- list()
for ( i in levels(anno$tissue)) {
  idx <- which(anno$tissue == i)
  m <- map_gene.m[,idx]
  data.m <- m[map.idx,]
  ssr.lv[[i]] <- QP(ref.m,data.m)
}
pdf("fitness_<.pdf",width = 15,height = 5)
par(mar=c(8,4,4,4))
boxplot(ssr.lv,ylab = "SSR",las = 3)
dev.off()

hepf.lv <- list()
for ( i in levels(anno$tissue)) {
  idx <- which(anno$tissue == i)
  m <- map_gene.m[,idx]
  data.m <- m[map.idx,]
  hepf.lv[[i]] <- QP(ref.m,data.m)
}
pdf("hep_f_qp.pdf",width = 10,height = 5)
par(mar=c(8,4,4,4))
boxplot(hepf.lv,ylab = "Fraction",las = 3,main = "Hepatocyte")
dev.off()
# for



