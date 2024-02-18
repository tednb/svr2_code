library(TCA)
library(glmnet)
source("~/code/Dmathy/code/ELN_clock.R")

tca_fuc <- function(m,ctf,C1,C2=NULL){
  tca.o <- tca(X=m, W=ctf, C1 = C1
               , C2 = C2, refit_W = FALSE
               , refit_W.features = NULL, refit_W.sparsity = 100
               , refit_W.sd_threshold = 0.02, parallel = TRUE
               , num_cores = 5, max_iters = 10, log_file = "TCA.log", debug = FALSE)
  #p_age.df <- tca.o$gammas_hat_pvals[,paste(rownames(ctf),".age",sep="")]
  postmuTCA.lm <- tensor(X = m, tca.mdl = tca.o, parallel = TRUE, num_cores = 30)
  return(list(tca.o,postmuTCA.lm))
}

setGeneric("glm_tca_cv",function(obj,...){
  standardGeneric("glm_tca_cv")
})

setMethod("glm_tca_cv","pp",function(obj,seqs,C1,nP=10,alpha=1){

m <- obj@m
age <- obj@s$age
ctf <- obj@ctf.o[[1]]
# base matrix
print("base TCA...")
# re_all<-tca_fuc(m,ctf,C1)
# save(re_all,file = "All_sample_hep.Rd")
load("~/Renv/TCA/All_sample_hep.Rd")
p_age_all <- re_all[[1]]$gammas_hat_pvals[,paste(colnames(ctf),".age",sep="")]
hep_ap<-p.adjust(p_age_all[,3],method = "BH")
hep_aidx <- which(hep_ap<=0.05)
# chol_ap<-p.adjust(p_age_all[,1],method = "BH")
# chol_aidx <- which(chol_ap<=0.05)
idx <- seq(1,ncol(m),1)
n.idx <- idx
#bags <- list();
# for(p in 1:nP){
#   bags[[p]] <- sample(idx,ncol(m)/nP,replace=FALSE)
#   idx <- setdiff(idx,bags[[p]])
# }

predH.m <- matrix(NA,nrow=length(seqs),ncol=ncol(m))
#predC.m <- matrix(NA,nrow=length(seqs),ncol=ncol(m))

hepDMCT.lm <- list()
cholDMCT.lm <- list()
test.li <- list()
train.li <- list()

for(p in 1:nP){
  load(paste0("~/Renv/TCA/",p,".Rd"))
  #test.idx <- bags[[p]]
  #train.idx <- setdiff(n.idx,test.idx)
  train.idx <- match(colnames(re[[2]][[1]]),colnames(m))
  test.idx <- setdiff(n.idx,train.idx)
  test.li[[p]] <- test.idx
  train.li[[p]] <- train.idx
  p_age.df <- re[[1]]$gammas_hat_pvals[,paste(colnames(ctf),".age",sep="")]
  hep_p<-p.adjust(p_age.df[,3],method = "BH")
  hepDMCT.lm[[p]] <- which(hep_p<=0.05)
  #cholDMCT.lm[[p]] <- rownames(dmcts[[p]])[which(dmcts[[p]][,1] != 0)]
  #cat("Got: ",length(hepDMCT.lm[[p]])," hep age-DMCs","\n","Got: ",length(cholDMCT.lm[[p]])," chol age-DMCs","\n")
}

#m_hep <- re_all[[2]][[3]]

for(p in 1:nP){
  print(paste0(p," TCA"))
  load(paste0("~/Renv/TCA/",p,".Rd"))
  hep.idx <- hepDMCT.lm[[p]]
  #chol.idx <- match(cholDMCT.lm[[p]],rownames(m_tr))
  
  #model_chol <- glmnet(t(m_tr[chol.idx,train.li[[p]]]),y=age_tr[train.li[[p]]],alpha=1,standardize=TRUE,lambda=seqs)
  model_hep <- glmnet(t(m[hep.idx,train.li[[p]]]),y=age[train.li[[p]]],alpha=1,standardize=TRUE,lambda=seqs)
  #model_hep <- glmnet(t(re[[2]][[3]]),y=age[train.li[[p]]],alpha=0.5,standardize=TRUE,lambda=seqs)
  match(test.li[[p]],n.idx) -> maptoN.idx
  predH.m[,maptoN.idx] <- t(predict(model_hep,newx=t(m[hep.idx,test.li[[p]]])))
  #predH.m[,maptoN.idx] <- t(predict(model_hep,newx=t(re_all[[2]][[3]][,test.li[[p]]])))
  # match(test.li[[p]],n.idx) -> maptoN.idx
  # predC.m[,maptoN.idx] <- t(predict(model_chol,newx=t(m_tr[chol.idx,test.li[[p]]])))
  
}




# for(p in 1:nP){
#   print(paste0(p," TCA"))
#   test.idx <- bags[[p]]
#   train.idx <- setdiff(n.idx,test.idx)
#   # re<-tca_fuc(m[,train.idx],ctf[train.idx,],C1[train.idx,])
#   # save(re,file = paste0(p,".Rd"))
#   #CTS age-DMC
#   load(paste0("~/Renv/TCA/",p,".Rd"))
#   p_age.df <- re[[1]]$gammas_hat_pvals[,paste(colnames(ctf),".age",sep="")]
#   hep_p<-p.adjust(p_age.df[,3],method = "BH")
#   hep_idx <- which(hep_p<=0.05)
#   model_hep <- glmnet(t(m[hep_idx,train.idx]),y=age[train.idx],alpha=alpha,standardize=TRUE,lambda=seqs)
#   
#   
#   # chol_p<-p.adjust(p_age.df[,1],method = "BH")
#   # chol_idx <- which(chol_p<=0.05)
#   # model_chol <- glmnet(t(m[chol_idx,train.idx]),y=age[train.idx],alpha=alpha,standardize=TRUE,lambda=seqs)
#   # 
#   predH.m[,train.idx] <- t(predict(model_hep,newx=t(m[hep_idx,test.idx]))) # prevent overfiting
#   #predC.m[,train.idx] <- t(predict(model_chol,newx=t(m[chol_idx,test.idx]))) # prevent overfiting
#   
# }

rmseH.v <- vector()
for(i in 1:nrow(predH.m)){
  rmseH.v[i] <- sqrt(mean((predH.m[i,] - age)^2))
}
la_hep <- rev(seqs)[which.min(rmseH.v)]

# rmseC.v <- vector()
# for(i in 1:nrow(predC.m)){
#   rmseC.v[i] <- sqrt(mean((predC.m[i,] - age)^2))
# }
# la_chol <- rev(seqs)[which.min(rmseC.v)]

pdf("OptimizationCurve-Hep.pdf",width=5,height=3);
par(mar=c(4,4,2,1))
plot(rev(seqs),rmseH.v,ylab="RMSE",xlab="Lambda",type="l",main = "Hepatocyte")
# add a vertical line at the optimal lambda
abline(v=la_hep,lty=2,col="red")
dev.off()

age_hep <- predH.m[which.min(rmseH.v),]
r <- cor(age,age_hep,method = "pearson")
# calculate the p-value
p <- cor.test(age,age_hep,method = "pearson")$p.value
# Calculate the Median Absolute Error between predicted age and real age
mae <- median(abs(age - age_hep))
# draw scatter plot of predicted age vs. chronological age
pdf("pre_vs_chro_hep.pdf",width=5,height=5)
plot(age,age_hep,pch = 16,col = "blue",ylab = "Predicted age (years)",xlab = "Chronological age (years)",main = "Cholangiocyte",xlim = c(10,95),ylim = c(10,95))
abline(0,1,col="red")
text(20,90,paste0("r = ",format(r, digits = 3),"\n","p = ",format(p, digits = 1),"\n","MedAE: ",format(mae, digits = 2),"years"))
dev.off()

clocks_hep <- glmnet(t(m[hep_aidx,]),age,family="gaussian",alpha = alpha,lambda = la_hep)

#clocks_try <- glmnet(t(re_all[[2]][[3]][hep_aidx,]),age,family="gaussian",alpha = alpha,lambda = la_hep)
return(list(clocks_hep))


})

