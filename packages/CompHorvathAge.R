#### ComputeHorvathAge.R

library(gdata);
library(foreign);

horv.df <- read.xls("~/code/packages/Horvath353agepred.xls")
horvageP.lv <- list();
horvageP.lv[[1]] <- as.numeric(horv.df[1,2]);
horvageP.lv[[2]] <- as.numeric(as.vector(horv.df[2:nrow(horv.df),2])); 
names(horvageP.lv[[2]]) <- as.vector(horv.df[2:nrow(horv.df),1]);# extract intersect and coefficients

### transform age functions
ageF <- function(age.v,adult.age=20){
  young.idx <- which(age.v <= adult.age);
  adult.idx <- which(age.v > adult.age);
  tage.v <- age.v;
  tage.v[adult.idx] <- (age.v[adult.idx]-adult.age)/(1+adult.age);
  tage.v[young.idx] <- log(age.v[young.idx]+1) - log(adult.age+1);
  return(tage.v);
}

iageF <- function(ptage.v,adult.age=20){

   y.idx <- which(ptage.v <= 0);
   a.idx <- which(ptage.v>0);
   mage.v <- ptage.v;
   mage.v[a.idx] <- ptage.v[a.idx]*(adult.age+1) + adult.age;
   mage.v[y.idx] <- exp(ptage.v[y.idx] + log(adult.age+1)) -1 ;
   return(mage.v);
 }

PredTage <- function(beta.m){

   match(names(horvageP.lv[[2]]),rownames(beta.m)) -> map.idx;
   rep.idx <- which(is.na(map.idx)==FALSE);
   tmpB.m <- beta.m[map.idx[rep.idx],];
   tage.v <- vector();
   for(s in 1:ncol(tmpB.m)){
     tage.v[s] <- horvageP.lv[[1]] + sum(tmpB.m[,s]*horvageP.lv[[2]][rep.idx]);
   }
   return(tage.v);
}

predTage.v <- PredTage(bmiqGY2.m[,selS.idx]);
predMage.v <- iageF(predTage.v);


pdf("DNAmAgeBCFD.pdf",width=7,height=2.5);
par(mar=c(4,4,2,1));
par(mfrow=c(1,2));
boxplot(predMage.v ~ status.v,pch=23,names=c("N","NADJ","BC"),col=c("green","darkgreen","red"),main="",ylab="DNAmAge(Horvath)",cex=0.5);
nspg.v <- summary(factor(status.v));
mtext(at=1:3,paste("n=",nspg.v,sep=""),side=1,line=1.75,cex=0.75);

ageaccel.v <- predMage.v - ages.v;
boxplot(ageaccel.v ~ status.v,pch=23,names=c("N","NADJ","BC"),col=c("green","darkgreen","red"),main="",ylab="AgeAccel(Horvath)",cex=0.5);
nspg.v <- summary(factor(status.v));
mtext(at=1:3,paste("n=",nspg.v,sep=""),side=1,line=1.75,cex=0.75);
pv01  <- wilcox.test(ageaccel.v[which(status.v==0)],ageaccel.v[which(status.v==1)],alt="less")$p.value;
pv12  <- wilcox.test(ageaccel.v[which(status.v==1)],ageaccel.v[which(status.v==2)],alt="less")$p.value;
text(x=1.5,y=-30,paste("P=",WriteOutPval(pv01,round.max=500),sep=""),font=2);
text(x=2.5,y=40,paste("P=",WriteOutPval(pv12,round.max=500),sep=""),font=2);
dev.off();

break;

library(Hmisc);
c.idx <- which(status.v>0);
rc.o <- rcorr.cens(predMage.v[c.idx],status.v[c.idx]-1);
auc <- rc.o[[1]];
dup <- rc.o[[2]]+1.96*rc.o[[3]]
ddn <- rc.o[[2]]-1.96*rc.o[[3]]
aucUP <- min(0.5*(dup+1),1);
aucDN <- 0.5*(ddn+1);

source("ComputeAUC.R");
auc.o <- ComputeAUC(predMage.v[c.idx],status.v[c.idx]-1);

plot(1-auc.o$sp,auc.o$se,type="l",col="red",xlab="1-Sp",ylab="Se",main="ROC(LCIS vs LCIS->LC)",lwd=2,xlim=c(0,1),ylim=c(0,1));
abline(a=0,b=1,col="green",lwd=2);

text(x=0.6,y=0.2,labels=paste("AUC=",round(auc,2)," (",round(aucDN,2),"-",round(aucUP,2),")",sep=""),cex=1,font=2);




