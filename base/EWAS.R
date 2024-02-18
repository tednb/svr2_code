setwd("Renv")
load("BFrac.rda")
load("/mnt/local-disk/data/guoxiaolong/data/BBC/dataBBC.Rd")
load("/mnt/local-disk/data/guoxiaolong/data/BBC/PhenoTypesBBC.Rd")
load("/mnt/local-disk/data/guoxiaolong/data/BBC/annoEPICv1B2.Rd")

library(qvalue) #FDR correction
library(ggplot2)
library(broom) #tidy() and glance()
library(pbapply) #display progress bar
library(ggpubr) #QQ plot
### EWAS of DNAm vs smoking adjusting for age and sex ###

## only three independent variables
sample_anno <- data.frame(PhenoTypesBBC.lv) #sample information
bmiqBBC.m <- t(bmiqBBC.m) # reverse
variables <- subset(sample_anno, select = c(Sex, Age, SMK)) #extract interested variables
pri_LR <- lm(bmiqBBC.m~Sex + Age + SMK, data = variables) #linear regression for each CpG
sum_LR <- summary(pri_LR)
# first filter : select significant results of LR
ftest_p <- pbsapply(sum_LR, function(x) glance(x)[[5]]) #extract p-value of f-test
adj_fp <- qvalue(ftest_p, fdr.level = 0.05) #FDR method corrects p-value > 0.05 of f-test
sig_sumLR <- sum_LR[adj_fp$significant] #filter CpG
# second filter : select CpGs related to smoking
SMK_p <- pbsapply(sig_sumLR, function(x) tidy(x)[[5]][4])
adj_smkp <- qvalue(SMK_p, fdr.level = 0.05)
SMK_sumLR <- sig_sumLR[adj_smkp$significant] # 17811

## put cell fractions in the model
cell_var <- subset(sample_anno, select = c(Neu, Mon, CD4T, CD8T, B, NK, Sex, Age, SMK)) #the values of Eos.bas are zero
sec_LR <- lm(bmiqBBC.m~Neu + Mon + CD4T + CD8T + B + NK + Sex + Age + SMK, data = cell_var)
sum_secLR <- summary(sec_LR)
# first filter : select significant results of LR
s_ftest_p <- pbsapply(sum_secLR, function(x) glance(x)[[5]])
s_adj_fp <- qvalue(s_ftest_p, fdr.level = 0.05)
sig_ssLR <- sum_secLR[s_adj_fp$significant]
# second filter : select CpGs related to smoking
s_SMK_p <- pbsapply(sig_ssLR, function(x) tidy(x)[[5]][10])
s_adj_smkp <- qvalue(s_SMK_p, fdr.level = 0.05)
SMK_ssLR <- sig_ssLR[s_adj_smkp$significant] # 129

## Quantile-quantile plot

# for result with cell fractions
qvals <- s_adj_smkp$qvalues[s_adj_smkp$significant]
qdf <- data.frame(obs=-log10(sort(qvals,decreasing=FALSE)),
                  exp=-log10(ppoints(length(qvals)))+1)
ggplot(data=qdf,aes(exp,obs))+
  geom_point(alpha=0.7,color="#7F7F7FFF")+
  geom_abline(color="#D62728FF")+
  xlab("Expected -log10")+
  ylab("Observed -log10")+
  scale_x_continuous(limits = c(0,4))+
  scale_y_continuous(limits = c(0,14))+
  theme(
    axis.title = element_text(size=12,face="bold"),
    axis.text = element_text(face="bold",size=8,color = "black"),
    axis.ticks= element_line(size=0.8,colour = "black"),
    panel.grid =element_blank(),
    panel.border = element_rect(fill=NA,size = 0.8),
    panel.background = element_blank())

# for result without cell fractions
fqvals <- adj_smkp$qvalues[adj_smkp$significant]
fqdf <- data.frame(obs=-log10(sort(fqvals,decreasing=FALSE)),
                  exp=-log10(ppoints(length(fqvals)))+1)
ggplot(data=fqdf,aes(exp,obs))+
  geom_point(alpha=0.7,color="#7F7F7FFF")+
  geom_abline(color="#D62728FF")+
  xlab("Expected -log10")+
  ylab("Observed -log10")+
  scale_x_continuous(limits = c(0,5))+
  scale_y_continuous(limits = c(0,5))+
  theme(
    axis.title = element_text(size=12,face="bold"),
    axis.text = element_text(face="bold",size=8,color = "black"),
    axis.ticks= element_line(size=0.8,colour = "black"),
    panel.grid =element_blank(),
    panel.border = element_rect(fill=NA,size = 0.8),
    panel.background = element_blank())
