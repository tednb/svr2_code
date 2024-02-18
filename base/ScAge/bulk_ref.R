source("~/code/Dmathy/code/Damthy.R")
source("~/code/base/ScAge/scAge.R")
load("~/data/infinium/liver_sample/341_liver_samples/RData/GSE180474_he.Rd")
load("~/data/infinium/liver_sample/40_normal_liver_samples/GSE107038_qc.Rd")
load("~/data/infinium/liver_sample/AA_hepatocytes/GSE123995_qc.Rd")
library(parallel)
## formula
LR <- lm(t(he.o@m)~ age,he.o@s)
LR_summary <- summary(LR)
coef.est <- mclapply(LR_summary,function(x) x$coefficients[, 1],mc.cores = 50)
names(coef.est) <- rownames(he.o@m)
save(coef.est,file = "formula_210.Rd")
## Hep aging DMC
load("~/result/scage++/celldmc_Re.Rd")
idx_hep<-match(rownames(dmct)[dmct[,4]!=0],rownames(he.o@m))
hep.m <- he.o@m[idx_hep,]
####################################################### prediction sets 1
load("formula_210.Rd")
cpg_share <- intersect(rownames(qc_107038.o@m),rownames(he.o@m))
formula <- coef.est[cpg_share]
pcc <- unlist(apply(he.o@m[cpg_share, ], 1,function(x) cor(x,he.o@s$age)))
m_obe <- qc_107038.o@m[cpg_share,]
### preform

age.lv<-compute_probabilities_bulk(m_obe=m_obe,
                                   m.e=formula,
                                   m.c=pcc,
                                   # selection_mode = "percentile",
                                   # CpG_parameter=0.8
                                   selection_mode = "numCpGs",
                                   CpG_parameter = 100
)

age.df <- data.frame(qc_107038.o@s$age,age.lv)
colnames(age.df) <- c("Chronological","scDNAm")
library(ggplot2)
pdf("liver_40_bulkscage.pdf",width = 5,height = 4)
ggplot(age.df, aes(x=Chronological, y=scDNAm)) + 
  geom_point() + 
  geom_smooth(method=lm, se=FALSE) +
  labs(x="Chronological Age", y="DNAm Age", title = "GSE107038: normal liver (n = 39)") +
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 100, by = 10),limits = c(0,100)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10),limits = c(0,100)) +
  # add a 斜率为1截距为0的虚线
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  # pcc
  geom_text(aes(label = paste("PCC = ", round(cor(Chronological, scDNAm), 2), sep = "")), x = 5, y = 90, hjust = 0, vjust = 0, size = 5) +
  # p value
  geom_text(aes(label = paste("P = ", format(cor.test(Chronological, scDNAm)$p.value, digits = 1, scientific = TRUE), sep = "")), x = 5, y = 85, hjust = 0, vjust = 0, size = 5) +
  # MedAE
  geom_text(aes(label = paste("MedAE = ", round(median(abs(Chronological - scDNAm)), 2), sep = "")), x = 5, y = 80, hjust = 0, vjust = 0, size = 5)
dev.off()
#############################################################
########################################################### prediction sets 2

load("formula_210.Rd")
cpg_share_hep <- intersect(rownames(qc.o@m),rownames(he.o@m))
formula_hep <- coef.est[cpg_share_hep]
pcc_hep <- unlist(apply(he.o@m[cpg_share_hep, ], 1,function(x) cor(x,he.o@s$age)))
m_obe_hep <- qc.o@m[cpg_share_hep,]
### preform

age.lv_hep<-compute_probabilities_bulk(m_obe=m_obe_hep,
                                   m.e=formula_hep,
                                   m.c=pcc_hep,
                                   # selection_mode = "percentile",
                                   # CpG_parameter=0.8
                                   selection_mode = "numCpGs",
                                   CpG_parameter = 100
)

age.df_raw <- data.frame(qc.o@s$age,age.lv_hep)
colnames(age.df_raw) <- c("Chronological","scDNAm")
library(ggplot2)
pdf("hep_55_bulkscage_100.pdf",width = 5,height = 4)
ggplot(age.df_raw, aes(x=Chronological, y=scDNAm)) + 
  geom_point() + 
  geom_smooth(method=lm, se=FALSE) +
  labs(x="Chronological Age", y="DNAm Age", title = "GSE123995: hepatocyte culture (n = 55)") +
  theme_bw()+
  scale_x_continuous(breaks = seq(-10, 80, by = 10),limits = c(-10,80)) +
  scale_y_continuous(breaks = seq(-10, 80, by = 10),limits = c(-10,80)) +
  # add a 斜率为1截距为0的虚线
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  # pcc
  geom_text(aes(label = paste("PCC = ", round(cor(Chronological, scDNAm), 2), sep = "")), x = -10, y = 80, hjust = 0, vjust = 0, size = 5) +
  # p value
  geom_text(aes(label = paste("P = ", format(cor.test(Chronological, scDNAm)$p.value, digits = 1, scientific = TRUE), sep = "")), x = -10, y = 75, hjust = 0, vjust = 0, size = 5) +
  # MedAE
  geom_text(aes(label = paste("MedAE = ", round(median(abs(Chronological - scDNAm)), 2), sep = "")), x = -10, y = 70, hjust = 0, vjust = 0, size = 5)
dev.off()
############################################################ prediction sets 1 + hep ADC
cpg_share <- intersect(rownames(qc_107038.o@m),rownames(hep.m))
formula <- coef.est[cpg_share]
pcc <- unlist(apply(hep.m[cpg_share, ], 1,function(x) cor(x,he.o@s$age)))
m_obe <- qc_107038.o@m[cpg_share,]

### preform

age.lv<-compute_probabilities_bulk(m_obe=m_obe,
                                   m.e=formula,
                                   m.c=pcc,
                                   # selection_mode = "percentile",
                                   # CpG_parameter=0.8
                                   selection_mode = "numCpGs",
                                   CpG_parameter = 100
)

age.df <- data.frame(qc_107038.o@s$age,age.lv)
colnames(age.df) <- c("Chronological","scDNAm")
library(ggplot2)
pdf("liver_40_bulkscage_enh.pdf",width = 5,height = 4)
ggplot(age.df, aes(x=Chronological, y=scDNAm)) + 
  geom_point() + 
  geom_smooth(method=lm, se=FALSE) +
  labs(x="Chronological Age", y="DNAm Age", title = "GSE107038: normal liver (n = 39)") +
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 100, by = 10),limits = c(0,100)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10),limits = c(0,100)) +
  # add a 斜率为1截距为0的虚线
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  # pcc
  geom_text(aes(label = paste("PCC = ", round(cor(Chronological, scDNAm), 2), sep = "")), x = 5, y = 90, hjust = 0, vjust = 0, size = 5) +
  # p value
  geom_text(aes(label = paste("P = ", format(cor.test(Chronological, scDNAm)$p.value, digits = 1, scientific = TRUE), sep = "")), x = 5, y = 85, hjust = 0, vjust = 0, size = 5) +
  # MedAE
  geom_text(aes(label = paste("MedAE = ", round(median(abs(Chronological - scDNAm)), 2), sep = "")), x = 5, y = 80, hjust = 0, vjust = 0, size = 5)
dev.off()


############################################################ prediction sets 2 + hep ADC
cpg_share_hep <- intersect(rownames(qc.o@m),rownames(hep.m))
formula_hep <- coef.est[cpg_share_hep]
pcc_hep <- unlist(apply(hep.m[cpg_share_hep, ], 1,function(x) cor(x,he.o@s$age)))
m_obe_hep <- qc.o@m[cpg_share_hep,]
### preform

age.lv_hep<-compute_probabilities_bulk(m_obe=m_obe_hep,
                                       m.e=formula_hep,
                                       m.c=pcc_hep,
                                       # selection_mode = "percentile",
                                       # CpG_parameter=0.8
                                       selection_mode = "numCpGs",
                                       CpG_parameter = 100
)

age.df_enh <- data.frame(qc.o@s$age,age.lv_hep)
colnames(age.df_enh) <- c("Chronological","scDNAm")
library(ggplot2)
pdf("hep_55_bulkscage_enh.pdf",width = 5,height = 4)
ggplot(age.df_enh, aes(x=Chronological, y=scDNAm)) + 
  geom_point() + 
  geom_smooth(method=lm, se=FALSE) +
  labs(x="Chronological Age", y="DNAm Age", title = "GSE123995: hepatocyte culture (n = 55)") +
  theme_bw()+
  scale_x_continuous(breaks = seq(-10, 80, by = 10),limits = c(-10,80)) +
  scale_y_continuous(breaks = seq(-10, 80, by = 10),limits = c(-10,80)) +
  # add a 斜率为1截距为0的虚线
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  # pcc
  geom_text(aes(label = paste("PCC = ", round(cor(Chronological, scDNAm), 2), sep = "")), x = -10, y = 80, hjust = 0, vjust = 0, size = 5) +
  # p value
  geom_text(aes(label = paste("P = ", format(cor.test(Chronological, scDNAm)$p.value, digits = 1, scientific = TRUE), sep = "")), x = -10, y = 75, hjust = 0, vjust = 0, size = 5) +
  # MedAE
  geom_text(aes(label = paste("MedAE = ", round(median(abs(Chronological - scDNAm)), 2), sep = "")), x = -10, y = 70, hjust = 0, vjust = 0, size = 5)
dev.off()
##################################################### A try #####################################

load("celldmc_210.Rd")
load("formula_210.Rd")
idx_hep<-match(rownames(dmct)[dmct[,4]!=0],rownames(he.o@m))
hep.m <- he.o@m[idx_hep,]
cpg_share_hep <- intersect(rownames(qc.o@m),rownames(hep.m))
formula_hep <- coef.est[cpg_share_hep]
m_obe_hep <- qc.o@m[cpg_share_hep,]
t_hep <- coefs$Hep$t
t_hep <- t_hep[match(cpg_share_hep,rownames(he.o@m))]
names(t_hep) <- rownames(m_obe_hep)
pcc_hep <- unlist(apply(hep.m[cpg_share_hep, ], 1,function(x) cor(x,he.o@s$age)))
efsz <- coefs$Hep$Estimate
efsz <- efsz[match(cpg_share_hep,rownames(he.o@m))]
names(efsz) <- rownames(m_obe_hep)

age.lv_try<-compute_probabilities_bulk(m_obe=m_obe_hep,
                                       m.e=formula_hep,
                                       m.c=pcc_hep,
                                       # selection_mode = "percentile",
                                       # CpG_parameter=0.8
                                       selection_mode = "numCpGs",
                                       CpG_parameter = 200
)

age.df_try <- data.frame(qc.o@s$age,age.lv_try)
colnames(age.df_try) <- c("Chronological","scDNAm")
library(ggplot2)
pdf("hep_55_bulkscage_enh_200.pdf",width = 5,height = 4)
ggplot(age.df_try, aes(x=Chronological, y=scDNAm)) + 
  geom_point() + 
  geom_smooth(method=lm, se=FALSE) +
  labs(x="Chronological Age", y="DNAm Age", title = "GSE123995: hepatocyte culture (n = 55)") +
  theme_bw()+
  scale_x_continuous(breaks = seq(-10, 80, by = 10),limits = c(-10,80)) +
  scale_y_continuous(breaks = seq(-10, 80, by = 10),limits = c(-10,80)) +
  # add a 斜率为1截距为0的虚线
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  # pcc
  geom_text(aes(label = paste("PCC = ", round(cor(Chronological, scDNAm), 2), sep = "")), x = -10, y = 80, hjust = 0, vjust = 0, size = 5) +
  # p value
  geom_text(aes(label = paste("P = ", format(cor.test(Chronological, scDNAm)$p.value, digits = 1, scientific = TRUE), sep = "")), x = -10, y = 75, hjust = 0, vjust = 0, size = 5) +
  # MedAE
  geom_text(aes(label = paste("MedAE = ", round(median(abs(Chronological - scDNAm)), 2), sep = "")), x = -10, y = 70, hjust = 0, vjust = 0, size = 5)
dev.off()







































