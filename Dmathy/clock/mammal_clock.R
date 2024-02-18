load("~/data/infinium/mammal_array/GSE199979_liver/raw_GSE199979.Rd")
load("~/data/MouseDNAmAtlas/ENCODE/tpm+location.Rdata")
source("~/code/Dmathy/code/hts_map_to_tss_mm10.R")
load("~/Renv/liver_EpiSCORE/Data/mrefMouseliver.Rd")
#load("~/data/infinium/mammal_array/GSE199979_liver/map_gene_GSE199979.Rd")
source("~/code/Dmathy/code/Damthy.R")
source("~/code/Dmathy/code/ELN_clock.R")
load("~/Renv/ref_for_mmal/mref_form40k.Rd")
load("~/data/scage/single_cell/hepatocyte/hep_26.Rd")
ph.lv <- read.csv("~/data/scage/single_cell/hepatocyte/SraRunTable.txt")
GPL <- fread("~/data/infinium/GPL_files/GPL28271/GPL28271-57075.txt")
load("~/data/scage/single_cell/hepatocyte/schep_age.Rd")
library(data.table)
library(EpiDISH)
{
#################################### cell type fraction #################################

est <- RPC(mice_ref,mf.m)
m_raw <- raw.o@raw.m
pheno <- raw.o@raw.s
pheno$sex <- NULL
colnames(m_raw) <- rownames(pheno)

#################################### Aging CpGs in hepatocyte #############################

qc.o <- new("qc.o")
qc.o@m <- m_raw
qc.o@ctf.o[[1]] <- est
qc.o@s <- pheno
# age dtr #
library(ggplot2)
data <- data.frame(
  value = qc.o@s$age
)
pdf("age_dtr.pdf",width = 6,height = 5)
ggplot(data, aes(x=value)) +
  geom_histogram(binwidth=10, color="black", fill="blue", aes(y=..count..))+
  labs(x="Age (days)", title="N = 97") +
   theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))

dev.off()

## SVD

qc.o@svd.o <- lm_svd(qc.o)
p_h(qc.o,pic = c(15,1))
save(qc.o,file = "mmal_liver_qc.Rd")
}
#1# CellDMC 
load("~/data/infinium/mammal_array/GSE199979_liver/mmal_liver_qc.Rd")

covs <- model.matrix(~qc.o@s$strain + qc.o@s$diet)
celldmc.o <- CellDMC(qc.o@m, qc.o@s$age ,qc.o@ctf.o[[1]],
                     cov.mod = covs,
                     adjPMethod = "fdr",
                     adjPThresh = 0.05,
                     sort = FALSE,
                     mc.cores = 50)
dmct<- celldmc.o$dmct
idx_cdc <- which(dmct[,2] != 0) # 2948




################################# ELN clock #############################
# trans coordinate
idx<- match(rownames(qc.o@m),GPL$ID)
mm10_odin <- paste(paste0("chr",GPL$Mouse.GRCm38.100_seqnames[idx]), GPL$Mouse.GRCm38.100_CGstart[idx], sep = "-")
#hep_dmc <- na.omit(mm10_odin[idx_cdc]) # 2874

idx <- which(!is.na(mm10_odin))
m_all <- qc.o@m[idx,]
rownames(m_all) <- mm10_odin[idx]
# 
qc_mm10.o <- new("qc.o")
qc_mm10.o@m <- m_all
qc_mm10.o@s <- qc.o@s
save(qc_mm10.o,file = "~/data/infinium/mammal_array/GSE199979_liver/mmal_liver_codinversion.Rd")
# map
load("~/data/infinium/mammal_array/GSE199979_liver/mmal_liver_codinversion.Rd")
mf_tss500 <- map_to_tss(qc_mm10.o@m,gene_info,region = 500)




#### Train ####
load("~/data/infinium/mammal_array/GSE199979_liver/mmal_liver_codinversion.Rd")

seqs <- seq(0.5,1.5,0.001)
result <- bench(qc_mm10.o,seqs,nP = 10)
clock_raw <- result[[1]]
rawclock_cpg <- clock_raw$All$beta@Dimnames[[1]][clock_raw$All$beta@i]
save(clock_raw,file = "mmal_rawclock.Rd")

# hep specific clock
covs <- model.matrix(~qc_mm10.o@s$strain + qc_mm10.o@s$diet)
result_hep <- clock(qc_mm10.o,covs=covs,seqs,nP = 10)
clock_hep <- result_hep[[1]]
hepclock_cpg <- clock_hep$Hep$beta@Dimnames[[1]][clock_hep$Hep$beta@i]
save(clock_hep,file = "mmal_hepclock.Rd")

#There are only 2 overlapped clock cpgs between hep and raw 
# test
############################################## HBP ########################
load("~/data/infinium/mammal_array/HBP/raw_GSE224361.Rd")
idx_liver <- which(raw.o@raw.s$tissue == "Liver")
liver.m <- raw.o@raw.m[,idx_liver]
liver.lv <- raw.o@raw.s[idx_liver,]

idx<- match(rownames(liver.m),GPL$ID)
mm10_odin <- GPL$Mouse.GRCm38.100_CGstart[idx]


### CTF with aging 
load("~/Renv/ref_for_mmal/mref_form40k.Rd")
mf.m<-mapping(liver.m,GPL)
hep_f<-RPC(mice_ref,mf.m)[,1]
group <- liver.lv$group
group <- gsub("\\d", "", group)
f_g<-data.frame(hep_f,factor(group))
colnames(f_g) <- c("fraction","group")
idx <- f_g$group %in% c("YoungIso","OldIso","OldHet","YoungIsoRec","OldIsoRec","OldHetRec")
f_g_hbp <- f_g[idx,]
pdf("hep_f_hbpgroup.pdf",width = 5,height = 4)
ggplot(f_g_hbp, aes(x=group, y=fraction))+
  scale_x_discrete(limits=c("YoungIso","OldIso","OldHet","YoungIsoRec","OldIsoRec","OldHetRec"))+
  geom_boxplot()+
  theme_bw() +geom_jitter(aes(color=group),alpha=0.5)+
  labs(title="Hepatocyte", x="", y="Fraction")+
  theme(legend.position="none")+
  geom_segment(aes(x=2,xend=3,y=0.59,yend=0.59))+
  annotate("text",x=2.5,y=0.592,label=paste("p = ",format(t.test(f_g$fraction[f_g$group=="OldIso"],f_g$fraction[f_g$group=="OldHet"],alternative="two.sided",var.equal=FALSE)$p.value,digits=1,scientific=TRUE),sep=""))+
  geom_segment(aes(x=5,xend=6,y=0.59,yend=0.59))+
  annotate("text",x=5.5,y=0.592,label=paste("p = ",format(t.test(f_g$fraction[f_g$group=="OldIsoRec"],f_g$fraction[f_g$group=="OldHetRec"],alternative="two.sided",var.equal=FALSE)$p.value,digits=1,scientific=TRUE),sep=""))
dev.off()

idx <- f_g$group %in% c("ST_YoungIso","ST_OldIso","ST_OldHet")
f_g_stbp <- f_g[idx,]
pdf("hep_f_stgroup.pdf",width = 5,height = 4)
ggplot(f_g_stbp, aes(x=group, y=fraction))+
  scale_x_discrete(limits=c("ST_YoungIso","ST_OldIso","ST_OldHet"))+
  geom_boxplot()+
  theme_bw() +geom_jitter(aes(color=group),alpha=0.5)+
  labs(title="Hepatocyte", x="", y="Fraction")+
  theme(legend.position="none")+
  geom_segment(aes(x=2,xend=3,y=0.59,yend=0.59))+
  annotate("text",x=2.5,y=0.592,label=paste("p = ",format(t.test(f_g$fraction[f_g$group=="OldIso"],f_g$fraction[f_g$group=="OldHet"],alternative="two.sided",var.equal=FALSE)$p.value,digits=1,scientific=TRUE),sep=""))
dev.off()  




##  prepare
idx <- which(!is.na(mm10_odin))
liver.m <- liver.m[idx,]
rownames(liver.m) <- mm10_odin[idx]
qc.liver <- new("qc.o")
qc.liver@m <- liver.m
qc.liver@s <- liver.lv
### age prediction
gse <- "GSE224361"
age_all <- validate(qc.liver,clock_raw,gse)
age_raw.df <- data.frame(qc.liver@s$age*365,age_all)
colnames(age_raw.df) <- c("chronological","DNAm")
group <- qc.liver@s$group
group <- gsub("\\d", "", group)
age_raw.df$group <- factor(group)
age_hepc <- validate(qc.liver,clock_hep,gse)
age_hep.df <- data.frame(qc.liver@s$age*365,age_hepc)
colnames(age_hep.df) <- c("chronological","DNAm")
group <- qc.liver@s$group
group <- gsub("\\d", "", group)
age_hep.df$group <- factor(group)

################ test clock
{#1 raw clock
  library(ggplot2)
idx <- grep("Iso",age_raw.df$group)
age_ISO.df <- age_raw.df[idx,]
pdf("raw_clock_IBP.pdf",width = 5,height = 4)
ggplot(age_ISO.df,aes(x=chronological,y=DNAm,color=group)) + 
  geom_point() + 
  theme_bw()+
  geom_abline(intercept=0,slope=1,linetype="dashed")+
  scale_x_continuous(limits=c(0,900),breaks=seq(0,900,100)) + 
  scale_y_continuous(limits=c(0,900),breaks=seq(0,900,100)) +
  xlab("Chronological age (days)") + ylab("DNAm age (days)") + ggtitle("IBP liver (n = 32)")+
  annotate("text",x=0,y=900,label=paste("PCC = ",round(cor(age_ISO.df$chronological,age_ISO.df$DNAm),2),sep=""),hjust=0)+
  annotate("text",x=0,y=850,label=paste("P = ",format(cor.test(age_ISO.df$chronological,age_ISO.df$DNAm)$p.value,digits=1,scientific=TRUE),sep=""),hjust=0)+
  annotate("text",x=0,y=800,label=paste("RMSE = ",format(sqrt(mean((age_ISO.df$chronological-age_ISO.df$DNAm)^2)),digits=1),sep=""),hjust=0)
dev.off()

#2 hep clock
idx <- grep("Iso",age_hep.df$group)
age_ISO.df <- age_hep.df[idx,]
pdf("hep_clock_IBP.pdf",width = 5,height = 4)
ggplot(age_ISO.df,aes(x=chronological,y=DNAm,color=group)) + 
  geom_point() + 
  theme_bw()+
  geom_abline(intercept=0,slope=1,linetype="dashed")+
  scale_x_continuous(limits=c(0,900),breaks=seq(0,900,100)) + 
  scale_y_continuous(limits=c(0,900),breaks=seq(0,900,100)) +
  xlab("Chronological age (days)") + ylab("DNAm age (days)") + ggtitle("IBP liver (n = 32)")+
  annotate("text",x=0,y=900,label=paste("PCC = ",round(cor(age_ISO.df$chronological,age_ISO.df$DNAm),2),sep=""),hjust=0)+
  annotate("text",x=0,y=850,label=paste("P = ",format(cor.test(age_ISO.df$chronological,age_ISO.df$DNAm)$p.value,digits=1,scientific=TRUE),sep=""),hjust=0)+
  annotate("text",x=0,y=800,label=paste("RMSE = ",format(sqrt(mean((age_ISO.df$chronological-age_ISO.df$DNAm)^2)),digits=1),sep=""),hjust=0)
dev.off()
}
################
{### compare different group
idx <- age_raw.df$group %in% c("YoungIso","OldIso","OldHet","YoungIsoRec","OldIsoRec","OldHetRec")
age_exp.df <- age_raw.df[idx,]
age_exp.df$chronological <- NULL

pdf("raw_clock_HBP_boxplot.pdf",width = 5,height = 4)
ggplot(age_exp.df,aes(x=group,y=DNAm)) + geom_boxplot() + theme_bw() + 
  scale_x_discrete(limits=c("YoungIso","OldIso","OldHet","YoungIsoRec","OldIsoRec","OldHetRec"))+
  xlab("") + ylab("DNAm age (days)") + ggtitle("The role of HBP, All clock (n = 29)")+
  geom_jitter(aes(color=group),width=0.2,size=2)+
  guides(color="none")+
  # 在第一列第二列之间添加一条线并标注p值
  geom_segment(aes(x=1,xend=2,y=820,yend=820))+
  annotate("text",x=1.5,y=850,label=paste("p = ",format(t.test(age_exp.df$DNAm[age_exp.df$group=="YoungIso"],age_exp.df$DNAm[age_exp.df$group=="OldIso"],alternative="less",var.equal=FALSE)$p.value,digits=1,scientific=TRUE),sep=""))+
  # 在第三列第四列之间添加一条线并标注p值
  geom_segment(aes(x=2,xend=3,y=850,yend=850))+
  annotate("text",x=2.5,y=870,label=paste("p = ",format(t.test(age_exp.df$DNAm[age_exp.df$group=="OldIso"],age_exp.df$DNAm[age_exp.df$group=="OldHet"],alternative="greater",var.equal=FALSE)$p.value,digits=1,scientific=TRUE),sep=""))+
  geom_segment(aes(x=4,xend=5,y=870,yend=870))+
  annotate("text",x=4.5,y=900,label=paste("p = ",format(t.test(age_exp.df$DNAm[age_exp.df$group=="YoungIsoRec"],age_exp.df$DNAm[age_exp.df$group=="OldIsoRec"],alternative="less",var.equal=FALSE)$p.value,digits=1,scientific=TRUE),sep=""))+
  # 在第五列第六列之间添加一条线并标注p值
  geom_segment(aes(x=5,xend=6,y=920,yend=920))+
  annotate("text",x=5.5,y=950,label=paste("p = ",format(t.test(age_exp.df$DNAm[age_exp.df$group=="OldIsoRec"],age_exp.df$DNAm[age_exp.df$group=="OldHetRec"],alternative="greater",var.equal=FALSE)$p.value,digits=1,scientific=TRUE),sep=""))
dev.off()

idx <- age_hep.df$group %in% c("YoungIso","OldIso","OldHet","YoungIsoRec","OldIsoRec","OldHetRec")
age_exp.df <- age_hep.df[idx,]
age_exp.df$chronological <- NULL

pdf("hep_clock_HBP_boxplot.pdf",width = 5,height = 4)
ggplot(age_exp.df,aes(x=group,y=DNAm)) + geom_boxplot() + theme_bw() + 
  scale_x_discrete(limits=c("YoungIso","OldIso","OldHet","YoungIsoRec","OldIsoRec","OldHetRec"))+
  xlab("") + ylab("DNAm age (days)") + ggtitle("The role of HBP, Hep clock (n = 29)")+
  geom_jitter(aes(color=group),width=0.2,size=2)+
  guides(color="none")+
  # 在第一列第二列之间添加一条线并标注p值
  geom_segment(aes(x=1,xend=2,y=820,yend=820))+
  annotate("text",x=1.5,y=850,label=paste("p = ",format(t.test(age_exp.df$DNAm[age_exp.df$group=="YoungIso"],age_exp.df$DNAm[age_exp.df$group=="OldIso"],alternative="less",var.equal=FALSE)$p.value,digits=1,scientific=TRUE),sep=""))+
  # 在第三列第四列之间添加一条线并标注p值
  geom_segment(aes(x=2,xend=3,y=850,yend=850))+
  annotate("text",x=2.5,y=870,label=paste("p = ",format(t.test(age_exp.df$DNAm[age_exp.df$group=="OldIso"],age_exp.df$DNAm[age_exp.df$group=="OldHet"],alternative="greater",var.equal=FALSE)$p.value,digits=1,scientific=TRUE),sep=""))+
  geom_segment(aes(x=4,xend=5,y=870,yend=870))+
  annotate("text",x=4.5,y=900,label=paste("p = ",format(t.test(age_exp.df$DNAm[age_exp.df$group=="YoungIsoRec"],age_exp.df$DNAm[age_exp.df$group=="OldIsoRec"],alternative="less",var.equal=FALSE)$p.value,digits=1,scientific=TRUE),sep=""))+
  # 在第五列第六列之间添加一条线并标注p值
  geom_segment(aes(x=5,xend=6,y=920,yend=920))+
  annotate("text",x=5.5,y=950,label=paste("p = ",format(t.test(age_exp.df$DNAm[age_exp.df$group=="OldIsoRec"],age_exp.df$DNAm[age_exp.df$group=="OldHetRec"],alternative="greater",var.equal=FALSE)$p.value,digits=1,scientific=TRUE),sep=""))
dev.off()
 }
### Short term
idx <- grep("ST",age_hep.df$group)
ST_hep <- age_hep.df[idx,]
ST_hep$chronological <- NULL
pdf("hep_clock_STHBP_boxplot.pdf",width = 5,height = 4)
ggplot(ST_hep,aes(x=group,y=DNAm)) + geom_boxplot() + theme_bw() + 
  scale_x_discrete(limits=c("ST_YoungIso","ST_OldIso","ST_OldHet"))+
  xlab("") + ylab("DNAm age (days)") + ggtitle("The role of short time HBP, Hep clock (n = 19)")+
  geom_jitter(aes(color=group),width=0.2,size=2)+
  guides(color="none")+
  
  geom_segment(aes(x=1,xend=2,y=820,yend=820))+
  annotate("text",x=1.5,y=850,label=paste("p = ",format(t.test(ST_hep$DNAm[ST_hep$group=="ST_YoungIso"],ST_hep$DNAm[ST_hep$group=="ST_OldIso"],alternative="less",var.equal=FALSE)$p.value,digits=1,scientific=TRUE),sep=""))+
  
  geom_segment(aes(x=2,xend=3,y=850,yend=850))+
  annotate("text",x=2.5,y=870,label=paste("p = ",format(t.test(ST_hep$DNAm[ST_hep$group=="ST_OldIso"],ST_hep$DNAm[ST_hep$group=="ST_OldHet"],alternative="greater",var.equal=FALSE)$p.value,digits=1,scientific=TRUE),sep=""))
dev.off()  

idx <- grep("ST",age_raw.df$group)
ST_raw <- age_raw.df[idx,]
ST_raw$chronological <- NULL
pdf("raw_clock_STHBP_boxplot.pdf",width = 5,height = 4)
ggplot(ST_raw,aes(x=group,y=DNAm)) + geom_boxplot() + theme_bw() + 
  scale_x_discrete(limits=c("ST_YoungIso","ST_OldIso","ST_OldHet"))+
  xlab("") + ylab("DNAm age (days)") + ggtitle("The role of short time HBP, All clock (n = 19)")+
  geom_jitter(aes(color=group),width=0.2,size=2)+
  guides(color="none")+
  
  geom_segment(aes(x=1,xend=2,y=820,yend=820))+
  annotate("text",x=1.5,y=850,label=paste("p = ",format(t.test(ST_raw$DNAm[ST_raw$group=="ST_YoungIso"],ST_raw$DNAm[ST_raw$group=="ST_OldIso"],alternative="less",var.equal=FALSE)$p.value,digits=1,scientific=TRUE),sep=""))+
  
  geom_segment(aes(x=2,xend=3,y=850,yend=850))+
  annotate("text",x=2.5,y=870,label=paste("p = ",format(t.test(ST_raw$DNAm[ST_raw$group=="ST_OldIso"],ST_raw$DNAm[ST_raw$group=="ST_OldHet"],alternative="greater",var.equal=FALSE)$p.value,digits=1,scientific=TRUE),sep=""))
dev.off() 

################################################ partial reprogramming #######
load("~/Renv/clock_mmal/mmal_hepclock.Rd")
load("clock_mmal/mmal_rawclock.Rd")
load("~/data/infinium/mammal_array/partial_reprogramming/gse190665/raw_GSE190665.Rd")
idx_liver <- which(raw.o@raw.s$tissue == "Liver")
liver.m <- raw.o@raw.m[,idx_liver]
liver.lv <- raw.o@raw.s[idx_liver,]

idx<- match(rownames(liver.m),GPL$ID)
mm10_odin <- GPL$Mouse.GRCm38.100_CGstart[idx]

idx <- which(!is.na(mm10_odin))
liver.m <- liver.m[idx,]
rownames(liver.m) <- mm10_odin[idx]


qc.repo_liver <- new("qc.o")
qc.repo_liver@m <- liver.m
qc.repo_liver@s <- liver.lv

gse <- "GSE190665"
age_all <- validate(qc.repo_liver,clock_raw,gse)
age_hepc <- validate(qc.repo_liver,clock_hep,gse)
### comparison
age_raw.df <- data.frame(round(qc.repo_liver@s$age*365,0),age_all)
colnames(age_raw.df) <- c("chronological","DNAm")
group <- liver.lv$group
age_raw.df$group <- factor(gsub("^\\d.", "", group))

age_raw.df$chronological <- NULL
agegroup <- rep(c(639, 820, 780, 91, 690,780), times = c(10, 3, 2, 4, 6,3))
age_raw.df$agegroup <- factor(agegroup)




age_hep.df <- data.frame(round(qc.repo_liver@s$age*365,0),age_hepc)
colnames(age_hep.df ) <- c("chronological","DNAm")
group <- liver.lv$group
age_hep.df $group <- factor(gsub("^\\d.", "", group))

age_hep.df $chronological <- NULL
agegroup <- rep(c(639, 820, 780, 91, 690,780), times = c(10, 3, 2, 4, 6,3))
age_hep.df $agegroup <- factor(agegroup)


## all samples
pdf("raw_clock_rpo_boxplot.pdf",width = 5,height = 4)
ggplot(age_raw.df,aes(x=group,y=DNAm)) + geom_boxplot() + theme_bw()+
  # 展示点，颜色按年龄区分，不要渐变图例，图例标题为chronological age
  geom_jitter(aes(color=agegroup),alpha=0.5) +
  guides(color=guide_legend(title="Chronological age"))+
  xlab("") + ylab("DNAm age (days)") + ggtitle("Partial reprogramming liver, All clock (n = 20)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
pdf("hep_clock_rpo_boxplot.pdf",width = 5,height = 4)
ggplot(age_hep.df,aes(x=group,y=DNAm)) + geom_boxplot() + theme_bw()+
  # 展示点，颜色按年龄区分，不要渐变图例，图例标题为chronological age
  geom_jitter(aes(color=agegroup),alpha=0.5) +
  guides(color=guide_legend(title="Chronological age"))+
  xlab("") + ylab("DNAm age (days)") + ggtitle("Partial reprogramming liver, hep clock (n = 20)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## part
idx <- which(age_raw.df$agegroup == 639)
age_raw_639 <- age_raw.df[idx,]
ggplot(age_raw_639,aes(x=group,y=DNAm)) + geom_boxplot() + theme_bw()+
  geom_jitter() +
  guides(color=guide_legend(title="Chronological age"))+
  xlab("") + ylab("DNAm age (days)") + ggtitle("Partial reprogramming liver, raw clock (n = 10)")

age_hep_639 <- age_hep.df[idx,]

ggplot(age_hep_639,aes(x=group,y=DNAm)) + geom_boxplot() + theme_bw()+
  geom_jitter() +
  guides(color=guide_legend(title="Chronological age"))+
  xlab("") + ylab("DNAm age (days)") + ggtitle("Partial reprogramming liver, raw clock (n = 10)")








pdf("raw_clock_repo.pdf",width = 6,height = 5)
ggplot(age_raw.df,aes(x=chronological,y=DNAm,color=group)) + 
  geom_point() + 
  theme_bw()+
  geom_abline(intercept=0,slope=1,linetype="dashed")+
  scale_x_continuous(limits=c(0,900),breaks=seq(0,900,100)) + 
  scale_y_continuous(limits=c(0,900),breaks=seq(0,900,100)) +
  xlab("Chronological age (days)") + ylab("DNAm age (days)") + ggtitle("Partial reprogramming liver (n = 28)")
dev.off()

age_hep.df <- data.frame(qc.repo_liver@s$age*365,age_hepc)
colnames(age_hep.df) <- c("chronological","DNAm")
group <- qc.repo_liver@s$condition
age_hep.df$group <- group
pdf("hep_clock_repo.pdf",width = 6,height = 5)
ggplot(age_hep.df,aes(x=chronological,y=DNAm,color=group)) + 
  geom_point() + 
  theme_bw()+
  geom_abline(intercept=0,slope=1,linetype="dashed")+
  scale_x_continuous(limits=c(0,900),breaks=seq(0,900,100)) + 
  scale_y_continuous(limits=c(0,900),breaks=seq(0,900,100)) +
  xlab("Chronological age (days)") + ylab("DNAm age (days)") + ggtitle("Partial reprogramming liver (n = 28)")
dev.off()












