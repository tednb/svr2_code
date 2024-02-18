load("~/data/infinium/MM285/Rdata/Gencode.Rd")
load("~/data/infinium/MM285/Rdata/MMA_liver_qc.Rd")

### 90 rrbs
# load("~/data/RRBS/mouse/GSE60012/GSE60012.Rd")
### 30 rrbs
# load("~/data/scage/REF/merge.Rd")
# load("~/data/scage/REF/sample_info.Rd")
idx <- which(pheno$tissue == "Liver")
m_all <- m_all[,idx]
pheno <- pheno[idx,]
###
source("~/code/Dmathy/code/Damthy.R")

### single cell
load("~/data/scage/single_cell/hepatocyte/hep_26.Rd")
ph.lv <- read.csv("~/data/scage/single_cell/hepatocyte/SraRunTable.txt")
age_hep <- as.numeric(sapply(strsplit(ph.lv$Age, " "), function(x) x[1]))
idx_liver <- match(names(hep.lv),ph.lv$Run)
age_hep <- age_hep[idx_liver]
age_hep[which(is.na(age_hep))] <- 4
###

source("~/code/base/ScAge/scAge.R")
library(parallel)


########################################## formula
LR <- lm(t(m_all)~ age, data = pheno)
LR_summary <- summary(LR)
coef.est <- mclapply(LR_summary,function(x) x$coefficients[, 1],mc.cores = 50)
names(coef.est) <- rownames(m_all)
### mm285
age <- qc.o@s$age
age <- round(age/30,digits = 3)
LR <- lm(t(qc.o@m)~ age)
LR_summary <- summary(LR)
coef.est <- mclapply(LR_summary,function(x) x$coefficients[, 1],mc.cores = 50)
idx <- match(rownames(qc.o@m),gencode$probeID)
names(coef.est) <- gencode$CpG_beg[idx]
m_tmp <- qc.o@m
rownames(m_tmp) <- gencode$CpG_beg[idx]
hep.pcc <- lapply(hep.lv, function(cell) {
  share <- intersect(names(cell),names(coef.est))
  c <- unlist(apply(m_tmp[share,],1,function(x) cor(x,age)))
  return(c)
})
### mm285

# hep.c <- lapply(hep.lv, function(cell) {
#   share <- intersect(names(cell),names(coef.est))
#   c <- coef.est[share]qc.o@mqc.o@m
#   efs <- sapply(c,function(x) x[2])
#   names(efs) <- share
#   return(efs)
# })

###### rrbs
hep.e <- lapply(hep.lv, function(cell) {
  share <- intersect(names(cell),names(coef.est))
  c <- coef.est[share]
  return(c)
})

# PCC
hep.pcc <- lapply(hep.lv, function(cell) {
  share <- intersect(names(cell),names(coef.est))
  c <- unlist(apply(m_all[share,],1,function(x) cor(x,pheno$age)))
  return(c)
})

###### rrbs
ages <- list()

# for(i in 1:length(hep.lv)){
#   age_pre<-compute_probabilities(cell=hep.lv[[i]],
#                                  cell.e=hep.e[[i]],
#                                  cell.c=hep.pcc[[i]],
#                                  # selection_mode = "percentile",
#                                  # CpG_parameter=0.9,
#                                  selection_mode = "numCpGs",
#                                  CpG_parameter = 100)
#   print(age_pre)
#   ages <- append(ages,age_pre)
# }



for(i in 1:length(hep.lv)){
  age_pre<-compute_probabilities_parallel(cell=hep.lv[[i]],
                                 cell.e=hep.e[[i]],
                                 cell.c=hep.pcc[[i]],
                                 # selection_mode = "percentile",
                                 # CpG_parameter=0.9,
                                 selection_mode = "numCpGs",
                                 CpG_parameter = 100)
  print(age_pre)
  ages <- append(ages,age_pre)
}

# 21 cell
age.df<-data.frame(age_hep,unlist(ages),factor(age_hep))
age.df <- na.omit(age.df)
colnames(age.df) <- c("Chronological","scDNAm","group")
library(ggplot2)
pdf("mm285_90ref.pdf",width = 6,height = 5)
ggplot(age.df, aes(x = Chronological, y = scDNAm)) +
  geom_boxplot(aes(fill = group), alpha = 0.7) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Chronological Age", y = "scDNAm Age") +
  geom_jitter(width = 0.5, height = 0) +
  scale_x_continuous(breaks = seq(0, 40, by = 10),limits = c(0,40)) +
  scale_y_continuous(breaks = seq(0, 40, by = 10),limits = c(0,40)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1)+
  annotate("text", x = 5, y = 35, label = paste("R = ", round(cor(age.df$Chronological,age.df$scDNAm), digits = 2), sep = ""), size = 5)+
  annotate("text", x = 5, y = 33, label = paste("p = ", format(cor.test(age.df$Chronological,age.df$scDNAm)$p.value, digits = 1, scientific = TRUE), sep = ""), size = 5)+
  annotate("text", x = 5, y = 30, label = paste("MedAE = ", round(median(abs(age.df$Chronological-age.df$scDNAm)), digits = 2), sep = ""), size = 5)
dev.off()
