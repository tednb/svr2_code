load("~/data/infinium/purified_blood/EGAD00010000850/blue_print_qc.Rd")
source("~/code/Dmathy/code/Damthy.R")
library(broom)
library(isva)

Donor <- names(which(table(qc.o@s$Donor) == 3))
idx <- as.numeric(sapply(Donor,grep,qc.o@s$Donor))
qc.o@s <- qc.o@s[idx,]
qc.o@m <- qc.o@m[,idx]

# batch correction
library(sva)

m_rbat <- do.call(cbind, lapply(levels(factor(qc.o@s$Cell)),function(i){
  idx <- which(qc.o@s$Cell == i)
  ComBat(qc.o@m[,idx], batch = qc.o@s$Batch[idx])
}))
### SVD ###
s <- qc.o@s[match(colnames(m_rbat),rownames(qc.o@s)),]
sam_var <- s[,c("Donor","Cell","Sentrix_ID","Sentrix_Position","Sex","Age")]
sam_var$Donor <- factor(sam_var$Donor)
sam_var$Cell <- factor(sam_var$Cell)
sam_var$Sentrix_ID <- factor(sam_var$Sentrix_ID)
sam_var$Sentrix_Position <- factor(sam_var$Sentrix_Position)
sam_var$Sex <- ifelse(sam_var$Sex == "M",1,0)


m <- m_rbat - rowMeans(m_rbat)
n <- EstDimRMT(m)$dim
svd.m <- svd(m)
v_top <- svd.m$v[,1:n]
pval.m <- matrix(ncol = n,nrow = ncol(sam_var))

for (i in 1:ncol(sam_var)){
  print(i)
  if (class(sam_var[[i]]) != "factor"){
    pval <- apply(v_top, 2,function(x) {tidy(summary(lm(x ~ sam_var[[i]])))[[5]][2]})
  }else {
    pval <- apply(v_top, 2,function(x) {glance(summary(lm(x ~ sam_var[[i]])))[[5]][1]})
  }
  pval.m[i,] <- pval
}
colnames(pval.m) <- paste0("PC-",1:n)
rownames(pval.m) <- colnames(sam_var)
qc.o@svd.o <- list(pval.m,svd.m)
p_h(qc.o)
#### PCA #####
qc.o@m <- m_rbat
qc.o@s <- s
qc.o@s$tissue <- qc.o@s$Cell
qc.o@s$tissue <- factor(qc.o@s$Batch)
qc.o@pca.o <- PCA(qc.o,3)
save(qc.o,file = "blue_print_qc_nobatch.Rd")

##################### set case  ################
load("~/Renv/simulation/blue_print_qc_nobatch.Rd")
source("~/code/Dmathy/code/Damthy.R")
library(broom)
library(qvalue)
idx <- which(qc.o@s$Cell == "CD4_N")
m_cd4 <- qc.o@m[,idx]
m_all <- qc.o@m
age <- qc.o@s$Age
# CD4 age-DMCs
cd4_md <- summary(lm(t(m_cd4)~age[idx]))
p_cd4 <- unlist(mclapply(cd4_md, function(x) glance(x)[[5]],mc.cores = 100))
q_cd4 <- qvalue(p_cd4, fdr.level = 0.05)
m_cd4_sig <- m_cd4[q_cd4$significant,] #26331 age_DMCs
# set.seed(1234)
# random_rows <- sample(1:nrow(m_cd4_sig), size = 1000, replace = FALSE)
# m_cd4_sig <- m_cd4_sig[random_rows,]
#19 cases
age_info<-unique(data.frame(donor = qc.o@s$Donor, age = qc.o@s$Age))
idx <- which(age_info$age <= 70 & age_info$age >= 60)
set.seed(1234)
idx <- sample(idx,19)
case_idv <- age_info$donor[idx]
cd4tro_idx <- which(qc.o@s$Cell == "CD4_N")
cd4tro_donors <- qc.o@s$Donor[cd4tro_idx]
matching_samples <- cd4tro_idx[cd4tro_donors %in% case_idv]
sample_name <- rownames(qc.o@s)[matching_samples]
idx_case <- match(sample_name,colnames(m_cd4_sig))
idx_control <- setdiff(1:ncol(m_cd4_sig),idx_case)
# manipulate
control_md <- lm(t(m_cd4_sig)~qc.o@s[colnames(m_cd4_sig),]$Age)
sum_cmd <- summary(control_md)
eftsize<-unlist(mclapply(sum_cmd, function(x) tidy(x)[[2]][2],mc.cores = 10))
intc <- unlist(mclapply(sum_cmd, function(x) tidy(x)[[2]][1],mc.cores = 10))
res.m <- t(residuals(control_md))
age <- qc.o@s[colnames(m_cd4_sig),]$Age
# build case
m_new <- matrix(NA, nrow = nrow(m_cd4_sig), ncol = ncol(m_cd4_sig))
m_new[,idx_control] <- as.matrix(eftsize) %*% age[idx_control]
m_new[,idx_case] <- as.matrix(2 * eftsize) %*% age[idx_case]
intc.m <- matrix(intc, nrow = length(intc), ncol = ncol(m_cd4_sig), byrow = FALSE)
m_new <- m_new + intc.m + res.m
colnames(m_new) <- colnames(m_cd4_sig)
rownames(m_new) <- rownames(m_cd4_sig)
# combine
idx_row <- match(rownames(m_new), rownames(m_all))
idx_col <- match(colnames(m_new), colnames(m_all))
m_all[idx_row,idx_col] <- m_new
pheno <- qc.o@s
save(m_all,pheno, file = "fake_matrix.Rd")
save(sample_name,case_idv,file = "case_env.Rd")
####################### mixtures #################
setwd("simulation/")
load("case_env.Rd")
load("fake_matrix.Rd")
identical(rownames(pheno),colnames(m_all)) # check
idx_case <- as.numeric(sapply(case_idv,grep,pheno$Donor))
m_case <- m_all[,idx_case]
length(intersect(colnames(m_case),sample_name)) #check
idv_case <- pheno$Donor[idx_case]
ct_case <- pheno$Cell[idx_case]


fracs <- function(n, length, mode = 1, min_val = 0.1, max_val = 0.8) {
  if (mode == 1) {
    # 模式1: 所有数相等
    vec <- rep(1/n, n) # 创建一个所有元素都是 1/n 的向量
    list_of_vecs <- replicate(length, vec, simplify = FALSE) # 创建一个列表，包含指定数量的向量
  } else if (mode == 2) {
    # 模式2: 所有数随机生成
    list_of_vecs <- lapply(seq_len(length), function(...) {
      repeat {
        vec <- runif(n - 1, min = .Machine$double.eps, max = 1 - .Machine$double.eps) # 在(0,1)范围内随机生成n-1个数
        vec <- c(vec, 0, 1) # 加入0和1，构成n+1个数
        vec <- sort(vec) # 排序
        res <- diff(vec) # 计算相邻数字的差值，构成n个数
        if (sum(res) == 1 && all(res != 0)) return(res) # 如果n个数的和为1且都不为0，则返回这个向量
      }
    })
  }else if (mode == 3) {
    # 模式3: 所有数平均，但每个数有随机噪声，且每个数是非负的，且每个数在[min_val, max_val]范围内
    list_of_vecs <- lapply(seq_len(length), function(...) {
      repeat {
        vec <- rep(1/n, n) + rnorm(n, 0, 0.2)
        vec <- pmax(pmin(vec, max_val), min_val) # 将所有大于max_val的数替换为max_val，所有小于min_val的数替换为min_val
        if (sum(vec) > 0) { 
          vec <- vec / sum(vec) # 归一化
          return(vec)
        }
      }
    })
  } else {
    stop("Invalid mode: mode must be 1,2 or 3")
  }
  
  return(list_of_vecs) # 返回结果列表
}
mixf <- function(mat, celltype, people,mode = 1) {
  # 检查输入参数的长度
  if (length(celltype) != ncol(mat) || length(people) != ncol(mat)) {
    stop("Length of celltype and people must be equal to the number of columns in mat")
  }
  fraction.lv <- fracs(length(levels(factor(celltype))),length(unique(people)),mode = mode)
  names(fraction.lv) <- unique(people)
  # 初始化结果矩阵
  mixed_mat <- matrix(0, nrow = nrow(mat), ncol = length(unique(people)))
  colnames(mixed_mat) <- unique(people)
  rownames(mixed_mat) <- rownames(mat)
  
  # 按照个体和细胞类型分组，计算混合样本
  for (i in seq_along(fraction.lv)) {
    fractions <- fraction.lv[[i]]
    person <- names(fraction.lv)[i]
    for (j in seq_along(fractions)) {
      cell <- levels(factor(celltype))[j]
      print(cell)
      # 找到与当前个体和细胞类型对应的列
      cols <- which(people == person & celltype == cell)
      mixed_mat[, person] <- mixed_mat[, person] + mat[, cols] * fractions[j]
    }
  }
  
  return(list(mixed_mat,fraction.lv))
}
# case mixture: same fraction
m_case_1<-mixf(m_case,ct_case,idv_case,mode = 1)[[1]]
any(is.na(m_case_1))
# case mixture: random fraction
m_case_2<-mixf(m_case,ct_case,idv_case,mode = 2)[[1]]
any(is.na(m_case_2))
### 20 CONTROL TESTSET
m_control <- m_all[,-idx_case]
age_info<-unique(data.frame(donor = pheno$Donor, age = pheno$Age))
age_info <- age_info[-match(case_idv,age_info$donor),]
idx <- which(age_info$age <= 70 & age_info$age >= 60) # many and same with case
set.seed(1234)
idx <- sample(idx,20)
control_idv <- age_info$donor[idx]
control_sample <- rownames(pheno)[as.numeric(sapply(control_idv,grep,pheno$Donor))]
m_control <- m_control[,control_sample]
ct_control <- pheno[control_sample,]$Cell
idv_control <- pheno[control_sample,]$Donor
# control mixture: same fraction
m_control_1<-mixf(m_control,ct_control,idv_control,mode = 1)[[1]]
any(is.na(m_control_1))
# control mixture: random fraction
m_control_2<-mixf(m_control,ct_control,idv_control,mode = 2)[[1]]
any(is.na(m_control_2))
save(m_control_1,m_control_2,m_case_1,m_case_2,file = "testset.Rd")
##100 train set
idx <- match(control_sample,colnames(m_all[,-idx_case]))
m_train <- m_all[,-idx_case][,-idx]
ct_train <- pheno[colnames(m_train),]$Cell
idv_train <- pheno[colnames(m_train),]$Donor
save(m_train, ct_train, idv_train, file = "Trainset.Rd")

############################ EpiDish + CellDMC ##########################
library(ggplot2)
source("~/code/Dmathy/code/ELN_clock.R")
load("Trainset.Rd")
re <- mixf(m_train,ct_train,idv_train,mode = 2)
m_train <- re[[1]]
frac.lv <-  re[[2]]

library(EpiDISH)
est <- as.data.frame(epidish(beta.m = m_train, ref.m = centDHSbloodDMC.m, method = "RPC")$estF)
est <- data.frame(CD4 = est$CD4T,Mono = est$Mono,Neu = est$Neutro)
cell <- levels(factor(ct_train))

for (l in 1:length(frac.lv[[1]])){
  pre_frac <- est[,l]
  name <- cell[l]
  true_frac <- sapply(frac.lv,"[[",l)
  # using ggplot2 to plot the esti frac vs true frac
  pdf(paste0(name,"_re.pdf"),width = 4,height = 4)
  p <- ggplot(data.frame(pre_frac,true_frac),aes(x = true_frac,y = pre_frac))+
    geom_point()+
    geom_abline(intercept = 0,slope = 1)+
    labs(x = "True Fraction",y = "Estimated Fraction",title = name)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.position = "none")+
    # 添加PCC
    annotate("text",x = 0.1,y = 0.8,label = paste0("PCC = ",round(cor(pre_frac,true_frac),2)))+
    # 添加RMSE
    annotate("text",x = 0.1,y = 0.7,label = paste0("RMSE = ",round(sqrt(mean((pre_frac-true_frac)^2)),2)))+
    # 添加P value
    annotate("text",x = 0.1,y = 0.75,label = paste0("P = ",format(cor.test(pre_frac,true_frac)$p.value,digits=2,scientific = T)))+
    coord_cartesian(xlim = c(0,1),ylim = c(0,1))
  print(p)
  dev.off()
}
# CellDMC
rownames(est) <- colnames(m_train)
sex <- ifelse(pheno[match(colnames(m_train),pheno$Donor),]$Sex == "M",1,0)
age <- pheno[match(colnames(m_train),pheno$Donor),]$Age
pdf("age_dtr.pdf",width = 4,height = 4)
hist(age,main = "n =100")
dev.off()
pheno.train <- data.frame(sex,age)
rownames(pheno.train) <- colnames(m_train)
qc.train <- new("qc.o")
qc.train@m <- m_train
qc.train@s <- pheno.train
qc.train@ctf.o[[1]] <- as.matrix(est)

re <- CTSAC(qc.train)
dmct <- re[[1]]
sum(dmct[,2] != 0)

### check
load("~/Renv/simulation/blue_print_qc_nobatch.Rd")
library(broom)
library(qvalue)
idx <- which(qc.o@s$Cell == "CD4_N")
m_cd4 <- qc.o@m[,idx]
age <- qc.o@s$Age
# cd4 age-DMCs
cd4_md <- summary(lm(t(m_cd4)~age[idx]))
p_cd4 <- unlist(mclapply(cd4_md, function(x) glance(x)[[5]],mc.cores = 100))
t_cd4 <- unlist(mclapply(cd4_md, function(x) tidy(x)[[4]][2],mc.cores = 100))
q_cd4 <- qvalue(p_cd4, fdr.level = 0.05)
m_cd4_sig <- m_cd4[q_cd4$significant,]

fdr_threshold_mix <- min(abs(re[[2]]$CD4$t[sapply(re[[2]]$CD4$adjP, function(x) {x<=0.05})]))
fdr_threshold_sorted <- min(abs(t_cd4[sapply(q_cd4$qvalues, function(x) {x<=0.05})]))
# #fdr_threshold_mix <- 1.96
# 
# fet <- function(t1,t2,...){
#   n1 <- sum(t1 >= fdr_threshold_sorted & t2 >= fdr_threshold_mix)
#   n2 <- sum(t1 <= -fdr_threshold_sorted & t2 >= fdr_threshold_mix)
#   n3 <- sum(t1 <= -fdr_threshold_sorted & t2 <= -fdr_threshold_mix)
#   n4 <- sum(t1 >= fdr_threshold_sorted & t2 <= -fdr_threshold_mix)
#   m_f <- matrix(c(n1,n2,n4,n3),ncol = 2,byrow = TRUE)
#   re <- fisher.test(m_f,alternative = "greater")
#   return(list(re, m_f))
# }
# re_or <- fet(t_cd4, re[[2]]$CD4$t)
# OR <- format(re[[1]]$estimate,digits = 2)
# pval <- format(re[[1]]$p.value,digits = 1)
# m_number <- re[[2]]

pdf("DMCTs_cd4.pdf",width = 5,height = 5)
my_palette <- colorRampPalette(colors = c("white","#F2F7FC", "#66A7D5", "#4477B9", "#2A4D9C", "#1B3189", "#0C2167"))
smoothScatter(t_cd4, re[[2]]$CD4$t
              , xlab = ""
              ,colramp = my_palette
              ,nrpoints = 500
              ,ylab = "", 
              xlim = c(-12,12), ylim = c(-12,12), cex = 1.5
              , transformation = function(x) x^.25,main = "")
abline(v = 0, lty = 2, col = "black")
abline(h = 0, lty = 2, col = "black")
abline(v = fdr_threshold_sorted, lty = 2, col = "red")
abline(v = -fdr_threshold_sorted, lty = 2, col = "red")
abline(h = fdr_threshold_mix, lty = 2, col = "red")
abline(h = -fdr_threshold_mix, lty = 2, col = "red")
title(ylab="t(CD4:CellDMC)", mgp=c(2.4,1,0), cex.lab=1.5)
title(xlab="t(CD4:sorted)", mgp=c(2.4,1,0), cex.lab=1.5)
title(main="Mixture vs sorted cells", mgp=c(2.4,1,0), cex.lab=1.5)
dev.off()

save(qc.train,file = "qc_train_mode2.Rd")

############################ CD4 clock ##########################\
setwd("simulation/")
load("qc_train_mode2.Rd")
source("~/code/Dmathy/code/ELN_clock.R")
# re_clcok<-clock(qc.train,covs = NULL,seqs = seq(0.001,1,0.001))

library(EpiDISH)
m_tr <- qc.train@m
age_tr <- qc.train@s$age
ctf <- qc.train@ctf.o[[1]]
idx <- seq(1,ncol(m_tr),1)
n.idx <- idx
nP <- 10
###### chose bags ######
bags <- list()
for(p in 1:nP){
  bags[[p]] <- sample(idx,ncol(m_tr)/nP,replace=FALSE)
  idx <- setdiff(idx,bags[[p]])
}
#for each 9 partitions run CellDMC
CD4DMCT.lm <- list()
for(p in 1:nP){
  test.idx <- na.omit(bags[[p]])
  train.idx <- setdiff(n.idx,test.idx)
  result <- CTSAC(qc.train,idx = train.idx)
  idx <- which(result[[2]]$CD4$adjP <= 0.05)
  CD4DMCT.lm[[p]] <- rownames(result[[1]][,c(1,2)])[idx]
  cat("Got: ",length(CD4DMCT.lm[[p]])," CD4 age-DMCs")
}

p <- which.max(sapply(CD4DMCT.lm,length))
train.idx <- setdiff(n.idx,bags[[p]])
CD4.idx <- match(CD4DMCT.lm[[p]],rownames(m_tr))
seqs <- seq(0.5,1.5,0.001)
model_CD4 <- glmnet(t(m_tr[CD4.idx,train.idx]),y=age_tr[train.idx],alpha=1,standardize=TRUE,lambda=seqs)
predCD4.m <- t(predict(model_CD4,newx=t(m_tr[CD4.idx,bags[[p]]])))

rmseH.v <- vector();
for(li in 1:nrow(predCD4.m)){
  rmseH.v[li] <- sqrt(mean((predCD4.m[li,] - age_tr[bags[[p]]])^2))
}
la_CD4 <- rev(seqs)[which.min(rmseH.v)]
age_CD4 <- predCD4.m[which.min(rmseH.v),]
pdf("OptimizationCurve.pdf",width=5,height=3);
par(mar=c(4,4,2,1))
plot(rev(seqs),rmseH.v,ylab="RMSE",xlab="Lambda",type="l",main = "naiveCD4T")
# add a vertical line at the optimal lambda
abline(v=la_CD4,lty=2,col="red")
dev.off()


r <- cor(age_tr[bags[[p]]],age_CD4,method = "pearson")
# calculate the p-value
p_val <- cor.test(age_tr[bags[[p]]],age_CD4,method = "pearson")$p.value
# Calculate the Median Absolute Error between predicted age and real age
mae <- median(abs(age_tr[bags[[p]]] - age_CD4))
# draw scatter plot of predicted age vs. chronological age

pdf("pre_vs_chro.pdf",width=6,height=6)
plot(age_tr[bags[[p]]],age_CD4,pch = 16,col = "blue",ylab = "Predicted age (years)",xlab = "Chronological age (years)",main = "NaiveCD4 n = 10")
abline(0,1,col="red")
text(40,60,paste0("r = ",format(r, digits = 3),"\n","p = ",format(p_val, digits = 1),"\n","MedAE: ",format(mae, digits = 2),"years"))
dev.off()

clocks <- glmnet(t(m_tr[CD4.idx,]),y=age_tr,alpha=1,standardize=TRUE,lambda=la_CD4)
save(clocks,file = "CD4T_clcok.Rd")
############################# Validation #######################
setwd("simulation/")
load("CD4T_clcok.Rd")
load("testset.Rd")
source("~/code/Dmathy/code/ELN_clock.R")
load("blue_print_qc_nobatch.Rd")
clocks <- list(clocks)
names(clocks) <- "CD4_N"

idx <- match(colnames(m_control_1),qc.o@s$Donor)
qc.control <- new("qc.o")
qc.control@m <- m_control_1
qc.control@s <- qc.o@s[idx,]
colnames(qc.control@s)[9] <- "age"
age_control_1 <- validate(qc.control,clocks,gse = "control_1")$CD4_N[,1]

idx <- match(colnames(m_case_1),qc.o@s$Donor)
qc.case <- new("qc.o")
qc.case@m <- m_case_1
qc.case@s <- qc.o@s[idx,]
colnames(qc.case@s)[9] <- "age"
age_case_1 <- validate(qc.case,clocks,gse = "case_1")$CD4_N[,1]

idx <- match(colnames(m_control_1),qc.o@s$Donor)
qc.control <- new("qc.o")
qc.control@m <- m_control_2
qc.control@s <- qc.o@s[idx,]
colnames(qc.control@s)[9] <- "age"
age_control_2 <- validate(qc.control,clocks,gse = "control_2")$CD4_N[,1]

idx <- match(colnames(m_case_1),qc.o@s$Donor)
qc.case <- new("qc.o")
qc.case@m <- m_case_2
qc.case@s <- qc.o@s[idx,]
colnames(qc.case@s)[9] <- "age"
age_case_2 <- validate(qc.case,clocks,gse = "case_2")$CD4_N[,1]



df_mix <- list(age_control_1,age_case_1,age_control_2,age_case_2)

df_mix_t <- list(qc.control@s$age,qc.case@s$age,qc.control@s$age,qc.case@s$age)


#EAA
df_mix_eaa <- list()
for (i in 1:length(df_mix)){
  df_mix_eaa[[i]] <- df_mix[[i]] - df_mix_t[[i]]
}

# EAA_RES
resdu_1 <- resid(lm(c(age_control_1,age_case_1)~c(qc.control@s$age,qc.case@s$age)))
resdu_2 <- resid(lm(c(age_control_2,age_case_2)~c(qc.control@s$age,qc.case@s$age)))
df_mix_raa <-  list(resdu_1[1:20],resdu_1[21:39],resdu_2[1:20],resdu_2[21:39])


# IAA
library(EpiDISH)
est1 <- epidish(beta.m = qc.control@m, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
est2 <- epidish(beta.m = qc.case@m, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
est <- rbind(est1,est2)
est <- est[,c(3,5,6)]


resdu_1 <- resid(lm(c(age_control_1,age_case_1)~c(qc.control@s$age,qc.case@s$age) + est))
resdu_2 <- resid(lm(c(age_control_2,age_case_2)~c(qc.control@s$age,qc.case@s$age) + est))
df_mix_iaa <-  list(resdu_1[1:20],resdu_1[21:39],resdu_2[1:20],resdu_2[21:39])


df_mix_age <- c(df_mix_eaa,df_mix_raa,df_mix_iaa)
names(df_mix_age) <- c("co_mo1_eaa","ca_mo1_eaa","co_mo2_eaa","ca_mo2_eaa","co_mo1_raa","ca_mo1_raa","co_mo2_raa","ca_mo2_raa","co_mo1_iaa","ca_mo1_iaa","co_mo2_iaa","ca_mo2_iaa")
library(purrr)
library(dplyr)
data_df <- df_mix_age %>%
  map_df(~data.frame(group = deparse(substitute(.x)), value = .x), .id = 'group')

# 设置 group 列为因子，并指定其水平的顺序
data_df$group <- factor(data_df$group, levels = names(df_mix_age))

# 使用 ggplot2 绘制箱线图
# 初始化 ggplot
pdf("AA.pdf",width = 6,height = 4)
p <- ggplot(data_df, aes(x = group, y = value)) + 
  geom_boxplot()+
  # 展示点
geom_jitter(width = 0.2, size = 1, alpha = 0.6)

# 添加颜色区域和文本标注
p <- p + annotate("rect", 
                  xmin = 0.5, xmax = 4.5, ymin = -Inf, ymax = Inf, 
                  fill = "red", alpha = 0.1) +
  annotate("text", x = 2.5, y = 25, label = "EAA", vjust = -1)

p <- p + annotate("rect", 
                  xmin = 4.5, xmax = 8.5, ymin = -Inf, ymax = Inf, 
                  fill = "green", alpha = 0.1) +
  annotate("text", x = 6.5, y = 25, label = "EAA(residuals)", vjust = -1)

p <- p + annotate("rect", 
                  xmin = 8.5, xmax = 12.5, ymin = -Inf, ymax = Inf, 
                  fill = "blue", alpha = 0.1) +
  annotate("text", x = 10.5, y = 25, label = "IAA", vjust = -1)

# 调整x轴标签的角度和主题设置
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Random mixtures vs equal mixtures")+
  xlab("")+
  ylab("Age acceleration")

# 显示图形
print(p)
dev.off()

################################# CTS markers vs CTS DMCTs #############
setwd("simulation/")
load("qc_train_mode2.Rd")
source("~/code/Dmathy/code/ELN_clock.R")
load("CD4T_clcok.Rd")

cd4dmct <- clocks$beta@Dimnames[[1]]
clockCpGs <- clocks$beta@Dimnames[[1]][clocks$beta@i]

load("~/Renv/simulation/blue_print_qc_nobatch.Rd")
library(broom)
library(qvalue)
idx <- which(qc.o@s$Cell == "CD4_N")
m_cd4 <- qc.o@m[,idx]
age <- qc.o@s$Age
# CD4 age-DMCs
cd4_md <- summary(lm(t(m_cd4)~age[idx]))
p_cd4 <- unlist(mclapply(cd4_md, function(x) glance(x)[[5]],mc.cores = 100))
q_cd4 <- qvalue(p_cd4, fdr.level = 0.05)
m_cd4_sig <- m_cd4[q_cd4$significant,]
cd4adc <- rownames(m_cd4_sig)

source("~/code/Dmathy/code/ref_functions.R")
celltypes <- factor(qc.o@s$Cell)
hy_c<-dmc_look(qc.o@m,celltypes,bar = T)
cd4_marker <- hy_c$CD4_N


library(VennDiagram)
my_list <- list(cd4dmct = cd4dmct, cd4adc = cd4adc, cd4_marker = cd4_marker, clockCpGs = clockCpGs)

col<-rainbow(4)
venn.diagram(x=my_list,
             scaled = TRUE, # 根据比例显示大小
             alpha= 0.7, #透明度
             lwd=2,lty=1,col=col, #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             cex = 1, # 数字大小
             fontface = "bold",  # 字体粗细；加粗bold
             fill=col, # 填充色 配色https://www.58pic.com/
             category.names = c("DMCTs", "Age_DMCs", "CTS_markers", "ClockCpGs") , #标签名
             cat.dist = 0.03, # 标签距离圆圈的远近
             cat.pos = -180, # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             cat.cex = 1, #标签字体大小
             cat.fontface = "bold",  # 标签字体加粗
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             cat.default.pos = "text",  # 标签位置, outer内;text 外
             output=TRUE,
             filename='cpg_venn.png',# 文件保存
             imagetype="png",  # 类型（tiff png svg）
             resolution = 400,  # 分辨率
             compression = "lzw"# 压缩算法
)











