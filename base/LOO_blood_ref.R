load("~/data/infinium/MM285/Rdata/blood/raw_sort.Rd")
source("~/code/Dmathy/code/Damthy.R")
source("~/code/base/ref_functions.R")
############################################## PCA ###########################
idx_ea <- grep("EasySep", raw_sort@raw.s$CellType)
idx_mono <- which(raw_sort@raw.s$CellType == "Monocytes EasySep")
idx <- setdiff(idx_ea,idx_mono)

qc.o <- new("qc.o")
qc.o@m <- raw_sort@raw.m[,-idx]
pheno <- raw_sort@raw.s[-idx,]
pheno$tissue <- pheno$CellType
pheno$sex <- gsub("Sex: Male","1",pheno$sex)
pheno$sex <- gsub("Sex: Female","0",pheno$sex)
qc.o@s <- pheno
# # outliers
# remove monocytes
idx <- which(pheno$CellType == "Monocytes")
qc.o@m <- qc.o@m[,-idx]
qc.o@s <- qc.o@s[-idx,]
pheno <- qc.o@s
pheno$tissue <- factor(gsub(" EasySep","",qc.o@s$tissue))
# # outliers
# outlier <- c("GSM5587742")
# idx_out <- match(outlier,colnames(qc.o@m))
# qc.o@m <- qc.o@m[,-idx_out]
# qc.o@s <- qc.o@s[-idx_out,]
# combine MDSC
idx_MDSC <- which(qc.o@s$tissue == "CD11b Gr1 MDSC")
qc.o@s$tissue[idx_MDSC] <- "myeloid cells"
qc.o@pca.o <- list()
# combine Neutrophils with Mono
idx_neu <- which(qc.o@s$tissue == "Neutrophils")
qc.o@s$tissue[idx_neu] <- "myeloid cells"
idx_mono <- which(qc.o@s$tissue == "Monocytes EasySep")
qc.o@s$tissue[idx_mono] <- "myeloid cells"
qc.o@s$tissue <- gsub(" EasySep","",qc.o@s$tissue)
qc.o@pca.o <- PCA(qc.o,2)
qc.o@pca.o <- list()
###################################### markers numbers ################

############### look for outliers ##################
celltypes <- factor(make.names(qc.o@s$tissue))
hy_c<-dmc_look(qc.o@m,celltypes,bar = T)
nt <- 100
number.lv <- matrix(NA,nrow = nt,ncol = length(levels(celltypes)))
outs <- list()
for (seed in 1:nt){
print(seed)
out <- c()
for (i in 1:length(levels(pheno$tissue))) {
  idx <- which(pheno$tissue == levels(pheno$tissue)[i])
  out[i] <- sample(idx, 1)
}
m <- qc.o@m[,-out]
DMC <- list()
for (l in 1:length(hy_c)){
  name <- names(hy_c)[l]
  DMC[[name]] <- by_diff_max(hy_c[[l]],m,celltypes[-out],name,THRE = 0)
  #print(length(DMC[[name]]))
}
number.lv[seed,] <- sapply(DMC,length)
outs[[seed]] <- out
}
number.m<-scale(number.lv)
idx <- which.max(rowMeans(number.m))
idx_out<-outs[[idx]]
qc.o@m <- qc.o@m[,-idx_out]
qc.o@s <- qc.o@s[-idx_out,]
qc.o@pca.o <- PCA(qc.o,2)
##################      check             ################
celltypes <- factor(make.names(qc.o@s$tissue))
hy_c<-dmc_look(qc.o@m,celltypes,bar = T)
DMC <- list()
for (l in 1:length(hy_c)){
  name <- names(hy_c)[l]
  DMC[[name]] <- by_diff_max(hy_c[[l]],qc.o@m,celltypes,name,THRE = 0)
  print(length(DMC[[name]]))
}
dmc_alln <- sapply(DMC,length)
mars <- as.character(unlist(DMC))
########heat map  ######
library(ComplexHeatmap)
m_heat <- qc.o@m[mars,]
col <-  list(celltypes = c("B.Cells" = "blue", 
                           "CD8.T.Cells" = "red", 
                           "CD49b.NK.Cells" = "green", 
                           "myeloid.cells" = "yellow", 
                           "CD4.T.Cells" = "purple"))
pdf("markers_heat.pdf",width = 6,height = 6)
Heatmap(m_heat,
        use_raster = TRUE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        #热图树状图下方的注释加一行，注释的内容为细胞类型：celltypes,每个细胞类型颜色单独设置
        top_annotation = HeatmapAnnotation(celltypes = celltypes, col = col)
)
dev.off()
DMC <- list()
for (l in 1:length(hy_c)){
  name <- names(hy_c)[l]
  DMC[[name]] <- by_diff_max_hypo(hy_c[[l]],qc.o@m,celltypes,name,THRE = 0)
  print(length(DMC[[name]]))
}
dmc_hypon <- sapply(DMC,length)
dmc.m <- matrix(NA,ncol = 2,nrow = 5)
colnames(dmc.m) <- c("Hypo","Hyper")
rownames(dmc.m) <- names(dmc_hypon)
dmc.m[,1] <- dmc_hypon
dmc.m[,2] <- dmc_alln - dmc_hypon

dmc.m <- as.data.frame(dmc.m)
dmc.m$celltypes <- rownames(dmc.m)
dmc.m <- reshape2::melt(dmc.m,id.vars = "celltypes")
pdf("delta_result.pdf",width = 6,height = 6)
ggplot(dmc.m,aes(x = celltypes,y = value,fill = variable)) + 
  geom_bar(stat = "identity",position = "dodge") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "",y = "Number of DMCs",title = "Delta > 0") + 
  scale_fill_manual(values = c("Hypo" = "blue","Hyper" = "red"))+
  geom_text(aes(label = value),position = position_dodge(width = 0.9),vjust = -0.5)+
  guides(fill = guide_legend(title = NULL))
dev.off()
################## Start #####################################################
celltypes <- factor(make.names(qc.o@s$tissue))
set.seed(1234)
sequences <- matrix(NA,nrow = 100,ncol = 5)
for (i in 1:100) {
  seq <- runif(5)
  total <- sum(seq)
  seq_normalized <- seq / total
  sequences[i,] <- seq_normalized
}
lens <- c()
rmse <- matrix(NA,nrow = 100,ncol = 5)
for (seed in 1:100){
  print(seed)
  out <- c()
  for (i in 1:length(levels(celltypes))) {
    idx <- which(celltypes == levels(celltypes)[i])
    out[i] <- sample(idx, 1)
  }
  m <- qc.o@m[,-out]
  ### limma
  hy_c <- list()
  # compare one cell type with another three cell types together
  for (i in seq_along(levels(celltypes))){
    med <- as.character(celltypes[-out])
    other<-setdiff(levels(celltypes), levels(celltypes)[i])
    med[med %in% other] <- "other"
    med <- factor(med)
    design <- model.matrix(~0+med)
    colnames(design) <- levels(med)
    # fit the linear model
    fit <- lmFit(m, design)
    # create a contrast matrix for specific comparisons
    contrast_formula <- paste(levels(celltypes)[i], "- other", sep="")
    contMatrix <- makeContrasts(contrast_formula, levels = design)
    # fit the contrasts
    fit2 <- contrasts.fit(fit, contMatrix)
    fit2 <- eBayes(fit2)
    # look at the numbers of DM CpGs at FDR < 0.05
    #print(summary(decideTests(fit2)))
    # filter hyM
    hy_c[[levels(celltypes)[i]]] <- hy(fit2)
  }
  
  
  DMC <- list()
  for (l in 1:length(hy_c)){
    name <- names(hy_c)[l]
    DMC[[name]] <- by_diff_max(hy_c[[l]],m,celltypes[-out],name,THRE = 0.1)
    #print(length(DMC[[name]]))
  }
  len <- min(sapply(DMC,length))
  if (len < 50){
    next
  }
  lens[seed] <- len
  # Initialize rmse vector
  mars <- as.character(unlist(lapply(DMC,function(x) x[1:len])))
  ref <- m[mars,]
  ref.mx <- matrix(data = NA,nrow = length(mars),ncol = length(levels(celltypes)))
  rownames(ref.mx) <- mars
  colnames(ref.mx) <- levels(celltypes)
  for ( i in levels(celltypes)){
    idx <- which(celltypes[-out] == i)
    if (length(idx) > 1){
      med.c<-apply(ref[,idx],1,median)}
    else {
      med.c<-ref[,idx]
    }
    ref.mx[,i] <- med.c
  }
  weight <- rep(NA,times = nrow(ref.mx))
  ref_blood <- cbind(ref.mx,weight)
  ### test sample
  m_test <- qc.o@m[,out]
  
  predicted <- matrix(NA,nrow=nrow(sequences),ncol = 5)
  for (w in 1:nrow(sequences)){
    testset <- m_test %*% sequences[w,]
    rownames(testset) <- rownames(qc.o@m)
    predicted[w,]<-RPC(ref_blood,testset)
  }
  for (b in 1:ncol(sequences)){
    rmse[seed,b] <- sqrt(mean((predicted[,b]-sequences[,b])^2))
  }
}
rmse <- na.omit(rmse)
lens <- na.omit(lens)
rownames(rmse) <- lens
colnames(rmse) <- c("B.Cells","CD4.T.Cells","CD49b.NK.Cells","CD8.T.Cells","myeloid.cells")
idx <- which(lens > 150)
rmse <- rmse[-idx,]
rmse <- aggregate(. ~ rownames(rmse), data = as.data.frame(rmse), FUN = mean)
rownames(rmse) <- rmse$`rownames(rmse)`
rmse <- as.matrix(rmse)[,-1]
pdf("rmse.pdf", width = 15, height =  3)
par(mfrow = c(1,5))
for(i in 1:5){
  plot(rownames(rmse),rmse[,i], xlab = "Number of markers", ylab = "RMSE", main = colnames(rmse)[i])
  # 最小的点标红
  points(rownames(rmse)[which.min(rmse[,i])], rmse[which.min(rmse[,i]),i], col = "red", pch = 19)
  # x 轴刻度值是对应的行名
  
}
dev.off()

idx.lv <- apply(rmse,2,which.min)
number.c <- as.numeric(rownames(rmse)[idx.lv])
################################################# over and build ##########
out <- c()
for (i in 1:length(levels(celltypes))) {
  idx <- which(celltypes == levels(celltypes)[i])
  out[i] <- sample(idx, 1)
}
m <- qc.o@m[,-out]
### limma
hy_c <- list()
# compare one cell type with another three cell types together
for (i in seq_along(levels(celltypes))){
  med <- as.character(celltypes[-out])
  other<-setdiff(levels(celltypes), levels(celltypes)[i])
  med[med %in% other] <- "other"
  med <- factor(med)
  design <- model.matrix(~0+med)
  colnames(design) <- levels(med)
  # fit the linear model
  fit <- lmFit(m, design)
  # create a contrast matrix for specific comparisons
  contrast_formula <- paste(levels(celltypes)[i], "- other", sep="")
  contMatrix <- makeContrasts(contrast_formula, levels = design)
  # fit the contrasts
  fit2 <- contrasts.fit(fit, contMatrix)
  fit2 <- eBayes(fit2)
  # look at the numbers of DM CpGs at FDR < 0.05
  print(summary(decideTests(fit2)))
  # filter hyM
  hy_c[[levels(celltypes)[i]]] <- hy(fit2)
}


DMC <- list()
for (l in 1:length(hy_c)){
  name <- names(hy_c)[l]
  DMC[[name]] <- by_diff_max(hy_c[[l]],m,celltypes[-out],name,THRE = 0.1)
  print(length(DMC[[name]]))
  # pdf(paste0(name,".pdf"),width = 10, height = 10)
  # par(mfrow=c(2,2))
  # sapply(DMC[[name]][c(1,2,length(DMC[[name]])-1,length(DMC[[name]]))], function(cpg){
  #   plotCpg(m, cpg=cpg, pheno=celltypes[-out], ylab = "methylation fraction")
  # })
  # dev.off()
}
for (i in 1:length(DMC)){
  DMC[[i]] <- DMC[[i]][1:number.c[i]]
}
#DMC <- lapply(DMC,function(x) x[1:min(sapply(DMC,length))])
mars <- as.character(unlist(DMC))
ref <- m[mars,]
ref.mx <- matrix(data = NA,nrow = length(mars),ncol = length(levels(celltypes)))
rownames(ref.mx) <- mars
colnames(ref.mx) <- levels(celltypes)
for ( i in levels(celltypes)){
  idx <- which(celltypes[-out] == i)
  if (length(idx) > 1){
    med.c<-apply(ref[,idx],1,median)}
  else {
    med.c<-ref[,idx]
  }
  ref.mx[,i] <- med.c
}
weight <- rep(NA,times = nrow(ref.mx))
ref_blood <- cbind(ref.mx,weight)
colnames(ref_blood) <- c("B","CD4T","NK","CD8T","Myeloid","weight")
### test sample
m_test <- qc.o@m[,out]
set.seed(4567)
sequences <- list()
for (i in 1:100) {
  seq <- runif(5)
  total <- sum(seq)
  seq_normalized <- seq / total
  sequences[[i]] <- seq_normalized
}
predicted <- list()
for (w in 1:length(sequences)){
  testset <- m_test %*% sequences[[w]]
  rownames(testset) <- rownames(qc.o@m)
  predicted[[w]]<-RPC(ref_blood,testset)
}
cell <- as.character(colnames(predicted[[1]]))

for (l in 1:length(predicted[[1]])){
  pre_frac <- as.numeric(sapply(predicted,function(x) x[,l]))
  name <- cell[l]
  true_frac <- sapply(sequences,function(x) x[l])
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
    coord_cartesian(xlim = c(0,0.8),ylim = c(0,0.8))
  print(p)
  dev.off()
}
############################################ validation ######################
load("~/data/infinium/MM285/GSE201923/raw_bloodcell.Rd")
load("~/data/infinium/MM285/GSE201923/raw.Rd") 
load("~/Renv/mousebloodref/mice_ref.Rd")
load("~/Renv/mousebloodref/onlyhypo/miceref_hypo.Rd")
m <- raw_bc.o@raw.m
pheno <- raw_bc.o@raw.s
est.m<-RPC(ref_blood,m)
est.m <- as.data.frame(est.m/3)
est.m$celltype <- pheno$celltype

df_long <- reshape2::melt(est.m, id.vars = "celltype", variable.name = "CellType_ref", value.name = "Percentage")
ggplot(df_long, aes(x = celltype, y = Percentage, fill = CellType_ref)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = "Purified samples", y = "Proportion", fill = "estimation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


