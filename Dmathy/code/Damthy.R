## ----packages-------------------------------------------------------
library(ggplot2)
library(factoextra)
library(dplyr)
library(pryr)
library(magrittr)
library(broom)
library(parallel)
library(data.table)
library(pbapply)

## ----intialize------------------------------------------------------
#GSEnumber <- readline(prompt = "What's the GSE number you are analyzing? ")


## ----create class structure-----------------------------------------
setClass("Damthy", slots = list(name = "character"),
         prototype = list(name = "For infinuim"))
  setClass("rawD",
           contains = "Damthy",
           slots = list(rgset = "list",
                        m_pval = "data.frame",
                        GPL = "data.frame",
                        series = "data.frame"))
  setClass("pp",
           contains = "Damthy",
           slots = list(name = "character"), 
           prototype = list(name = "processed data"))
    setClass("raw.o",
             contains = "pp",
             slots = list(raw.m = "matrix",
                          raw.p = "matrix",
                          raw.s = "data.frame",
                          raw.g = "data.frame",
                          beta_dtr = "list"))
    setClass("qc.o",
             contains = "pp",
             slots = list(name = "character",
                          m = "matrix",
                          s = "data.frame",
                          beta_dtr = "list",
                          svd.o = "list",
                          pca.o = "list",
                          ctf.o = "list"),
             prototype = list(name = "quality control"))
    setClass("he.o",
               contains = "pp",
               slots = list(name = "character",
                            m = "matrix",
                            s = "data.frame",
                            svd.o = "list",
                            ctf.o = "list",
                            pca.o = "list"))


## ----instantiation function-----------------------------------------
f.rawD <- function(rgset,GPL,m_pval = data.frame(),...){
  if (ncol(m_pval) != 0){
    #print(paste0("you got the processed data with p values of ",GSEnumber))
    rawD <- new("rawD", m_pval = m_pval) 
    rawD@GPL <- GPL[match(rawD@m_pval[,1],GPL$Name),] #match GPL to m
  }else if(ncol(m_pval) == 0) {
    rgset <- rgset
    rawD <- new("rawD", rgset = list(rgset))
    rawD@GPL <- GPL #match GPL to m
  
    }
  return(rawD)
}

f.raw.o <- function(rawD,...){
  
  library(parallel)

  raw.o <- new("raw.o")
  if (length(rawD@m_pval) != 0){
  print("Tidy matrix with p values, you need to check the first column of the matrix (CpGs)!")
  # extract beta and P
  DNAm.beta <- rawD@m_pval[, seq(from=2, to=ncol(rawD@m_pval)-1, by= 2)]
  dim(DNAm.beta)
  DNAm.Pval <- rawD@m_pval[, seq(from=3, to=ncol(rawD@m_pval), by= 2)]
  colnames(DNAm.Pval) <- colnames(DNAm.beta)
  rownames(DNAm.Pval) <- rawD@m_pval[, 1]
  dim(DNAm.Pval)
  # p>0.01 = NA if illumina array
  DNAm.beta <- mcmapply(function(x,y) {ifelse(x>0.05,NA,y)}, DNAm.Pval, DNAm.beta, mc.cores = 100) #
  rownames(DNAm.beta) <- rawD@m_pval[, 1]
  cho <- readline(prompt = "do you want to use minfi to process it? T/F")
  if(cho){
  #require(minfiData)
  require(minfi)
  #-probe qc----------------------------------------------------------------------------------------------
  GRset =makeGenomicRatioSetFromMatrix(DNAm.beta) #  ,array ="IlluminaHumanMethylationEPIC",annotation = "ilm10b4.hg19"
  # delete sex chromosome CpGs ###
  annotation <- getAnnotation(GRset)
  sex_probe <- rownames(annotation)[annotation$chr %in% c("chrX","chrY")]
  keep <- !(featureNames(GRset) %in% sex_probe)
  GRset <- GRset[keep,]
  m <- getBeta(GRset)
  p.m <- as.matrix(DNAm.Pval[rownames(m),]) # filter CpGs with SNPs
  if(any(is.na(p.m))){
    p.m[is.na(p.m)] <- 1
  }
  GRset <- dropLociWithSnps(GRset,snps = c("CpG", "SBE"))
  }else{
    # delete sex chromosome CpGs ###
    keep <- !(rawD@GPL$chr %in% c("chrY","chrX"))
    m <- DNAm.beta[keep,]
    p.m <- as.matrix(DNAm.Pval[rownames(m),])
  }
  }
  if (length(rawD@rgset) != 0){
    rgset <- rawD@rgset[[1]]
    getManifest(rgset)
    detP <- detectionP(rgset)
    MSet <- preprocessRaw(rgset)
    GMset <- mapToGenome(MSet)
    GRset <- ratioConvert(GMset, what = "beta", keepCN = TRUE)
    # filter CpGs with SNPs
    GRset <- dropLociWithSnps(GRset,snps = c("CpG", "SBE"))
    # delete sex chromosome CpGs
    annotation <- getAnnotation(GRset)
    sex_probe <- rownames(annotation)[annotation$chr %in% c("chrX","chrY")]
    keep <- !(featureNames(GRset) %in% sex_probe)
    GRset <- GRset[keep,]
    m_r <- getBeta(GRset)
    # p>0.01 = NA
    p.m <- detP[rownames(m_r),]
    m <- matrix(mcmapply(function(x,y) {ifelse(x>0.01,NA,y)}, p.m, m_r, mc.cores = 50), nrow = nrow(m_r), ncol = ncol(m_r))
    rownames(m) <- rownames(m_r)
    colnames(m) <- colnames(m_r)
 }
  # over
  raw.o@raw.m <- m
  raw.o@raw.p <- p.m
  raw.o@raw.g <- rawD@GPL[match(rownames(m),rawD@GPL$Name),]

  gc()
  return(raw.o)
  
}
#coverage and imputation----------------------------------------------------------------------------
imp <- function(raw.o,cutoff=0.01){
  m <- raw.o@raw.m
  s <- raw.o@raw.s
  if (length(raw.o@raw.p) != 0){
    p.m <- raw.o@raw.p
    
    #-coverage filter-----------------------------------------------------------------------------------------------
    #samples
    tf <- readline(prompt = "Filter samples ?(T/F)")
    tf <- ifelse(tf == "F", FALSE, TRUE)
    if(tf){
    thre <- readline(prompt = "what's the threshold for sample ?")
    idx <- pbapply(p.m,2,function(x) {sum(x < cutoff)/length(x) >= thre})
    m <- m[,idx]
    p.m <- p.m[,idx]
    cat("Getting: ",sum(idx),"samples","\n")
    raw.o@raw.m <- m
    raw.o@raw.p <- p.m
    raw.o@raw.s <- s[idx,]
    coverage(raw.o,cutoff)
    }
    #probes
    thre <- readline(prompt = "what's the threshold for probes ?")
    idx <- pbapply(p.m,1,function(x) {sum(x < cutoff)/length(x) >= thre})
    p.m <- p.m[which(idx),]
    m <- m[which(idx),]
    cat("Getting: ",sum(idx),"probes","\n")
    raw.o@raw.p <- p.m
    
  }
  else if(any(is.na(m))){
  # no p-values but has NA in beta matrix
  tf <- readline(prompt = "Filter samples ?(T/F)")
  tf <- ifelse(tf == "F", FALSE, TRUE)
  if(tf){
    thre <- readline(prompt = "what's the threshold for sample ?")
    idx <- pbapply(m,2,function(x) {sum(!is.na(x))/length(x) >= thre})
    m <- m[,idx]
    cat("Getting: ",sum(idx),"samples","\n")
    raw.o@raw.m <- m
    raw.o@raw.s <- s[idx,]
    coverage(raw.o)
  } 
  thre <- readline(prompt = "what's the threshold for probes ?")
  idx <- pbapply(m,1,function(x) {sum(!is.na(x))/length(x) >= thre})
  m <- m[idx,]
  cat("Getting: ",sum(idx),"probes","\n")
  raw.o@raw.m <- m
  }
  tf <- readline(prompt = "Can I carry on T/F ?(T/F)")
  tf <- ifelse(tf == "F", FALSE, TRUE)
  #-imputation--------------------------------------------------------------------------------------------
  if (tf){
    if (any(is.na(m) & ncol(m) >= 100)){
      library(impute)
      print("starting impute")
      m <- impute.knn(m,k=5)$data
    }else{
      print("No NA")
    }
  }else{
    stop("over")
  }
  raw.o@raw.m <- m
  raw.o@raw.g <- raw.o@raw.g[match(rownames(m),raw.o@raw.g$Name),]
  print("OK")
  return(raw.o)
}

f.qc.o <- function(raw.o,...){
qc.o <- new("qc.o")
#bmiq-------------------------------------------------
choice <- readline(prompt = "Do you need to adjust type 2 probe bias? y/n")
if (choice == "y"){
  result <-bmiq(raw.o)
  qc.o@m <- result[[1]]
  design.v <- result[[2]]
}else if (choice == "n"){
  design.v <- tryCatch({
    raw.o@raw.g$Infinium_Design_Type %>% gsub("II",2,.) %>% gsub("I",1,.) %>% as.numeric()
  }, error = function(e) { 
    print("Check the GPL file:GPL$Infinium_Design_Type,II,I") 
  })
  qc.o@m <-raw.o@raw.m
}else{
  stop("please input a right word:y/n")
}
qc.o@beta_dtr <- beta_dtp(qc.o,design.v,raw.o)
pdf("BMIQ.pdf", width=6,height=4)
for (i in 1:length(qc.o@beta_dtr)){
  plot(qc.o@beta_dtr[[i]])}
dev.off()
#average dups-------------------------------------------
choice <- readline(prompt = "Remove duplicates? 0/1 ")
if (choice == 1){
duplicates <- table(factor(rawData@series$`individual id:ch1`))[table(factor(rawData@series$`individual id:ch1`)) > 1]
if (length(duplicates)>=1){
  print("Dealing dups ...")
nm <- qc.o@m
ns <- raw.o@raw.s
for (i in names(duplicates)){
  idx <- which(rawData@series$`individual id:ch1` %in% i)
  mean_row <- rowMeans(qc.o@m[,idx])
  nm <- cbind(nm,mean_row)
  colnames(nm)[ncol(nm)] <- colnames(qc.o@m)[idx[1]]
  ns <- rbind(ns,raw.o@raw.s[idx[1],])
}
ddx <- which(rawData@series$`individual id:ch1` %in% names(duplicates))
qc.o@m <- nm[,-ddx]
qc.o@s <- ns[-ddx,]
}
}else{
  qc.o@s <- raw.o@raw.s
}
return(qc.o)
}

f.he.o <- function(qc.o,...){
  he.o <- new("he.o")
  idx <- which(qc.o@s$disease == 0)
  he.o@m <- qc.o@m[,idx]
  he.o@s <- qc.o@s[idx,-which(colnames(qc.o@s) == "disease")]
  return(he.o)
}

## ----coverage of probe--------------------------------------------------
setGeneric("coverage",function(obj,...){
  standardGeneric("coverage")
})
setMethod("coverage","raw.o",function(obj,cutoff=0.01,...){
  library(ggplot2)
  library(pbapply)
  library(gridExtra)
 
  if (length(obj@raw.p) == 0){
    cg <- pbapply(obj@raw.m,1,function(x) {sum(!is.na(x))/length(x)})   # coverage in each probe
    sg <- pbapply(obj@raw.m,2,function(x) {sum(!is.na(x))/length(x)})   # coverage in each sample
  }else{
    cg <- pbapply(obj@raw.p,1,function(x) {sum(x < cutoff)/length(x)})
    sg <- pbapply(obj@raw.p,2,function(x) {sum(x < cutoff)/length(x)})
  }
  pc.df <- data.frame(rownames(obj@raw.m),cg)
  colnames(pc.df)<-c("probe","coverage")
  ps.df <- data.frame(colnames(obj@raw.m),sg)
  colnames(ps.df)<-c("sample","coverage")
  plot1 <- ggplot(pc.df, aes(x = probe, y = coverage)) + 
    geom_point() +
    theme(axis.text.x = element_blank())+
    xlab("Probe") + ylab("Coverage")+
    ggtitle("Probe Coverage Plot")
  plot2 <- ggplot(ps.df, aes(x = sample, y = coverage)) + 
    geom_point() +
    theme(axis.text.x = element_blank())+
    xlab("Sample") + ylab("Coverage") + 
    ggtitle("Sample Coverage Plot")
  my_grid <- grid.arrange(plot1, plot2, layout_matrix=rbind(c(1,1),c(2,2)))
  ggsave("coverage.pdf",my_grid,width = 8,height = 8)
  
  })


## ----beta_dtp-----------------------------------------------------------
setGeneric("beta_dtp",function(obj,...){
  standardGeneric("beta_dtp")
})
setMethod("beta_dtp","raw.o",function(obj,...){
  plots <- list()
  for (i in 1:ncol(obj@raw.m)){
  density2 <- density(obj@raw.m[which(obj@raw.g$Infinium_Design_Type == "II"),i])
  density1 <- density(obj@raw.m[which(obj@raw.g$Infinium_Design_Type == "I"),i])
  density_data <- data.frame(
  beta1 = density1$x,
  beta2 = density2$x,
  density1 = density1$y,
  density2 = density2$y
  )
  fplot <- ggplot(density_data, aes(x = beta1, y = density1)) + 
  geom_line(aes(x = beta1, y = density1, color = "Type I")) + 
  geom_line(aes(x = beta2, y = density2, color = "Type II")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "Beta", y = "Density", title = "Density plot of two types of probe") +
  scale_color_manual(name = "Probe Type", values = c("Type I" = "blue", "Type II" = "red"), labels = c("Type I", "Type II")) +
  theme(legend.position = c(0.95, 0.95),
        legend.justification = c("right", "top"),
        legend.margin = margin(-10, 0, 0, -10),
        legend.box.margin = margin(0, 0, 10, 10),
        legend.direction = "vertical",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
  plots[[i]] <- fplot
  }
  return(plots)
})

setMethod("beta_dtp","qc.o",function(obj,design.v,raw.o,...){
  plots <- list()
  for (i in 1:ncol(obj@m)){
  density2 <- density(obj@m[which(design.v == 2),i])
  density_data <- data.frame(
    beta = density2$x,
    density = density2$y
  )
  fplot <- raw.o@beta_dtr[[i]] + 
    geom_line(data = density_data,aes(x = beta, y = density, color = "Type II-BMIQ")) + 
    scale_color_manual(name = "Probe Type", values = c("Type I" = "blue", "Type II" = "red", "Type II-BMIQ" = "black"), labels = c("Type I", "Type II", "Type II-BMIQ")) + labs(x = "Beta", y = "Density", title = paste("Density plot of two types of probe of sample",i  ,"after BMIQ",sep = " "))
  plots[[i]] <- fplot
  }
  return(plots)
})
#1# ----Dealing with the duplicates------------------------------------
setGeneric("Del_dup",function(obj,...){
  standardGeneric("Del_dup")
})
setMethod("Del_dup","raw.o",function(obj,...){
  
  c <- readline(prompt = "Have you dealt the duplicates? T/F")
  library(broom)
  #library(isva)
  # outliers
  # identify_outliers <- function(x) {
  #   # 计算所有数的绝对偏差
  #   abs_dev <- abs(x - median(x))
  #   # 计算绝对中位差
  #   mad <- median(abs_dev)
  #   # 根据绝对中位差，将可能的离群值定义为离群区域的2倍以上的值
  #   # 输出所有可能的离群值
  #   outliers <- x[abs_dev / mad > 2]
  #   return(outliers)
  # }
  #----------------------------------
  #m <- obj@raw.m - rowMeans(obj@raw.m)
  #print("Estimating dimensionality")
  # n <- EstDimRMT(m)$dim
  # print("Run SVD")
  # svd.m <- svd(m)
  # v_top <- svd.m$v[,1:n]
  # #  rownames(v_top) <- colnames(obj@m)
  # cor_m <- as.dist(1-cor(t(v_top)))
  # clu.o <- hclust(cor_m, method = "complete")
  nm <- obj@raw.m
  ns <- obj@raw.s
  duplicates <- table(factor(obj@raw.s$`individual id:ch1`))[table(factor(obj@raw.s$`individual id:ch1`)) > 1]
  
  # par(mfrow=c(1,length(duplicates)),mar=c(4,4,4,4))
  # pdf("dup_clust.pdf", width=10,height=5)
  for (i in names(duplicates)){
    idx <- which(obj@raw.s$`individual id:ch1` %in% i)
    #pos <- as.numeric(clu.o$order) %in% idx
    #delete and average----------
    # if (length(idx) == 2){
    #   dis <- abs(match(idx[1],clu.o$order)-match(idx[2],clu.o$order))
    #   if(dis == 1){
        mean_col <- colMeans(obj@raw.m[,idx])
        nm <- cbind(nm,mean_col)
        colnames(nm)[ncol(nm)] <- rownames(obj@raw.m)[idx[1]]
        ns <- rbind(ns,obj@raw.s[idx[1],])
      # }else{
      #   s1<-sum(obj@raw.p[,ids[1]]<0.01)
      #   s2<-sum(obj@raw.p[,ids[2]]<0.01)
      #   if (s1>s2){
      #     nm <- cbind(nm,obj@raw.m[,ids[1]])
      #     ns <- rbind(ns,obj@raw.s[idx[1],])
      #   }else {
      #     nm <- cbind(nm,obj@raw.m[,ids[2]])
      #     ns <- rbind(ns,obj@raw.s[idx[2],])
      #   }
      # }
    # }else {
    #   a <- identify_outliers(sort(match(idx,clu.o$order)))
    #   dx <- clu.o$order[a]
    #   px <- setdiff(idx,dx)
    #   mean_col <- colMeans(obj@raw.m[,px])
    #   nm <- cbind(nm,mean_col)
    #   colnames(nm)[ncol(nm)] <- rownames(obj@raw.m)[px[1]]
    #   ns <- rbind(ns,obj@raw.s[px[1],])
    # }
    #---------------------------------------
    # labels <- rep("", length(clu.o$order))
    # labels[pos] <- as.character(idx)
    # dend <- as.dendrogram(clu.o)
    # dend <- color_branches(dend,k=13)
    # dend %>% set("labels", labels) %>%  set("labels_cex",0.7) %>% plot
    
  }
  ddx <- which(obj@raw.s$`individual id:ch1` %in% names(duplicates))
  nm <- nm[,-ddx]
  ns <- ns[-ddx,]
  #dev.off()
  return(list(nm,ns))
  gc()
})

#2# ----bmiq-----------------------------------------------------------
setGeneric("bmiq",function(obj,...){
  standardGeneric("bmiq")
})
setMethod("bmiq","raw.o",function(obj,...){
source("~/code/packages/BMIQ_1.6.R")
tmp.m <- obj@raw.m
num_cpus <- detectCores()-2
tmp.lv <- split(tmp.m, col(tmp.m))
design.v <- tryCatch({
  obj@raw.g$Infinium_Design_Type %>% gsub("II",2,.) %>% gsub("I",1,.) %>% as.numeric()
}, error = function(e) { 
  print("Check the GPL file:GPL$Infinium_Design_Type,II,I") 
})
print("Start BMIQ")
# tmp.m <- mclapply(tmp.lv, function(x) {
#   bmiq.o <- BMIQ(x, design.v)
#   return(bmiq.o$nbeta)
# }, mc.cores = num_cpus)
# print(dim(tmp.m))
tmp.m <-tryCatch({
result_list <- mclapply(tmp.lv, function(x) {
  bmiq.o <- BMIQ(x, design.v)
  return(bmiq.o$nbeta)
}, mc.cores = num_cpus)
do.call(cbind, result_list)
}, warning = function(w) {
  message("Warning:", w$message)
  return(NULL)
}, error = function(e) {
  message("Error:", e$message)
  return(NULL)
}
)
# for(s in 1:ncol(tmp.m)){
#   beta.v <- obj@raw.m[,s]
#   bmiq.o <- BMIQ(beta.v, design.v, sampleID=s)
#   tmp.m[,s] <- bmiq.o$nbeta
#   print(paste("Done BMIQ for sample ",s,sep=""))
# }
rownames(tmp.m) <- rownames(obj@raw.m)
colnames(tmp.m) <- colnames(obj@raw.m)
print(dim(tmp.m))
print("BMIQ is over")
lv <- list(tmp.m,design.v)
return(lv)
})

#3# ----CTF------------------------------------------------------------
setGeneric("CTF",function(obj,...){
  standardGeneric("CTF")
})
setMethod("CTF","pp",function(obj,ref,wth=0.4,useW=TRUE,type="850k",...){
  if (readline(prompt = "EpiSCORE or EpiDISH? T/F")){
    library(EpiSCORE)
    avSIM.m <- constAvBetaTSS(obj@m, type=type)
    estF.o <- wRPC(avSIM.m, ref=ref, useW=useW, wth=wth, maxit=500)
    fplot <- boxplot(estF.o$estF, main = "Estimation of cell-type fractions ", xlab = "Cell types", ylab = "Fraction")
    lv <- list(estF.o$estF,fplot)
    return(lv)
  } else{
    library(EpiDISH) 
    est <- epidish(beta.m = obj@m, ref.m = ref, method = "RPC")$estF
    fplot <- boxplot(est, main = "Estimation of cell-type fractions ", xlab = "Cell types", ylab = "Fraction")
    lv <- list(est,fplot)
    return(lv)
  }
})



#4# ----svd_lm---------------------------------------------------------
setGeneric("lm_svd",function(obj,...){
  standardGeneric("lm_svd")
})
setMethod("lm_svd","pp",function(obj,...){
  if (readline(prompt = "Are you ready to rum lm with a neat tibble of sample features? T/F")){
  library(broom)
  library(isva)
  m <- obj@m - rowMeans(obj@m)
  print("Estimating dimensionality")
  n <- EstDimRMT(m)$dim
  sam_var <- bind_cols(obj@ctf.o[[1]],obj@s)
  print("Run SVD")
  if (length(obj@svd.o)==0){
  svd.m <- svd(m)
  }else {
  svd.m <- obj@svd.o[[2]] 
  }
  v_top <- svd.m$v[,1:n]
  pval.m <- matrix(ncol = n,nrow = ncol(sam_var))
  print("Run lm")
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
  list <- list(pval.m,svd.m)
  return(list)
  }
})


#5# ----p_heatmap------------------------------------------------------
setGeneric("p_h",function(obj,...){
  standardGeneric("p_h")
})
setMethod("p_h","pp",function(obj,pic = c(24,2),...){
  fV.v <- obj@svd.o[[2]]$d^2/sum(obj@svd.o[[2]]$d^2)
  plot.new()
  name <- paste0("SVDsummary_",obj@name)
  pdf(paste0(name,".pdf"),width=8,height=8)
  layout(matrix(1:2,nrow=2),heights=c(1,1))
  par(mar=c(4,6,2,1))
  plot(fV.v[1:ncol(obj@svd.o[[1]])],pch=23,type="b",col="red",ylab="fracV",xlab="Top-PC",main="")
  my_colors <- c("#3C2E2B", "#C92742", "#ECF659", "#FAC7EE", "#FFFFFF")
  # Define breaks for color legend
  my_breaks <- c(-300,-50,-15,-5,log10(0.05),0); #EPiSCORE  chosed by plot(m)
  # breaks.v <- c(-82,-25,-5,-3,log10(0.05),0); #EpiDISH
  
  # Create heat map without clustering
#  try(heatmap.2(log10(obj@svd.o[[1]]), dendrogram = "none", Rowv = FALSE, Colv = FALSE, col = my_colors, breaks = my_breaks, key = FALSE, key.title = NA, margins = c(5,10),trace = "none",labRow = "",labCol = ""))
  par(mar = c(4, pic[1], 2, pic[2]))
  image(x=1:ncol(obj@svd.o[[1]]),y=1:nrow(obj@svd.o[[1]]),z=log10(t(obj@svd.o[[1]])),col=my_colors,breaks=my_breaks,xlab="",ylab="",axes=FALSE,asp=2);
  axis(1,at=1:ncol(obj@svd.o[[1]]),labels=paste("PC-",1:ncol(obj@svd.o[[1]]),sep=""),las=2)
  axis(2,at=1:nrow(obj@svd.o[[1]]),labels=rownames(obj@svd.o[[1]]),las=2)
  # Add legend
  legend("bottomright", legend = c("p<1e-50", "p<1e-15", "p<1e-5", "p<0.05", "p=ns"), fill = my_colors, bty = "o",
         cex = 0.4,
         pt.cex = 1)
  dev.off()
})


#6# ------PCA ----------------------------------------------------------

setGeneric("PCA",function(obj,...){
  standardGeneric("PCA")
})

setMethod("PCA", "pp", function(obj, s, ...) {
  library(RColorBrewer)
  if (length(obj@pca.o)==0){
#    matrix <- scale(t(obj@m), center = TRUE, scale = TRUE)
    matrix <- t(obj@m)
    pca_o <- prcomp(matrix, scale = FALSE)
    pca_data <- as.data.frame(pca_o$x[, 1:s])
    importance <- summary(pca_o)$importance
  } else {
    pca_o <- obj@pca.o[[1]]
    pca_data <- as.data.frame(obj@pca.o[[1]]$x[, 1:s])
    importance <- summary(pca_o)$importance
  }
  n <- readline(prompt = "Which feature do you want to mark?")
  idx <- NULL
  fill_values <- NULL
  test <- 0
  par(mfrow=c(1,(s-1)),mar=c(4,4,4,4))
  name <- paste("PCAsummary_", obj@name, "_", n, sep = "")
  pdf(paste0(name, ".pdf"), width=6,height=6)
  if (n == "sex") {
    idx <- which(obj@s$sex == 0)
    fill_values <- c("male" = "gray", "Female" = "brown")
  } else if (n == "disease") {
    idx <- which(obj@s$disease == 0)
    fill_values <- c("Bad" = "gray", "Healthy" = "blue")
  } else if(n == "age"){
    test <- 1
    idx1 <- which(obj@s$age <= 40)
    idx2 <- which(obj@s$age > 40 & obj@s$age <=60)
    fill_values <- c("Old" = "gray", "Young" = "red","Middle" = "blue" )
    for (i in 1:(ncol(pca_data)-1)) {
      
      g <- ggplot(data = pca_data, aes_string(x = paste0("PC", i), y = paste0("PC", i+1))) +
        geom_point(shape = 21, aes(fill = names(fill_values)[1]), color = "black", size = 2) +
        geom_point(data = pca_data[idx1, ], shape = 21, aes(fill = names(fill_values)[2]), color = "black", size = 2) +
        geom_point(data = pca_data[idx2, ], shape = 21, aes(fill = names(fill_values)[3]), color = "black", size = 2) +
        labs(x = paste0(colnames(pca_data)[i],",",format(100*importance[2,i],digits=1),"%"), y = paste0(colnames(pca_data)[i],",",format(100*importance[2,i+1],digits=1),"%"), fill = "") +
        scale_fill_manual(values = fill_values) +
        theme_classic() +
        theme(legend.position = "top", legend.justification = "right")
      
      plot(g)
    }
   }else if (n == "tissue") {
    test <- 1
    tissues <- levels(factor(obj@s$tissue))
    idx_lv <- list()
    colors <- rainbow(length(tissues))
    fill_values <- c()
    for( i in 1:length(tissues)){
      idx_lv[[i]] <- which(obj@s$tissue == tissues[i])
      fill_values[tissues[i]] <- colors[i]
    }
    pca_data$tissue_color <- obj@s$tissue
    number_of_colors <- length(unique(pca_data$tissue_color))
    distinct_colors <- brewer.pal(number_of_colors, "Set1")
    for (i in 1:(ncol(pca_data)-2)) {
      
      g <- ggplot(data = pca_data, aes(x = .data[[paste0("PC", i)]], y = .data[[paste0("PC", i+1)]])) +
        geom_point(aes(fill = tissue_color), shape = 21, color = "black", size = 4) +
        scale_fill_manual(values = distinct_colors) + # Use the distinct color palette
        labs(x = paste0(colnames(pca_data)[i],", ",format(100*importance[2,i],digits=1),"%"), 
             y = paste0(colnames(pca_data)[i+1],", ",format(100*importance[2,i+1],digits=1),"%"),
             title = paste0("n = ", ncol(obj@m)),
             fill = "") +
        theme_classic() +
        theme(legend.position = "top", 
              legend.justification = "right",
              axis.text.x = element_text(size=rel(1.5)), # Increase x-axis text size
              axis.text.y = element_text(size=rel(1.5)), # Increase y-axis text size
              axis.title.x = element_text(size=rel(1.5)), # Increase x-axis title size
              axis.title.y = element_text(size=rel(1.5))) # Increase y-axis title size
      
      plot(g)
      }
  }else {
    stop("Please input a valid feature!")
    test = 1
  }
if (test == 0){
  
  
  for (i in 1:(ncol(pca_data)-1)) {
    
    g <- ggplot(data = pca_data, aes_string(x = paste0("PC", i), y = paste0("PC", i+1))) +
      geom_point(shape = 21, aes(fill = names(fill_values)[1]), color = "black", size = 2) +
      geom_point(data = pca_data[idx, ], shape = 21, aes(fill = names(fill_values)[2]), color = "black", size = 2) +
      labs(x = paste0(colnames(pca_data)[i],",",format(100*importance[2,i],digits=1),"%"), y = paste0(colnames(pca_data)[i],",",format(100*importance[2,i+1],digits=1),"%"), fill = "") +
      scale_fill_manual(values = fill_values) +
      theme_classic() +
      theme(legend.position = "top", legend.justification = "right")
    
    plot(g)
  }
}
  dev.off()
  list(pca_o)
  })  




































































