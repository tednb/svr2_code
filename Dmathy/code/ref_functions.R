library(minfi)
library(ggplot2)


# deconvolution method ----------------------------------------------------
RPC <- function(mice_ref,mf.m,no=1,maxit){
  library(MASS)
  common.v <- intersect(rownames(mice_ref),rownames(mf.m))
  map.idx <- match(common.v,rownames(mf.m))
  rep.idx <- match(common.v,rownames(mice_ref))
  data.m <- as.matrix(mf.m[map.idx,])
  ref.m <- mice_ref[rep.idx,-ncol(mice_ref)]
  est.m <- matrix(data = NA,nrow=ncol(data.m),ncol=ncol(ref.m)) 
  colnames(est.m) <- colnames(ref.m)
  rownames(est.m) <- colnames(data.m)
  sum_squared_residuals <- vector("numeric", ncol(data.m))
  for (s in 1:ncol(data.m)) {
    rlm.o <- rlm(data.m[, s] ~ ref.m, maxit = maxit)
    if(nrow(summary(rlm.o)$coef) < ncol(ref.m)+1){
      coef.v <- rep(0,ncol(ref.m))
    }else{
      coef.v <- summary(rlm.o)$coef[2:(ncol(ref.m)+1),1];
    }
    if (no == 1){
      coef.v[which(coef.v<0)] <- 0;
      total <- sum(coef.v);
      coef.v <- coef.v/total; 
    }
    est.m[s,] <- coef.v;
    residuals <- rlm.o$residuals
    #sum_squared_residuals[s] <- sum(residuals^2)  # 平均平方残差
  }
  # pdf("liver_fraction.pdf",width = 6, height = 4)
  # boxplot(est.m, main = "Estimation of liver fractions", xlab = "Cell types", ylab = "Fraction")
  # dev.off()
  #mean(sum_squared_residuals)
  est.m
}

QP <- function(mice_ref,mf.m,choice = "<="){
  mice_ref <- na.omit(mice_ref)
  common.v <- intersect(rownames(mice_ref),rownames(mf.m))
  map.idx <- match(common.v,rownames(mf.m))
  rep.idx <- match(common.v,rownames(mice_ref))
  data.m <- mf.m[map.idx,]
  ref.m <- mice_ref[rep.idx,-ncol(mice_ref)]
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
    
    
    residual_sum_of_squares_liver[s] <- sum(residuals^2)/nrow(data.m)
    #print(sum(residuals^2))
  }
  return(residual_sum_of_squares_liver)
  #return(est.m)
}

library(EpiDISH)
hepidish_bm <- function(beta.m,refMouseBoneM.list){
  frac1 <- hepidish(beta.m = beta.m,ref1.m = refMouseBoneM.list$refLayer1.m, ref2.m = refMouseBoneM.list$Layer2_T.m, h.CT.idx = 5, method = "RPC", maxit = 500)
  frac2 <- hepidish(beta.m = beta.m,ref1.m = refMouseBoneM.list$refLayer1.m, ref2.m = refMouseBoneM.list$Layer2_LC.m, h.CT.idx = 2, method = "RPC", maxit = 500)
  frac3 <- hepidish(beta.m = beta.m,ref1.m = refMouseBoneM.list$refLayer1.m, ref2.m = refMouseBoneM.list$Layer2_MN.m, h.CT.idx = 4, method = "RPC", maxit = 500)
  frac.m <- cbind(frac1[,c(1,3,5,6)],frac2[,c(5,6)],frac3[,c(5,6)])
  return(frac.m)
}

bmRPC <- function(ref1,ref2,ref3,ref4,mat){ # ncol(ref1) < ncol(ref2)
  est1 <- as.data.frame(RPC(ref1,mat))
  est2 <- as.data.frame(RPC(ref2,mat))
  est3 <- as.data.frame(RPC(ref3,mat))
  est4 <- as.data.frame(RPC(ref4,mat))
  est1$CGP <-  est1$LC * est2$CGP
  est1$LSK <-  est1$LC * est2$LSK
  
  est1$CD8 <-  est1$T * est4$CD8
  est1$CD4 <-  est1$T * est4$CD4
  
  est1$Neu <-  est1$MN * est3$Neu
  est1$Mono <-  est1$MN * est3$Mono
  
  est1$LC <- NULL
  est1$MN <- NULL
  est1$T <- NULL
  est1 <- as.matrix(est1)
  # est.m<-apply(est2,1,function(x) {
  #   x[which(x<0)] <- 0;
  #   total <- sum(x);
  #   x <- x/total; 
  # })
  # rownames(est.m) <- colnames(est1)
  return(est1)
}
# limma + delta -----------------------------------------------------------
library(limma)
## Extract hyM genes
hy <- function(fit){
  marker <- list()
  DMPs <- topTable(fit, num = Inf, coef = 1, sort.by = "p")
  DMPs_hy <- DMPs[which(DMPs$adj.P.Val <= 0.05),]
  markers <- rownames(DMPs_hy)
  #print(length(markers))
  return(markers)
}
## filter sig markers

by_diff_max <- function(mark,ref.m,CT,cell,THRE){
  deltas <- function(x,idx_other,idx_cell){
    if(mean(x[idx_other]) > mean(x[idx_cell])){
      delta <- min(x[idx_other])-max(x[idx_cell])
    }
    else{
      delta <- min(x[idx_cell])-max(x[idx_other])
    }
    return(delta)
  }
  m <- ref.m[mark,]
  idx_cell <- which(CT == cell)
  idx_other <- setdiff(1:length(CT),idx_cell)
  g_diff.c <- apply(m,1,function(x) {deltas(x,idx_other,idx_cell)}) # rank
  idx_diff <- which(g_diff.c > THRE)
  if (length(idx_diff) == 0){
    print(paste0("failed for ",cell))
    return(c())
  }
  g_diff.c <- g_diff.c[idx_diff]
  mark <- mark[idx_diff]
  names(g_diff.c) <- 1:length(g_diff.c)
  result<-sort(g_diff.c, decreasing = TRUE)
  idx <- as.numeric(names(result))
  #print(length(idx))
  return(mark[idx])
}
by_diff_max_hypo <- function(mark,ref.m,CT,cell,THRE){
  deltas <- function(x,idx_other,idx_cell){
    delta <- min(x[idx_other])-max(x[idx_cell])
    return(delta)
  }
  m <- ref.m[mark,]
  idx_cell <- which(CT == cell)
  idx_other <- setdiff(1:length(CT),idx_cell)
  g_diff.c <- apply(m,1,function(x) {deltas(x,idx_other,idx_cell)}) # rank
  idx_diff <- which(g_diff.c > THRE)
  if (length(idx_diff) == 0){
    print(paste0("failed for ",cell))
    return(c())
  }
  g_diff.c <- g_diff.c[idx_diff]
  mark <- mark[idx_diff]
  names(g_diff.c) <- 1:length(g_diff.c)
  result<-sort(g_diff.c, decreasing = TRUE)
  idx <- as.numeric(names(result))
  #print(length(idx))
  return(mark[idx])
}

by_diff_max_hyper <- function(mark,ref.m,CT,cell,THRE){
  deltas <- function(x,idx_other,idx_cell){
    delta <- min(x[idx_cell])-max(x[idx_other])
    return(delta)
  }
  m <- ref.m[mark,]
  idx_cell <- which(CT == cell)
  idx_other <- setdiff(1:length(CT),idx_cell)
  g_diff.c <- apply(m,1,function(x) {deltas(x,idx_other,idx_cell)}) # rank
  idx_diff <- which(g_diff.c > THRE)
  if (length(idx_diff) == 0){
    print(paste0("failed for ",cell))
    return(c())
  }
  g_diff.c <- g_diff.c[idx_diff]
  mark <- mark[idx_diff]
  names(g_diff.c) <- 1:length(g_diff.c)
  result<-sort(g_diff.c, decreasing = TRUE)
  idx <- as.numeric(names(result))
  #print(length(idx))
  return(mark[idx])
}

dmc_look <- function(m,celltypes,bar = FALSE,name){ #
  hy_c <- list()
  if (!is.factor(celltypes)) {
    stop("Input is not a factor")
  }
  marker.m <- matrix(NA,ncol = 3,nrow = length(levels(celltypes)))
  for (i in seq_along(levels(celltypes))){
    med <- as.character(celltypes)
    other<-setdiff(levels(celltypes), levels(celltypes)[i])
    med[med %in% other] <- "other"
    med <- factor(med)
    design <- model.matrix(~0+med)
    colnames(design) <- levels(med)
    # fit the linear model
    fit <- lmFit(m, design)
    # create a contrast matrix for specific comparisons
    contrast_formula <- paste(levels(celltypes)[i], "- other", sep="")
    contMatrix <- makeContrasts(contrasts = contrast_formula, levels = design)
    # fit the contrasts
    fit2 <- contrasts.fit(fit, contMatrix)
    fit2 <- eBayes(fit2)
    # look at the numbers of DM CpGs at FDR < 0.05
    cout<-summary(decideTests(fit2))
    marker.m[i,] <- cout[,1]
    # filter hyM
    hy_c[[levels(celltypes)[i]]] <- hy(fit2)
  }
  
  if (bar){
    print("Drawing...")
    rownames(marker.m) <- levels(celltypes)
    colnames(marker.m) <- c("Down","NotSig","Up")
    marker.m <- marker.m[,-2]
    marker.m <- as.data.frame(marker.m)
    marker.m$celltypes <- rownames(marker.m)
    marker.m <- reshape2::melt(marker.m,id.vars = "celltypes")
    pdf(paste0(name,"_limma.pdf"),width = 6,height = 6)
    p <- ggplot(marker.m,aes(x = celltypes,y = value,fill = variable)) + 
      geom_bar(stat = "identity",position = "dodge") + 
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      labs(x = "",y = "Number of DMCs",title = "Limma: FDR < 0.05") + 
      scale_fill_manual(values = c("Down" = "blue","Up" = "red"))+
      geom_text(aes(label = value),position = position_dodge(width = 0.9),vjust = -0.5)+
      guides(fill = guide_legend(title = NULL))
    print(p)
    dev.off()
  }
  return(hy_c)
}

dmc_count <- function(dmc_hypo,dmc_hyper){
  dmc.m <- matrix(NA,ncol = 2,nrow = length(dmc_hypo))
  colnames(dmc.m) <- c("Hypo","Hyper")
  rownames(dmc.m) <- names(dmc_hypo)
  dmc.m[,1] <- dmc_hypo
  dmc.m[,2] <- dmc_hyper
  
  dmc.m <- as.data.frame(dmc.m)
  dmc.m$celltypes <- rownames(dmc.m)
  dmc.m <- reshape2::melt(dmc.m,id.vars = "celltypes")
  # distribution of hypo and hyper
  pdf("delta_result.pdf",width = 6,height = 6)
  p <- ggplot(dmc.m,aes(x = celltypes,y = value,fill = variable)) + 
    geom_bar(stat = "identity",position = "dodge") + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(x = "",y = "Number of DMCs",title = "Delta > 0") + 
    scale_fill_manual(values = c("Hypo" = "blue","Hyper" = "red"))+
    geom_text(aes(label = value),position = position_dodge(width = 0.9),vjust = -0.5)+
    guides(fill = guide_legend(title = NULL))
  print(p)
  dev.off()
}
# build reference matrix --------------------------------------------------
gener_ref_med <- function(DMC,m,celltypes){
  mars <- as.character(unlist(DMC))
  ref <- m[mars,]
  ref.mx <- matrix(data = NA,nrow = length(mars),ncol = length(levels(celltypes)))
  rownames(ref.mx) <- mars
  colnames(ref.mx) <- levels(celltypes)
  for ( i in levels(celltypes)){
    idx <- which(celltypes == i)
    if (length(idx) > 1){
      med.c<-apply(ref[,idx],1,median)} # or mean
    else {
      med.c<-ref[,idx]
    }
    ref.mx[,i] <- med.c
  }
  weight <- rep(NA,times = nrow(ref.mx))
  ref <- cbind(ref.mx,weight)
  return(ref)
}

# end to end --------------------------------------------------------------
find_markers <- function(m,celltypes,THRE,class="hypo",name,draw = 1){
  hy <- dmc_look(m,celltypes,bar = T,name)
  DMC <- list()
  if (class == "hypo"){
  for (l in 1:length(hy)){
    na <- names(hy)[l]
    DMC[[na]] <- by_diff_max_hypo(hy[[l]],m,celltypes,na,THRE = THRE)
    print(length(DMC[[na]]))
  }
  }else{
    for (l in 1:length(hy)){
      na <- names(hy)[l]
      DMC[[na]] <- by_diff_max_hyper(hy[[l]],m,celltypes,na,THRE = THRE)
      print(length(DMC[[na]]))
    }  
  }
  choice <- readline(prompt = ("Continue?"))
  if(choice!=1){
    stop("Input new delta threshold!")
  }
if (draw == 1){
DMC <- lapply(DMC,function(x) x[1:min(sapply(DMC,length))])
  library(ComplexHeatmap)
  m_heat <- m[as.character(unlist(DMC)),]
  pdf(paste0(name,"_hm_",class,".pdf"),width = 6,height = 6)
  rainbow_colors <- rainbow(length(levels(celltypes)))
  color_dict <- setNames(rainbow_colors, levels(celltypes))
  col <-  list(celltypes = color_dict)
  split_vector <- unlist(lapply(DMC, function(group) rep(names(DMC)[which(lapply(DMC, function(x) all(x %in% group)) == TRUE)], length(group))))
  p <- Heatmap(m_heat,
          name = "beta",
          column_title = paste0(class,"methylated CpGs of ",name," delta >", THRE),
          use_raster = TRUE,
          show_row_names = FALSE,
          show_column_names = FALSE,
          top_annotation = HeatmapAnnotation(celltypes = celltypes, col = col),
          row_split = split_vector
  )
  print(p)
  dev.off()}
  return(DMC)
}


# boxplot for markers -----------------------------------------------------
markers_overview <- function(m,celltypes,name){
  df <- as.data.frame(t(m))
  df$celltypes <- celltypes
  df_long <- tidyr::pivot_longer(df, -celltypes, names_to = "CpG", values_to = "Value")
  library(ggplot2)
  pdf(paste0(name,"_vl.pdf"),width = 6,height = 6)
  p <- ggplot(df_long, aes(x = celltypes, y = Value, fill = celltypes)) +
    geom_boxplot() +
    #geom_jitter(shape=16, position=position_jitter(0.2)) +
    labs(x = "Cell Type", y = "Beta") +
    theme_minimal()
  print(p)
  dev.off()
  }


# EpiSCORE pro ------------------------------------------------------------
## MSS control -------------------------------------------------------------
mss_filter <- function(markers, mat, cell_type,relax =0) {
  # 将markers和cell_type转换为data.table以加速操作
  markers_dt <- as.data.table(markers)
  cell_type_dt <- as.data.table(cell_type, keep.rownames = "sample")
  
  # 提取唯一的细胞类型
  unique_cell_types <- unique(cell_type)
  k <- length(unique_cell_types) - 1
  # 定义一个内部函数来检查每个基因
  check_gene <- function(gene, celltype) {
    current_cells_idx <- which(cell_type == celltype)
    
    # 计算给定细胞类型中该基因表达值的中位数
    median_in_celltype <- median(mat[gene, current_cells_idx])
    
    if (median_in_celltype != 0) {
      medians_in_others <- sapply(setdiff(unique_cell_types, celltype), function(ctype) {
        other_cells_idx <- which(cell_type == ctype)
        median(mat[gene, other_cells_idx])
      })
      
      # MSS
      MSS <- sum(medians_in_others == 0)
      if (MSS >= k-relax) {
        return(list(gene = gene, cluster = celltype))
      }
    }
    
    return(NULL)
  }
  
  # 使用mclapply并行处理每个基因
  results_list <- mclapply(seq_len(nrow(markers_dt)), function(i) {
    gene <- markers_dt$gene[i]
    celltype <- markers_dt$cluster[i]
    
    if (gene %in% rownames(mat)) { 
      check_gene(gene, celltype)
    } else {
      return(NULL)
    }
    
  }, mc.cores = 20)
  
  
  results_list <- Filter(Negate(is.null), results_list)
  
  if (length(results_list) > 0) {
    filtered_markers_dt <- rbindlist(results_list)
    return(filtered_markers_dt)
  } else {
    return(data.table(gene = character(), cluster = factor(levels(cell_type))))
  }
}

mss_overlap_imp <- function(expref.o,impg.df){
  markes_mss <- rownames(expref.o$ref$med)
  overlap <- intersect(markes_mss,rownames(impg.df))
  #1
  expref.o$ref$med <- expref.o$ref$med[overlap,]
  expref.o$ref$av <- expref.o$ref$av[overlap,]
  #2
  for(i in 1:length(expref.o$markers)){
    overlap_i <- intersect(overlap,rownames(expref.o$markers[[i]]))
    expref.o$markers[[i]] <- expref.o$markers[[i]][overlap_i,]
  }
  return(expref.o)
}


## Impute DNAm reference ---------------------------------------------------
ImputeDNAmRef_mouse <- function(expref.o,pEgX.m,beta.m,tpm.m){
  markers_k <- list()
  for (i in 1:length(expref.o$markers)){
    cell <- names(expref.o$markers)[[i]]
    markers_k[[cell]] <- rownames(expref.o$markers[[i]])
  }
  markers_k <- data.frame(gene = unlist(markers_k), cluster = rep(names(markers_k), times = sapply(markers_k, length)))
  idx <- match(markers_k$gene,rownames(tpm.m))
  tpm_marker <- tpm.m[idx,]
  p_marker <- pEgX.m[rownames(tpm_marker),]
  mf_marker <- beta.m[rownames(tpm_marker),]
  lung_ref <- matrix(data = NA,nrow = nrow(markers_k),ncol = length(unique(markers_k$cluster)))
  colnames(lung_ref) <- unique(markers_k$cluster)
  for (i in 1:nrow(mf_marker)){
    if (all(p_marker[i,] < 0.2)){
      next
    }
    idx <- which(p_marker[i,] < 0.2)
    beta<-median(mf_marker[idx,])
    lung_ref[i,] <- rep(beta,ncol(lung_ref))
    cell <- as.character(markers_k[i,2])
    lung_ref[i,cell] <- 0
  }
  rownames(lung_ref) <- rownames(mf_marker)
  idx <- which(apply(lung_ref,1,function(x) any(is.na(x))))
  lung_ref <- lung_ref[-idx,]
  return(lung_ref)
}


