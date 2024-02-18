
load("~/chain/MM285_map/Tss200_gene_to_cpgs.Rdata")
MM285_map_tss <- function(beta_matrix, gene_to_cpgs) {
  gene_matrix <- data.frame(matrix(ncol = ncol(beta_matrix), nrow = 0))
  colnames(gene_matrix) <- colnames(beta_matrix)
  
  # 初始化进度条
  pb <- txtProgressBar(min = 0, max = length(names(gene_to_cpgs)), style = 3)
  
  for(i in seq_along(names(gene_to_cpgs))) {
    gene <- names(gene_to_cpgs)[i]
    
    # 更新进度条
    setTxtProgressBar(pb, i)
    
    # 确保只选择存在于beta_matrix中的CpG位点
    cpgs_for_gene <- gene_to_cpgs[[gene]]
    valid_cpgs <- cpgs_for_gene %in% rownames(beta_matrix)
    
    # 如果没有有效的CpG位点，则跳过这个基因
    if(!any(valid_cpgs)) next
    
    # 提取这个基因对应所有有效CpG位点在beta矩阵中的值
    cpg_values <- beta_matrix[cpgs_for_gene[valid_cpgs], , drop = FALSE]
    
    # 检查是否至少有一个CpG位点被选中
    if(is.null(dim(cpg_values))) {
      cpg_values <- t(as.matrix(cpg_values))
    }
    
    avg_values <- apply(cpg_values, 2, mean, na.rm = TRUE)
    
    gene_matrix[gene, ] <- avg_values
  }
  
  # 关闭进度条
  close(pb)
  
  gene_matrix <- gene_matrix[complete.cases(gene_matrix), ]
  
  return(gene_matrix)
}

