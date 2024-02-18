
# 0 -----------------------------------------------------------------------

library(data.table)
library(dplyr)
gencode <- fread("~/data/infinium/GPL_files/GPL30650/MM285.mm10.manifest.gencode.vM25.tsv.gz",sep = "\t")

dist <- lapply(strsplit(gencode$distToTSS,";"),as.numeric)
strand <- gencode$probe_strand
ID <- gencode$probeID
genesymbol <- strsplit(gencode$geneNames,";")
# 
# gtf_vM25 <- rtracklayer::import('~/chain/gencode.vM25.annotation.gtf.gz') %>% as.data.frame()
# gtf_mrna_vM25 <- dplyr::select(gtf_vM25,c("gene_type", "gene_name"))%>%
#   subset(., gene_type == "protein_coding") %>% unique()

# 1 -----------------------------------------------------------------------
pb <- txtProgressBar(min = 0, max = length(dist), style = 3)
group <- list()
for (i in 1:length(dist)) {
  if(strand[i] == "+") {
    group[[i]] <- sapply(dist[[i]], function(x) 
      ifelse(x < 0 & x >= -200, 1,
             ifelse(x < 0 & x >= -500, 2,
                    ifelse(x < 0 & x >= -1500, 3, 0))))
    names(group[[i]]) <- genesymbol[[i]]
  } else {
    group[[i]] <- sapply(dist[[i]], function(x) 
      ifelse(x > 0 & x <= 200, 1,
             ifelse(x > 0 & x <= 500, 2,
                    ifelse(x > 0 & x <= 1500, 3, 0))))
    names(group[[i]]) <- genesymbol[[i]]
  }
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 2 -----------------------------------------------------------------------
# 定义一个函数来去除列表中向量的重复项，并且排除任何包含NA的名字或值
unique_list_elements <- function(lst) {
  # 遍历列表中的每个元素
  lapply(lst, function(vec) {
    # 排除名字或值为NA的元素
    vec <- vec[!is.na(names(vec)) & !is.na(vec)]
    # 将向量转换为数据框（保留名字和值）
    df <- data.frame(name = names(vec), value = vec, stringsAsFactors = FALSE)
    # 去除数据框中的重复行
    unique_df <- unique(df)
    # 将唯一的数据框转换回向量，并设置名字
    setNames(unique_df$value, unique_df$name)
  })
}

group_unique <- unique_list_elements(group)

remove_empty_or_all_zero <- function(lst) {
  indices_to_remove <- sapply(lst, function(vec) {
    is_empty <- is.logical(vec) && length(vec) == 0
    all_zeros <- all(vec == 0)
    return(is_empty | all_zeros)
  })
  indices_to_remove
}
idx <- remove_empty_or_all_zero(group_unique)
group_filter <- group_unique[!idx]
ID_filter <- ID[!idx]
names(group_filter) <- ID_filter
save(group_filter,file = "mm285_mapinfo.Rdata")

# 3 -----------------------------------------------------------------------
# 初始化一个空列表来收集每个基因对应的CpG位点
gene_to_cpgs <- list()
# 遍历group_filter列表
for(cpg_id in names(group_filter)) {
  # 获取当前CpG位点映射到的所有基因及其值
  genes <- names(group_filter[[cpg_id]])
  values <- group_filter[[cpg_id]]
  
  # 对于每个基因，检查其值是否为1或2，并相应地更新gene_to_cpgs列表
  for(i in seq_along(genes)) {
    gene <- genes[i]
    value <- values[i]
    # 只关注值为1或2的基因
    #if(value == 1 || value == 2) { #TSS500
    if(value == 1) { #TSS200
      # 如果当前基因已经在gene_to_cpgs中，则添加这个CpG位点到基因对应的向量中
      if(gene %in% names(gene_to_cpgs)) {
        gene_to_cpgs[[gene]] <- c(gene_to_cpgs[[gene]], cpg_id)
      } else { # 如果当前基因不在gene_to_cpgs中，则创建一个新向量并添加这个CpG位点
        gene_to_cpgs[[gene]] <- cpg_id
      }
    }
  }
}
gtf_vM25 <- rtracklayer::import('~/chain/gencode.vM25.annotation.gtf.gz') %>% as.data.frame()
gtf_mrna_vM25 <- dplyr::select(gtf_vM25,c("gene_type", "gene_name"))%>%
  subset(., gene_type == "protein_coding") %>% unique()
shared <- intersect(gtf_mrna_vM25$gene_name,names(gene_to_cpgs))
gene_to_cpgs <- gene_to_cpgs[shared]
#save(gene_to_cpgs, file = "Tss500_gene_to_cpgs.Rdata")
save(gene_to_cpgs, file = "Tss200_gene_to_cpgs.Rdata")


# 4 -----------------------------------------------------------------------

# 初始化一个空的数据框来存储最终的基因beta值矩阵
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
  
  gene_matrix <- as.matrix(gene_matrix[complete.cases(gene_matrix), ])
  
  return(gene_matrix)
}




