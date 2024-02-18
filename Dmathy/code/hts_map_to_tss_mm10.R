library(parallel)
library(dplyr)
library(tibble)
depth_mask <- function(mat, dep.m, thre) {
  if (inherits(mat, "matrix")) {
    if (!inherits(dep.m, "matrix")) {
      stop("Both 'mat' and 'dep.m' must be matrices or both must be lists.")
    }
    mat[dep.m < thre] <- NA
  } else if (inherits(mat, "list")) {
    if (!inherits(dep.m, "list")) {
      stop("Both 'mat' and 'dep.m' must be matrices or both must be lists.")
    }
    len <- length(mat)
    for (n in seq_len(len)) {
      dep <- as.numeric(as.matrix(dep.m[[n]])[,2])
      mat[[n]][dep < thre, 2] <- NA
    }
  } else {
    stop("Input 'mat' must be either a matrix or a list of matrices.")
  }
  return(mat)
}



map_to_tss <- function(m, gene_info, region, cores = 1) {
  is_matrix <- inherits(m, "matrix")
  is_list <- inherits(m, "list")
  
  if (!is_matrix && !is_list) {
    stop("The input 'm' must be a matrix or list of matrices.")
  }
  
  process_locations <- function(locations) {
    split_locs <- strsplit(locations, "-")
    list(chr = vapply(split_locs, function(x) x[1], character(1)),
         codin = vapply(split_locs, function(x) as.numeric(x[2]), numeric(1)))
  }
  
  process_gene <- function(i, chr, codin, data,region) {
    gene_chr <- gene_info[i, 2]
    gene_strand <- gene_info[i, 5]
    tss <- if (gene_strand == "+") gene_info[i, 3] else gene_info[i, 4]
    idx <- if (gene_strand == "+") {
      which(chr == gene_chr & codin <= tss & codin >= tss - region)
    } else {
      which(chr == gene_chr & codin >= tss & codin <= tss + region)
    }
    if (length(idx) > 0) {
      re <- colMeans(data[idx, , drop = FALSE], na.rm = TRUE)
      re[is.nan(re)] <- NA
      return(re)
    } else {
      re <- rep(NA, ncol(data))
      return(re)
    }
  }
  min_na <- function(x) {
    if (all(is.na(x))) {
      return(NA)
    } else {
      return(min(x, na.rm = TRUE))
    }
  }
  mf.m <- matrix(NA, ncol = if (is_matrix) ncol(m) else length(m), nrow = nrow(gene_info))
  colnames(mf.m) <- if (is_matrix) colnames(m) else sapply(m, function(x) names(x)[2])
  if (is_matrix) {
    locs <- process_locations(rownames(m))
    results_list <- pbapply::pblapply(seq_len(nrow(gene_info)), function(i) process_gene(i, locs$chr, locs$codin, m,region), cl = cores)
    for (i in seq_along(results_list)) {
      mf.m[i, ] <- results_list[[i]]
    }
  }
  
  if (is_list) {
    for (j in seq_along(m)) {
      locs <- process_locations(m[[j]]$location)
      mf.m[, j] <- unlist(pbapply::pblapply(seq_len(nrow(gene_info)), function(i) process_gene(i, locs$chr, locs$codin, m[[j]][, 2,drop = F],region), cl = cores))
    }
  }
  # Remove rows with only NAs
  idx_na_rows <- which(rowSums(is.na(mf.m)) == ncol(mf.m))
  mf.m <- as.data.frame(mf.m[-idx_na_rows, , drop = FALSE])
  mf.m$mgi_symbol <- gene_info$mgi_symbol[-idx_na_rows]
  library(dplyr)
  mf.m_cleaned <- mf.m %>%
    group_by(mgi_symbol) %>%
    summarise(across(everything(), \(x) min_na(x))) %>%
    ungroup()
  rown <- mf.m_cleaned$mgi_symbol
  mf.m_cleaned$mgi_symbol <- NULL
  mf.m_cleaned <- as.matrix(mf.m_cleaned)
  rownames(mf.m_cleaned) <- rown
  return(as.matrix(mf.m_cleaned))
}


library(ggplot2)

get_tss <- function(start_position, end_position, strand) {
  if (strand == "+") {
    return(as.numeric(start_position))
  } else {
    return(as.numeric(end_position))
  }
}

get_tts <- function(start_position, end_position, strand) {
  if (strand == "+") {
    return(as.numeric(end_position))
  } else {
    return(as.numeric(start_position))
  }
}

dtr_tss <- function(gene_info, mat, sample_idx, gene_idx) {
  
  # 获取样本名称
  sample_name <- colnames(mat)[sample_idx]
  
  # 获取基因信息
  gene <- gene_info[gene_idx, ]
  
  # 获取 TSS 和 TTS
  tss <- get_tss(gene$start_position, gene$end_position, gene$strand)
  tts <- get_tts(gene$start_position, gene$end_position, gene$strand)
  
  # 计算范围
  start_range <- tss - 500
  end_range <- tts + 500
  
  # 提取甲基化水平数据，确保位置存在于mat中
  position_vec <- seq(from = min(start_range, tts - 500), to = max(end_range, tss + 500))
  valid_positions <- position_vec %in% rownames(mat)
  
  # 如果没有有效位置直接返回错误信息
  if (!any(valid_positions)) {
    stop("No valid methylation data found for the specified positions.")
  }
  
  position_vec <- position_vec[valid_positions]
  meth_data <- mat[as.character(position_vec), sample_idx, drop = FALSE]
  
  # 移除包含NA值的行
  valid_indices <- !is.na(meth_data)
  meth_data <- meth_data[valid_indices, , drop = FALSE]
  position_vec <- position_vec[valid_indices]
  
  # 创建数据框
  methylation_data_frame <- data.frame(Position = as.numeric(position_vec),
                                       Methylation = as.numeric(meth_data),
                                       Gene = rep(gene$mgi_symbol, length(meth_data)))
  
  # 绘图
  p <- ggplot(methylation_data_frame, aes(x = Position, y = Methylation)) +
    geom_line(aes(group = Gene), color="blue") +
    theme_minimal() +
    labs(title = paste("Methylation levels around", gene$mgi_symbol,
                       "for", sample_name),
         x = "Position",
         y = "Methylation") +
    scale_x_continuous(breaks = c(ifelse(tss > tts, tss + 500, tts - 500), tss, tts, ifelse(tss > tts, tts - 500, tts + 500)),
                       labels = c("TSS500", "TSS", "TTS", "TTS500"))
  
  print(p)
}

# 假设gene_info和mat已经被正确定义，并且sample_idx和gene_idx是单个数字。
# 这里不包括这些定义，您需要在调用dtr_tss之前定义它们。

calculate_methylation_tss <- function(gene, data, upstream = 10000, downstream = 10000,region) {
  # 根据strand决定如何计算upstream和downstream
  if (gene$strand == "-") {
    TSS <- gene$end_position
    TTS <- gene$start_position
    region_upstream_end <- TSS + upstream
    region_tss_start <- TSS + region
    region_tss_end <- TSS
    region_tss_back <- TSS - region
    
    region_tts_start <- TTS
    region_downstream_end <- TTS - downstream
    
  } else {
    TSS <- gene$start_position
    TTS <- gene$end_position
    region_upstream_end <- TSS - upstream
    region_tss_start <- TSS - region
    region_tss_end <- TSS
    region_tss_back <- TSS + region
    
    region_tts_start <-  TTS
    region_downstream_end <- TTS + downstream
  }
  
  # 初始化存储每个区域甲基化水平的向量
  methylation_levels <- list()
  
  # 计算每个区域内的CpG位点的甲基化水平平均值
  methylation_levels$upstream_to_TSS_200 <- calculate_region_methylation(region_upstream_end, region_tss_start, gene, data)
  
  methylation_levels$TSS_200_to_TSS <- calculate_region_methylation(region_tss_start, region_tss_end, gene, data)
  
  methylation_levels$TSS_to_TSS_n200 <- calculate_region_methylation(region_tss_end, region_tss_back, gene, data)
  
  methylation_levels$TSS_n200_to_TTS <- calculate_region_methylation(region_tss_back, region_tts_start, gene, data)
  
  methylation_levels$TTS_to_downstream <- calculate_region_methylation(region_tts_start, region_downstream_end, gene, data)
  
  return(unlist(methylation_levels))
}

calculate_region_methylation <- function(start_pos, end_pos, gene, data) {
  # 根据基因所在的链来选择正确的比较运算符
  if (gene$strand == "-") {
    cpgs_in_region <- data[data$chr == as.character(gene$chromosome_name) &
                             data$Start <= start_pos &
                             data$End >= end_pos, ]
  } else {
    cpgs_in_region <- data[data$chr == as.character(gene$chromosome_name) &
                             data$Start >= start_pos &
                             data$End <= end_pos, ]
  }
  if(nrow(cpgs_in_region) != 0){
    mean_cpgs_in_region <- mean(cpgs_in_region$MF, na.rm = TRUE)
  }else {
    mean_cpgs_in_region <- NA
  }
  return(mean_cpgs_in_region)
}
