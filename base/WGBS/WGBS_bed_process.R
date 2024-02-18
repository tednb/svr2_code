setwd("~/data/scage/REF/raw_data/")
library(data.table)
library(dplyr)

files <- list.dirs(recursive = FALSE)
chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                 "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chrY","chrX")

process_bed_file <- function(file_path) {
  # 读取文件并创建data.table
  data <- fread(file_path)
  setDT(data)
  # 创建染色体名(chr)和坐标(coordinate)列，并转换cov/MF值为数值类型
  # data[, c("chr", "coordinate", "cov", "MF") := .(
  #   chr = as.character(V1),                      # 使用V1作为染色体名
  #   coordinate = as.numeric(gsub("[^0-9]", "", V2)), # 假设V2包含非数字字符需要清理
  #   cov = as.numeric(V10),
  #   MF = as.numeric(gsub("[^0-9.]", "", V11))/100   # 假设V11包含非数字字符需要清理，且保留小数点
  # )]
  data[, c("chr", "coordinate", "cov", "MF") := .(
    chr = V1,   # 使用V1作为染色体名 paste0("chr",as.character(V1))
    coordinate = as.numeric(gsub("[^0-9]", "", V3)), # 假设V2包含非数字字符需要清理
    cov = as.numeric(V5),
    MF = as.numeric(V6)  # 假设V11包含非数字字符需要清理，且保留小数点
  )]
  
  # 删除chr *** random,chrUn 和chrM  和 V1-V14列
  # data[,1:14] <- NULL
  data[,1:7] <- NULL
  idx <- which(data$chr %in% chromosomes)
  data <- data[idx,]
  # 合并chr和坐标
  data[, location := paste(chr,coordinate,sep = "-")][,c("chr","coordinate") := NULL]
  return(data)
}


all_final_results_MF <- list()
all_final_results_cov <- list()

for (i in 1:length(files)) {
  print(i)
  setwd(files[i])
  file_names <- list.files(pattern = "\\.txt$", full.names = TRUE)
  all_results <- list()
  # progerss bar
  total <- length(file_names)
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  for (j in 1:length(file_names)) {
    # 直接传递文件路径给 process_bed_file 函数
    processed_data <- process_bed_file(file_names[j])
    # 选择并重命名列
    result_MF <- processed_data[, .(location, MF)]
    result_cov <- processed_data[, .(location, cov)]
    setnames(result_MF, "MF", paste0(gsub("\\./","",files[i]),"_rep",j))
    setnames(result_cov, "cov", paste0(gsub("\\./","",files[i]),"_rep",j))
    # 将结果添加到列表
    all_final_results_MF[[length(all_final_results_MF) + 1]] <- result_MF
    all_final_results_cov[[length(all_final_results_cov) + 1]] <- result_cov
    gc()
    setTxtProgressBar(pb, j)
  }
  close(pb)
  setwd("..")
}

# merged_result_MF <- Reduce(function(x, y) {merge(x, y, by = "location", all = FALSE)}, all_final_results_MF)
# merged_result_cov <- Reduce(function(x, y) {merge(x, y, by = "location", all = FALSE)}, all_final_results_cov)
# 
# write.table(merged_result_MF, file = "~/data/scage/REF/raw_data/merged_liver_MF.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(merged_result_cov, file = "~/data/scage/REF/raw_data/merged_liver_cov.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save(all_final_results_MF,all_final_results_cov,file = "~/data/scage/REF/raw_data/liver30_list.Rdata")
