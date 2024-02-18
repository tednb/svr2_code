library(data.table)

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)
file1_path <- args[1]
file2_path <- args[2]

# 读取文件
file1 <- fread(file1_path, header = FALSE)
file2 <- fread(file2_path, header = FALSE)

# 设置列名
colnames(file1) <- c("chr", "pos", "ref", "fraction1", "coverage1", "coverage2")
colnames(file2) <- c("chr", "pos", "ref", "fraction2", "coverage3", "coverage4")

# 合并两个文件
merged_data <- merge(file1, file2, by = c("chr", "pos", "ref"), all = TRUE)

# 计算fraction平均值
merged_data[, avg_fraction := (fraction1 + fraction2) / 2, by = pos]

# 计算各自coverage总和
merged_data[, coverage_sum1 := coverage1 + coverage2, by = pos]
merged_data[, coverage_sum2 := coverage3 + coverage4, by = pos]

# 计算coverage平均值
merged_data[, avg_coverage := (coverage_sum1 + coverage_sum2) / 2, by = pos]

# 选择结果列
result <- merged_data[, .(pos, avg_fraction, avg_coverage)]

# 保存结果到txt文件
fwrite(result, "WGBS.txt")

print("处理完成！")
