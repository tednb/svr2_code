setwd("~/data/MouseDNAmAtlas/GSE42836")
files <- list.dirs(recursive = FALSE)
library(GenomicRanges)
library(rtracklayer)
library(data.table)
chain <- import.chain("~/chain/mm9ToMm10.over.chain")
mf <- list()
for (i in 1:length(files)) {
  print(i)
  setwd(files[i])
  file_names <- list.files(pattern = "\\.calls.txt$", full.names = TRUE)
  raw.m <- fread(file_names)
  # Keep only the rows where all needed columns are not NA
  complete_cases <- complete.cases(raw.m$chromosome, as.numeric(raw.m$`CpG location left`), as.numeric(raw.m$`CpG location right`), raw.m$`#CpG sequenced`, raw.m$`%mCG`)
  raw.m <- raw.m[complete_cases, ]
  gr <- GRanges(seqnames = Rle(raw.m$chromosome),
                ranges = IRanges(start = as.numeric(raw.m$`CpG location left`), end = as.numeric(raw.m$`CpG location right`)),
                mcols = DataFrame(cov=raw.m$`#CpG sequenced`, beta=raw.m$`%mCG`))
  
  gr_lifted <- liftOver(gr, chain)
  gr_lifted <- unlist(gr_lifted)
  # Remove NA values after liftOver
  gr_lifted <- gr_lifted[!is.na(gr_lifted)]
  
  odin <- start(gr_lifted)
  cov <- mcols(gr_lifted)[,1]
  beta <- mcols(gr_lifted)[,2]
  
  mf[[i]] <- as.matrix(data.frame(odin, beta, cov))
  setwd("..")
}



shared_elements <- Reduce(intersect, lapply(mf, function(x) x[,1]))
idx <- match(shared_elements,)
sample_name <- c()
for (i in 1:length(files)) {
  print(i)
  setwd(files[i])
  file_names <- list.files(pattern = "\\.calls.txt$", full.names = TRUE)
  pattern <- "GSM(\\d+)"
  match <- regexpr(pattern, file_names)
  sample_name[i]<-regmatches(file_names, match)
  setwd("..")
  }
# 提取每个矩阵的第二列和第三列，并构建两个新的大矩阵
matrix_2nd_col <- do.call(cbind, lapply(mf, function(x) x[match(shared_elements,x[,1]), 2]))
matrix_3rd_col <- do.call(cbind, lapply(mf, function(x) x[match(shared_elements,x[,1]), 3]))




# 设置行名和列名
rownames(matrix_2nd_col) <- as.character(shared_elements)
colnames(matrix_2nd_col) <- sample_name
rownames(matrix_3rd_col) <- as.character(shared_elements)
colnames(matrix_3rd_col) <- sample_name

# 输出第二列大矩阵为txt文件
write.table(matrix_2nd_col, file = "beta.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# 输出第三列大矩阵为txt文件
write.table(matrix_3rd_col, file = "coverage.txt", sep = "\t", quote = FALSE, row.names = TRUE)
