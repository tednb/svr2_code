setwd("~/data/RRBS/GSE120132/liver")
library(pbapply)
files <- list.files()
input <- list()
for (i in files){
  name <- sub("^(GSM\\d+)_.*$", "\\1", i)  #以文件名中的GSM编号作为nameda's'd
  input[[name]] <- read.table(i)[,c(1,2,3,6)]  #使用name给每个文件读取的数据框命名
}
df3 <- pblapply(input, function(x) {data.frame(paste(x$V1,"_","cg",x$V3,sep = ""),x$V6)})
data <- as.data.frame(matrix(nrow = 1992813, ncol = 60))
rownames(data) <- df1$GSM3394354$paste0.x.V2..x.V3.
