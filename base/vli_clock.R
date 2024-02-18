setwd("/mnt/local-disk/data/guoxiaolong/Renv/")
load("~/data/infinium/liver_sample/341_liver_samples/RData/GSE180474_he.Rd")
source("~/code/Dmathy/ELN_clock.R")
# 1. test
dir.create("/mnt/local-disk/data/guoxiaolong/Renv/train/")
new_folder_path <-"/mnt/local-disk/data/guoxiaolong/Renv/train/"
he.o@name <- "GSE180474"
sets <- extract(he.o,seed = 1234)
trainset <- sets[[1]]
testset <- sets[[2]]
covs <- model.matrix(~ sex + t2d, data = trainset@s)
seqs <- seq(0.001,0.7,0.001)
clocks <- clock(trainset,covs,seqs)
testset@GSE <- "GSE180474"
age_pre <- validate(testset,clocks)
file.copy(from = list.files(path = "/mnt/local-disk/data/guoxiaolong/Renv/", pattern = ".pdf"), to = new_folder_path, overwrite = TRUE)
file.remove(list.files("/mnt/local-disk/data/guoxiaolong/Renv", pattern = ".pdf"))

# 2. GSE123995
dir.create("/mnt/local-disk/data/guoxiaolong/Renv/HEP/")
new_folder_path <-"/mnt/local-disk/data/guoxiaolong/Renv/HEP/"
load("~/data/infinium/liver_sample/AA_hepatocytes/GSE123995_qc.Rd")
age_pre <- pip_pre(qc.o,trainset,covs,seqs)
file.copy(from = list.files(path = "/mnt/local-disk/data/guoxiaolong/Renv/", pattern = ".pdf"), to = new_folder_path, overwrite = TRUE)
file.remove(list.files("/mnt/local-disk/data/guoxiaolong/Renv", pattern = ".pdf"))

# 3. GSE103078
dir.create("/mnt/local-disk/data/guoxiaolong/Renv/40/")
new_folder_path <-"/mnt/local-disk/data/guoxiaolong/Renv/40/"
load("~/data/infinium/liver_sample/40_normal_liver_samples/GSE107038_qc.Rd")
pip_pre(qc_107038.o,trainset,covs,seqs)
file.copy(from = list.files(path = "/mnt/local-disk/data/guoxiaolong/Renv/", pattern = ".pdf"), to = new_folder_path, overwrite = TRUE)
file.remove(list.files("/mnt/local-disk/data/guoxiaolong/Renv", pattern = ".pdf"))

# 4. GSE61446
dir.create("/mnt/local-disk/data/guoxiaolong/Renv/obese/")
new_folder_path <-"/mnt/local-disk/data/guoxiaolong/Renv/obese/"
load("~/data/infinium/liver_sample/67_obese_livers/GSE61446/GSE61446_qc.Rd")
pip_pre(qc.o,trainset,covs,seqs)
file.copy(from = list.files(path = "/mnt/local-disk/data/guoxiaolong/Renv/", pattern = ".pdf"), to = new_folder_path, overwrite = TRUE)
file.remove(list.files("/mnt/local-disk/data/guoxiaolong/Renv", pattern = ".pdf"))

# 5. GSE61258
load("~/data/infinium/liver_sample/79_livers(25_healthy)/GSE61258/GSE61258_qc.Rd")
dir.create("/mnt/local-disk/data/guoxiaolong/Renv/79/")
new_folder_path <-"/mnt/local-disk/data/guoxiaolong/Renv/79/"
age_pre <- pip_pre(qc.o,trainset,covs,seqs)
file.copy(from = list.files(path = "/mnt/local-disk/data/guoxiaolong/Renv/", pattern = ".pdf"), to = new_folder_path, overwrite = TRUE)
file.remove(list.files("/mnt/local-disk/data/guoxiaolong/Renv", pattern = ".pdf"))