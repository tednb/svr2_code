load("~/result/clock/train/GSE180474_clock.Rd")
load("~/data/infinium/GPL_files/GPL21145/annoEPICv1B4.Rd")
library(data.table)
clock_Cpg <- re_epic[[1]]$Hep$beta@Dimnames[[1]][re_epic[[1]]$Hep$beta@i]

GPL <- as.data.frame(annoEPICv1B4.m)
idx <- match(clock_Cpg,GPL$IlmnID)
info <- GPL[idx,]
codin <- as.numeric(info$MAPINFO)