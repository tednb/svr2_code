GPL <- fread("~/data/infinium/GPL_files/GPL28271/GPL28271-57075.txt")
load("~/data/infinium/mammal_array/GSE199979_liver/mmal_liver_qc.Rd")
# trans coordinate
idx<- match(rownames(qc.o@m),GPL$ID)
mm10_odin <- GPL$Mouse.GRCm38.100_CGstart[idx]
#hep_dmc <- na.omit(mm10_odin[idx_cdc]) # 2874

idx <- which(!is.na(mm10_odin))
m_all <- qc.o@m[idx,]
rownames(m_all) <- mm10_odin[idx]
load("~/Renv/clock_mmal/mmal_hepclock.Rd")
load("clock_mmal/mmal_rawclock.Rd")
rawclock_cpg <- clock_raw$All$beta@Dimnames[[1]][clock_raw$All$beta@i]
hepclock_cpg <- clock_hep$Hep$beta@Dimnames[[1]][clock_hep$Hep$beta@i]

idx<- match(rawclock_cpg,GPL$Mouse.GRCm38.100_CGstart)
rawclock_cpg <- GPL$ID[idx]


clock_cpgs <- rawclock_cpg
databulk <- qc.o@m
allcpg <- rownames(databulk)
databulk <- as.data.frame(apply(databulk, 2, as.numeric))
sapply(databulk, class)
rownames(databulk) <- allcpg



clock.vector <- databulk[clock_cpgs,]


# 1. create matrix with euclidean distances
print("Start creating Matrix of euclidean distances.")


ClockDist <- list()

j=0
for (i in 1:dim(clock.vector)[1]){

  if (length(ClockDist) == 0) {
    meth2 <- sweep(databulk, 2, as.matrix(clock.vector[i,])) # default: "-"
    ClockDist <- sqrt(rowSums(meth2^2))
    j = j+1
  } else {
    meth2 <- sweep(databulk, 2, as.matrix(clock.vector[i,])) # default: "-"
    a <- sqrt(rowSums(meth2^2))
    ClockDist <- cbind(ClockDist, a)
    j = j+1
  }
  if (j == 10) {
    print(paste("Calculated euclidean distances for ", i, "/", length(clock_cpgs), " clock CpGs.", sep = ""))
    j = 0
  } else { next }
}
rownames(ClockDist) <- rownames(databulk)
colnames(ClockDist) <- rownames(clock.vector)

# 2. create a matrix with ranked indices according to the euclidean distances

print("Order CpGs according to euclidean distances.")

sortedindi.Clock <- list()
d=0
for (i in 1:dim(ClockDist)[2])
{
  if (length(sortedindi.Clock)==0)
  {
    sortedindi.Clock <- order(as.numeric(ClockDist[,i]))
    d=d+1
  }
  else
  {
    sortedindi.Clock <- cbind(sortedindi.Clock, order(as.numeric(ClockDist[,i])))
    d=d+1
  }
  if (d == 20)
  {
    print(paste("Sorted euclidean distances for ", i, "/", length(clock_cpgs), " clock CpGs.", sep = ""))
    d=0
  }
}
dim(sortedindi.Clock)

colnames(sortedindi.Clock) <- colnames(ClockDist)


idx<- match(colnames(sortedindi.Clock),GPL$ID)
colnames(sortedindi.Clock) <- GPL$Mouse.GRCm38.100_CGstart[idx]
head(sortedindi.Clock)
# 3. prepare singlecell data

load("~/data/scage/single_cell/hepatocyte/hep_26.Rd")
databulk <- m_all
clock_cpgs <- rawclock_cpg

cpg_row <- lapply(hep.lv,function(cell) {
  share <- intersect(names(cell),rownames(m_all))
  return(share)
})
cpgs <- Reduce(union,cpg_row)

cpg_val <- lapply(hep.lv,function(cell) {
  idx <- match(cpgs,names(cell))
  return(cell[idx])
})

samples <- Reduce(cbind,cpg_val)
colnames(samples) <- names(cpg_val)

samples <- samples[match(rownames(databulk), rownames(samples)),]
head(samples)
rownames(samples) <- rownames(databulk)
################### start
indi.clock <- match(clock_cpgs, rownames(databulk)) # rownumbers of Clock CpGs

if (anyNA(samples[indi.clock,]) == FALSE){
  print("No missing values for clock CpGs in your sample. No imputation needed.")
}
if (anyNA(samples[indi.clock,]) == TRUE){
  vectorclock <- samples[clock_cpgs,]
  a <- sum(is.na(vectorclock[,1]))
  b <- which(is.na(vectorclock))
  c <- sum(is.na(vectorclock))
  b <- b[1:a]
  clockcpgwithNA <- clock_cpgs[b]
  print(paste("Found ",c, " missing clock CpGs in your samples. Imputation starts. This can take some time, depending on the amount of missing values.", sep = ""))
  d = 1
  for (i in 1:dim(samples)[2]){ # loop over experiments
    for (j in 1:length(clock_cpgs)){ #loop over Clock CpGs
      if (is.na(samples[indi.clock[j],i]) == TRUE){ # missing value for Clock CpG
        a = rownames(samples[indi.clock[j],]) # name of Clock CpG of interest
        sortedclockcpg <- samples[,i][sortedindi.Clock[,a]] # get vector with single cell experiment, sorted according to ordered indices of distances
        b = which(!is.na(sortedclockcpg)) # vector with indices of not NA values
        samples[indi.clock[j],i] <- sortedclockcpg[b[1]] # replace missing value with first not NA value in sorted vector
      } else { next}
    }
    print(paste("imputed exp ", d, "/", dim(samples)[2], sep = ""))
    d = d+1
  }
  print("Imputation finished.")
}
impClock <- samples[clock_cpgs,]
