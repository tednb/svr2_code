#### Calculation and ordering of euclidean distance ####
###
### Zoe Sassmannshausen 
### zoe.sassmannshausen@dkfz-heidelberg.de
###
### 2023.02.06
###
### 
### 
##
#
# R --vanilla --slave --args TrainingData ClockCpGs Pathsave < estiMAge_2_CalculateDistances.R 

args           <- commandArgs()
TrainingData   <- args[5]
ClockCpGs      <- args[6]
Pathsave       <- args[7]

Sys.time()
Sys.Date()

print("--------------------------------")
print("Load required packages and data.")
print("--------------------------------")

# load libraries
library(glmnet)

Pathsave <- paste(Pathsave, "/", sep = "")
if (!dir.exists(Pathsave)) {dir.create(Pathsave)}
print(paste("Save results in directory ", Pathsave, sep = ""))

# load required data
print("Load training data.")
databulk      <- read.csv2(TrainingData, sep = ",", row.names = 1) 
databulk      <- databulk[-1,]
print(paste("Number of samples: ", dim(databulk)[2], ", Number of CpG sites: ", dim(databulk)[1], sep = ""))

clock_cpgs    <- read.csv(ClockCpGs, sep = ",") # vector of your clock CpGs
clock_cpgs    <- clock_cpgs[,2]
print("loaded clock CpGs.")
print(paste("Number of your Clock CpGs:", length(clock_cpgs), sep = " "))

######
#teste:geht auch ohne?
allcpg <- rownames(databulk)
databulk <- as.data.frame(apply(databulk, 2, as.numeric))
sapply(databulk, class)
rownames(databulk) <- allcpg

clock.vector <- databulk[clock_cpgs,] 

print("--------------------------------")
print("Calculate euclidean distances.")
print("--------------------------------")

# 1. create matrix with euclidean distances
print("Start creating Matrix of euclidean distances.")


ClockDist <- list()

j=0
for (i in 1:dim(clock.vector)[1])
{
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

save(ClockDist, file = paste(Pathsave, "EUC_Dists.RData", sep = ""))
print("Euclidean distances are saved.")


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
head(sortedindi.Clock)

save(sortedindi.Clock, file = paste(Pathsave, "EUC_sorted.RData", sep = ""))
print("Sorted indices are saved.")

print("--------------------------------")
print("Script finished.")
print("--------------------------------")

Sys.time()
Sys.Date()