library(PCBS)
library(ggplot2)

setwd("C:\\PROJECTS\\Shane\\Harding_combined\\local\\PCBS_res")

eigen_mutation <- read.delim("coverage_matrix_mutation.tab", header = T, sep = "\t")
dim(eigen_mutation)
eigen_mutation <- eigen_mutation[, -c(2:5)]
eigen_radiation <- read.delim("coverage_matrix_radiation.tab", header = T, sep = "\t")
dim(eigen_radiation)
eigen_radiation <- eigen_radiation[, -c(2:5)]

pca_mutation <- DefineBestPC(eigen_mutation, IDs = c("treat", "control")) 
pca_radiation <- DefineBestPC(eigen_radiation, IDs = c("treat", "control")) # IDs segregate two conditions based
# on a common identifier in the column names of the eigen dataframe, rather than 
# by column number. You may have to rename the columns of your input object 
# if no common name exists.

getPCRanks(eigen_mutation, IDs = c("treat", "control"), PC = 1) -> ranks_mutation # Get an eigenvector score for each CpG based on principal component 1
getPCRanks(eigen_radiation, IDs = c("treat", "control"), PC = 2) -> ranks_radiation

rd_mutation <- rankDist(ranks_mutation, mode="intersect") # Two modes "intersect" and "strict"
rd_radiation <- rankDist(ranks_radiation, mode="intersect") 

ranks <- ranks_radiation

DMLs <- addRanks(ranks) # add rank order to our CpGs
DMLs$significant <- DMLs$abs.order <= 511 # significant CpG cut-off defined by rankDist() is: 
head(DMLs)

dim(DMLs[DMLs$significant,])

test_5000 <- checkRank(ranks, 5000) # set cut-off to 5000
print(test_5000)
test_10000 <- checkRank(ranks, 35000) # set cut-off to 35000
print(test_10000)
test_200k <- checkRank(ranks, 200000) # set cut-off to 200k
print(test_200k)

chromDictObj <- chromDict(ranks) # Optional step to create a chromDict object, 
# this step is computationally intensive and run internally in many functions, 
# so it is recommended to generate it one time only.


c() -> sig
seeds <- c(100000, 150000, 200000, 250000, 300000)
for(seed in seeds){
  DMR <- Get_Novel_DMRs(ranks, seed, chromDictObj = chromDictObj, QueryLimit=1000, DMR_resolution=100, minCpGs=10)
  c(sig, nrow(DMR[DMR$FDR <= 0.05,])) -> sig
}

p <- ggplot(data.frame(nSeeds=seeds,
                  nDMRs=sig), aes(x=nSeeds, y=nDMRs))+
  geom_line()+ geom_point(size=2)+ Ol_Reliable()

print(p)


CheckOvercompression(ranks, 200000) # For larger genomes, the target seed number will usually be about 1-2% of the number of rows in the ranks file. 
# But in this small test set, we are just testing values around the DML cut-off.

DMRs <- Get_Novel_DMRs(ranks, 200000, chromDictObj=chromDictObj, QueryLimit=1000, DMR_resolution=100, minCpGs=10) # No overcompression was detected, 
# so we will just use the largest seed value.

head(DMRs[order(DMRs$FDR, decreasing = F),])
dim(DMRs[DMRs$FDR<0.1,])
