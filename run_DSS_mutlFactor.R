library(DSS)
library(bsseq)
library(comethyl)
library(dplyr)

bsdir <- "C:\\PROJECTS\\Shane\\Harding_combined\\local\\DMRcate_res"

bs <- readRDS(file.path(bsdir, "bs.RDS"))

colnames(bs) <- sapply((strsplit(colnames(bs), split=".", fixed = T)), function(x)x[1])
standardChr <- c(paste0("chr", seq.int(22)), "chrX", "chrY")
bsf <- bs[seqnames(bs) %in% standardChr, ]
seqlevels(bsf) <- standardChr

setwd("C:\\PROJECTS\\Shane\\Harding_combined\\local\\DSS_res")

design <- as.data.frame(Reduce(rbind, strsplit(colnames(bsf), split="-")))
colnames(design) <- c("pheno", "treat", "rep")
rownames(design) <- colnames(bsf)
design$batch <- rep(c(1,2,2), by=4)

design <- design |>
  mutate(treat=factor(treat, levels=c("NIR", "IR"))) |>
  mutate(pheno=factor(pheno, levels= c("WT", "R172K"))) |>
  mutate(batch=factor(batch, levels=c(1,2))) |>
  mutate(group=factor(paste(pheno, treat, sep="-"))) |>
  mutate(group=relevel(group, "WT-NIR"))


# three factor model ####
DMLfit= DMLfit.multiFactor(bsf, design = design, formula = ~batch+treat+pheno+treat:pheno)
   
DMLfit$X                 

DMLtestTreat <- DMLtest.multiFactor(DMLfit, coef = 2)
DMLtestPheno <- DMLtest.multiFactor(DMLfit, coef = 3)
DMLtestBatch <- DMLtest.multiFactor(DMLfit, coef = 4)
DMLtestInter <- DMLtest.multiFactor(DMLfit, coef = 5)

DMLtestTreat[DMLtestTreat$fdrs<0.2,]
DMLtestPheno[DMLtestPheno$fdrs<0.2,]
DMLtestBatch[DMLtestBatch$fdrs<0.05,]
DMLtestInter[DMLtestInter$fdrs<0.5,]

dmrsTreat = callDMR(DMLtestTreat, p.threshold = 0.05)
dmrsPheno = callDMR(DMLtestPheno, p.threshold = 0.05)
dmrsBatch = callDMR(DMLtestBatch, p.threshold = 0.05)
dmrsInter = callDMR(DMLtestInter, p.threshold = 0.05)

dmrs_list <- list(dmrsTreat, dmrsPheno, dmrsBatch, dmrsInter)
names(dmrs_list) <- c("treat", "pheno", "batch", "interaction")

for(aname in names(dmrs_list)){
  dmrs <- dmrs_list[[aname]]
  dmrs$name <- rownames(dmrs)
  dmrs$strand <- rep("*", nrow(dmrs))
  dmrs <- dmrs[, c("chr", "start", "end", "name", "areaStat", "strand", "length", "nCG")]
  
  write.table(dmrs, paste0("DMR_", aname, ".bed"), row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
}

write.table(dmrsTreat, "DMR_treat.tsv", row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
write.table(dmrsPheno, "DMR_pheno.tsv", row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
write.table(dmrsBatch, "DMR_batch.tsv", row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
write.table(dmrsInter, "DMR_interaction.tsv", row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)

showOneDMR(dmrsTreat[1,], bsf)

# two factor model ####
DMLfit= DMLfit.multiFactor(bsf, design = design, formula = ~batch+group)

DMLfit$X                 

DMLtestTreatWT <- DMLtest.multiFactor(DMLfit, coef = 5)
DMLtestPhenoNIR <- DMLtest.multiFactor(DMLfit, coef = 4)
DMLtestBatch2 <- DMLtest.multiFactor(DMLfit, coef = 2)
DMLtestInter2 <- DMLtest.multiFactor(DMLfit, coef = 3)

DMLtestTreatWT <- na.omit(DMLtestTreatWT)
DMLtestPhenoNIR <- na.omit(DMLtestPhenoNIR)
DMLtestBatch2 <- na.omit(DMLtestBatch2)
DMLtestInter2 <- na.omit(DMLtestInter2)

dmrsTreatWT = callDMR(DMLtestTreatWT, p.threshold = 0.001)
dmrsPhenoNIR = callDMR(DMLtestPhenoNIR, p.threshold = 0.001)
dmrsBatch2 = callDMR(DMLtestBatch2, p.threshold = 0.001)
dmrsInter2 = callDMR(DMLtestInter2, p.threshold = 0.001)

dmrs_list <- list(dmrsTreatWT, dmrsPhenoNIR, dmrsBatch2, dmrsInter2)
names(dmrs_list) <- c("treatWT", "phenoNIR", "batch2", "interaction2")

for(aname in names(dmrs_list)){
  dmrs <- dmrs_list[[aname]]
  dmrs$name <- rownames(dmrs)
  dmrs$strand <- rep("*", nrow(dmrs))
  dmrs <- dmrs[, c("chr", "start", "end", "name", "areaStat", "strand", "length", "nCG")]
  
  write.table(dmrs, paste0("DMR_", aname, ".bed"), row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
}

write.table(dmrsTreatWT, "DMR_treatWT.tsv", row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
write.table(dmrsPhenoNIR, "DMR_phenoNIR.tsv", row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
write.table(dmrsBatch2, "DMR_batch2.tsv", row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
write.table(dmrsInter2, "DMR_interaction2.tsv", row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)

showOneDMR(dmrsTreatWT[1,], bsf)




