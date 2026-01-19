## ----bioconductor, message=FALSE, warning=FALSE, eval=FALSE-------------------
if (!require("BiocManager"))
 	install.packages("BiocManager")
 BiocManager::install("DMRcate")

## ----libr, message=FALSE, warning=FALSE---------------------------------------
library(DMRcate)

## ----loadeh, message=FALSE----------------------------------------------------
library(ExperimentHub)
eh <- ExperimentHub()
bis_1072 <- eh[["EH1072"]]
bis_1072
colnames(bis_1072)

## ----bisphen------------------------------------------------------------------
bsseq::pData(bis_1072) <- data.frame(replicate=gsub(".*-", "", colnames(bis_1072)),
                                     tissue=substr(colnames(bis_1072), 1, 
                                                   nchar(colnames(bis_1072))-3), 
                                     row.names=colnames(bis_1072))
colData(bis_1072)$tissue <- gsub("-", "_", colData(bis_1072)$tissue)
as.data.frame(colData(bis_1072))

## ----changeseqlevs------------------------------------------------------------
bis_1072 <- renameSeqlevels(bis_1072, mapSeqlevels(seqlevels(bis_1072), "UCSC"))

## ----chr19filter--------------------------------------------------------------
bis_1072 <- bis_1072[seqnames(bis_1072)=="chr19",]
bis_1072

## ----bsdesign, message=FALSE--------------------------------------------------
tissue <- factor(pData(bis_1072)$tissue)
tissue <- relevel(tissue, "Liver_Treg")

#Regular matrix design
design <- model.matrix(~tissue)
colnames(design) <- gsub("tissue", "", colnames(design))
colnames(design)[1] <- "Intercept"
rownames(design) <- colnames(bis_1072)
design

#Methylation matrix design
methdesign <- edgeR::modelMatrixMeth(design)
methdesign

## ----fitBSseq-----------------------------------------------------------------
cont.mat <- limma::makeContrasts(treg_vs_tcon=Lymph_N_Treg-Lymph_N_Tcon,
                                 fat_vs_ln=Fat_Treg-Lymph_N_Treg,
                                 skin_vs_ln=Skin_Treg-Lymph_N_Treg,
                                 fat_vs_skin=Fat_Treg-Skin_Treg,
                                 levels=methdesign)
cont.mat

## ----sequencingannotate-------------------------------------------------------
seq_annot <- sequencing.annotate(bis_1072, methdesign, all.cov = TRUE, 
                                 contrasts = TRUE, cont.matrix = cont.mat, 
                                 coef = "treg_vs_tcon", fdr=0.05)
seq_annot

## ----seqdmrcate---------------------------------------------------------------
dmrcate.res <- dmrcate(seq_annot, C=2, min.cpgs = 5)
dmrcate.res
treg_vs_tcon.ranges <- extractRanges(dmrcate.res, genome="mm10")
treg_vs_tcon.ranges

## ----seqDMRplot1, message=FALSE-----------------------------------------------
cols <- as.character(plyr::mapvalues(tissue, unique(tissue), 
                                     c("darkorange", "maroon", "blue", 
                                       "black", "magenta")))
names(cols) <- tissue

DMR.plot(treg_vs_tcon.ranges, dmr = 1, 
         CpGs=bis_1072[,tissue %in% c("Lymph_N_Tcon", "Lymph_N_Treg")], 
         phen.col = cols[tissue %in% c("Lymph_N_Tcon", "Lymph_N_Treg")], 
         genome="mm10")

## ----fatskin------------------------------------------------------------------
seq_annot <- sequencing.annotate(bis_1072, methdesign, all.cov = TRUE, 
                                 contrasts = TRUE, cont.matrix = cont.mat, 
                                 coef = "fat_vs_skin", fdr=0.05)

## ----redefinethresh-----------------------------------------------------------
seq_annot <- changeFDR(seq_annot, 0.25)

## ----dmrsfatskin--------------------------------------------------------------
dmrcate.res <- dmrcate(seq_annot, C=2, min.cpgs = 5)
fat_vs_skin.ranges <- extractRanges(dmrcate.res, genome="mm10")

## ----seqDMRplot2, message=FALSE-----------------------------------------------
cols
DMR.plot(fat_vs_skin.ranges, dmr = 1, CpGs=bis_1072, phen.col = cols, genome="mm10")

## ----sessionInfo--------------------------------------------------------------
sessionInfo()
