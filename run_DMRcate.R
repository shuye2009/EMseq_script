library(DMRcate)
library(bsseq)
library(comethyl)
library(dplyr)

setwd("C:\\PROJECTS\\Shane\\Harding_combined\\local\\DMRcate_res")
cores <- 1
files = list.files(path = getwd(), pattern = "*.CG_report.txt.gz")

print(glue::glue("Reading cytosine reports..."))
if(!file.exists("bs.RDS")){
  bs <- bsseq::read.bismark(files = files,
                            #colData = names,
                            rmZeroCov = FALSE,
                            strandCollapse = TRUE,
                            verbose = TRUE,
                            BPPARAM = BiocParallel::SnowParam(workers = cores, progressbar = TRUE), # BPPARAM # bpparam() # MulticoreParam(workers = cores, progressbar = TRUE)
                            nThread = 1) # 1L # nThread
  #bs <- bis_1072
  
  saveRDS(bs, "bs.RDS")
}else{
  bs <- readRDS("bs.RDS")
}


colnames(bs) <- sapply((strsplit(colnames(bs), split=".", fixed = T)), function(x)x[1])
#bs <- bs[, !grepl("-0", colnames(bs))]

## ----changeseqlevs------------------------------------------------------------
#bs <- renameSeqlevels(bs, mapSeqlevels(seqlevels(bs), "UCSC"))

## ----chrfilter--------------------------------------------------------------
standardChr <- c(paste0("chr", seq.int(22)), "chrX", "chrY")
bs <- bs[seqnames(bs) %in% standardChr, ]
seqlevels(bs) <- standardChr


if(!file.exists("bss.RDS")){
  bss <- BSmooth(bs, BPPARAM = BiocParallel::SnowParam(workers = cores, progressbar = T))
  
  saveRDS(bss, "bss.RDS")
}else{
  bss <- readRDS("bss.RDS")
}

bss <- bss[seqnames(bss) %in% standardChr, ]
seqlevels(bss) <- standardChr

rm(bs)
gc()
## ----filtering by coverage -------------------------------------------------------
if(0){
bss_df <- comethyl::getCpGtotals(
  bss,
  cov = seq(10, 20, 1),
  perSample = seq(0.5, 1, 0.05),
  save = TRUE,
  file = "CpG_Totals.txt",
  verbose = TRUE
)

comethyl::plotCpGtotals(
  bss_df,
  nBreaks = 4,
  legend.position = c(1.08, 0.73),
  save = TRUE,
  file = "CpG_Totals.pdf",
  width = 11,
  height = 4.25,
  verbose = TRUE
)


bss_filtered <- comethyl::filterCpGs(bss, cov=15, perSample=1.0, save = TRUE,
                           file = "Filtered_bss_cov10_pct100.RDS",)

cov <- bsseq::getCoverage(bss)
meth <- bsseq::getMeth(bss)
}
## ----bsdesign, message=FALSE--------------------------------------------------

if(1){
  sample_info <- data.frame(replicate=gsub(".*-", "", colnames(bss)),
                            group=factor(sapply((strsplit(colnames(bss), split="-")), 
                                             function(x)paste(x[1:2], collapse="_"))),
                            batch=factor(ifelse(grepl("-0", colnames(bss)), 1, 2)),
                            row.names=colnames(bss))
  
  pData(bss) <- sample_info
  colData(bss)
  pData(bss)
  #Regular matrix design
  design <- model.matrix(~0+group+batch, data=pData(bss))
  #design1 <- model.matrix(~as.factor(group), data=pData(bs))
  colnames(design) <- gsub("group", "", colnames(design))
  #colnames(design)[1] <- "Intercept"
  rownames(design) <- colnames(bss)
  colnames(design) <- gsub(":", "_", colnames(design))
  design
  
  #Methylation matrix design
  methdesign <- edgeR::modelMatrixMeth(design)
  methdesign
  
  ## ----fitBSseq-----------------------------------------------------------------
  cont.mat <- limma::makeContrasts(Irradiation=(WT_IR+R172K_IR)/2-(WT_NIR+R172K_NIR)/2,
                                   Mutation=(R172K_IR+R172K_NIR)/2-(WT_IR+WT_NIR)/2,
                                   Rad_in_Mut=R172K_IR-R172K_NIR,
                                   Rad_in_WT=WT_IR-WT_NIR,
                                   Mut_in_Rad=R172K_IR-WT_IR,
                                   Mut_in_NIR=R172K_NIR-WT_NIR,
                                   levels=methdesign)
  cont.mat
}

# use ir/mut design ####
if(0){
  sample_info <- data.frame(replicate=gsub(".*-", "", colnames(bs)),
                           ir=factor(sapply((strsplit(colnames(bs), split="-")), 
                                     function(x)x[2]), levels=c("NIR", "IR")),
                           mut=factor(sapply((strsplit(colnames(bs), split="-")), 
                                     function(x)x[1]), levels=c("WT", "R172K")),
                           batch=factor(ifelse(grepl("-0", colnames(bs)), 1, 2)),
                           row.names=colnames(bs))
  
  pData(bss) <- sample_info
  colData(bss)
  pData(bss)
  
  design <- model.matrix(~ir+mut+batch+ir*mut, data=pData(bss))
  colnames(design) <- gsub("ir|mut", "", colnames(design))
  colnames(design)[1] <- "Intercept"
  colnames(design) <- gsub(":", "_", colnames(design))
  design
  
  #Methylation matrix design
  methdesign <- edgeR::modelMatrixMeth(design)
  methdesign
  
  ## ----fitBSseq-----------------------------------------------------------------
  cont.mat <- limma::makeContrasts(Irradiation=IR,
                                   Mutation=R172K,
                                   Interaction=IR_R172K,
                                   Batch=batch2,
                                   levels=methdesign)
  cont.mat
}

## Irradiation DMR --------------------------------------------------------------
if(file.exists("irradiation_DMR.RDS")){
  irradiation.ranges <- readRDS("irradiation_DMR.RDS")
}else{
  seq_annotI <- sequencing.annotate(bss, methdesign, all.cov = TRUE, 
                                    contrasts = TRUE, cont.matrix = cont.mat, 
                                    coef = "Irradiation", fdr=0.05)
  
  min(seq_annotI@ranges$rawpval)
  
  length(seq_annotI@ranges[seq_annotI@ranges$rawpval<0.001])
  
  dmrcate.res <- dmrcate(seq_annotI, C=2, min.cpgs = 5, pcutoff = 0.9, lambda = 300)
  dmrcate.res
  irradiation.ranges <- extractRanges(dmrcate.res, genome="hg38")
  saveRDS(irradiation.ranges, "irradiation_DMR.RDS")
}

ir_filtered_dmr <- irradiation.ranges[irradiation.ranges$min_smoothed_fdr<0.001 & abs(irradiation.ranges$meandiff) > 0.25]

group <- pData(bss)$group
cols <- as.character(plyr::mapvalues(group, unique(group), 
                                     c("darkorange", "maroon", "blue", 
                                       "magenta")))
names(cols) <- group

pdf("IR_DMRs1_10.pdf")
for (i in 1:10){
  DMR.plot(ir_filtered_dmr, dmr = i, 
         CpGs=bss, flank = 2500,
         phen.col = cols, 
         genome="hg38")
}
dev.off()

hist(ir_filtered_dmr$meandiff)

## ----Mutation ------------------------------------------------------------------
if(file.exists("mutation_DMR.RDS")){
  mutation.ranges <- readRDS("mutation_DMR.RDS")
}else{
  seq_annotM <- sequencing.annotate(bss, methdesign, all.cov = TRUE, 
                                   contrasts = TRUE, cont.matrix = cont.mat, 
                                   coef = "Mutation", fdr=0.05)
  
  dmrcate.resM <- dmrcate(seq_annotM, C=2, min.cpgs = 5,  pcutoff = 0.9, lambda = 300)
  mutation.ranges <- extractRanges(dmrcate.resM, genome="hg38")
  
  saveRDS(mutation.ranges, "mutation_DMR.RDS")
}


mut_filtered_dmr <- mutation.ranges[mutation.ranges$min_smoothed_fdr<0.001 & abs(mutation.ranges$meandiff) > 0.8]
#names(mut_filtered_dmr) <- paste0("mut", seq_along(mut_filtered_dmr))
#brca2 <- mut_filtered_dmr[grepl("BRCA2", mut_filtered_dmr$overlapping.genes)]
pdf("Mutation_DMRs1_10.pdf")
for (i in 1:10){
  DMR.plot(mut_filtered_dmr, dmr = i, 
         CpGs=bss,
         phen.col = cols, 
         genome="hg38")
}
dev.off()
hist(mut_filtered_dmr$meandiff)

## ----Mutation in NIR ------------------------------------------------------------------
if(file.exists("mutation_in_NIR_DMR.RDS")){
  mutation.rangesNIR <- readRDS("mutation_in_NIR_DMR.RDS")
}else{
  seq_annotM_NIR <- sequencing.annotate(bss, methdesign, all.cov = TRUE, 
                                    contrasts = TRUE, cont.matrix = cont.mat, 
                                    coef = "Mut_in_NIR", fdr=0.05)
  
  
  dmrcate.resM_NIR <- dmrcate(seq_annotM_NIR, C=2, min.cpgs = 5,  pcutoff = 0.9, lambda = 300)
  mutation.rangesNIR <- extractRanges(dmrcate.resM_NIR, genome="hg38")
  
  saveRDS(mutation.rangesNIR, "mutation_in_NIR_DMR.RDS")
}

mutNIR_filtered_dmr <- mutation.rangesNIR[mutation.rangesNIR$min_smoothed_fdr<0.001 & abs(mutation.rangesNIR$meandiff) > 0.5]

pdf("Mutation_NIR_DMRs1_10.pdf")
for (i in 1:10){
  DMR.plot(mutNIR_filtered_dmr, dmr = i, 
           CpGs=bss[,group %in% c("R172K_NIR", "WT_NIR")], 
           phen.col = cols[group %in% c("R172K_NIR", "WT_NIR")], 
           genome="hg38")
}
dev.off()
## ----Irradiation in WT -----------------------------------------------
if(file.exists("irradiation_in_WT_DMR.RDS")){
  irradiation.rangesWT <- readRDS("irradiation_in_WT_DMR.RDS")
}else{
  seq_annotI_WT <- sequencing.annotate(bss, methdesign, all.cov = TRUE, 
                                    contrasts = TRUE, cont.matrix = cont.mat, 
                                    coef = "Rad_in_WT", fdr=0.05)
  
  
  dmrcate.resI_WT <- dmrcate(seq_annotI_WT, C=2, min.cpgs = 5, pcutoff = 0.9, lambda = 300)
  dmrcate.resI_WT
  irradiation.rangesWT <- extractRanges(dmrcate.resI_WT, genome="hg38")
  saveRDS(irradiation.rangesWT, "irradiation_in_WT_DMR.RDS")
}

irWT_filtered_dmr <- irradiation.rangesWT[irradiation.rangesWT$min_smoothed_fdr<0.001 & abs(irradiation.rangesWT$meandiff) > 0.25]

pdf("Irradiation_WT_DMRs1_10.pdf")
for (i in 1:10){
  DMR.plot(irWT_filtered_dmr, dmr = i, 
           CpGs=bss[,group %in% c("WT_IR", "WT_NIR")], 
           phen.col = cols[group %in% c("WT_IR", "WT_NIR")], 
           genome="hg38")
}
dev.off()

# functional enrichment #### --------------------------------------------------------
dmg <- function(filtered_dmr){
  
  hyper <- filtered_dmr[filtered_dmr$meandiff>0]
  hyper_gene <- hyper$overlapping.genes
  hyper_gene <- hyper_gene[!is.na(hyper_gene)]
  hyper_gene <- unlist(strsplit(hyper_gene, split=", "))
  hyper_gene <- unique(noquote(hyper_gene))
  
  hypo <- filtered_dmr[filtered_dmr$meandiff<0]
  hypo_gene <- hypo$overlapping.genes
  hypo_gene <- hypo_gene[!is.na(hypo_gene)]
  hypo_gene <- unlist(strsplit(hypo_gene, split=", "))
  hypo_gene <- unique(noquote(hypo_gene))
  
  return(list(hyper=hyper_gene, hypo=hypo_gene))
}
source("C:/PROJECTS/Shane/script/GO_KEGG_enrichment.R")
sample_pairs <- c("IR-NIR", "R172K-WT", "WT_IR-WT_NIR", "R172K_NIR-WT_NIR")
comparisons <- c("Irradiation", "Mutation", "IrradiationWT", "MutationNIR")
directions <- c("hyper", "hypo")

gene_list <- list(dmg(ir_filtered_dmr), dmg(mut_filtered_dmr), dmg(irWT_filtered_dmr), dmg(mutNIR_filtered_dmr))
names(gene_list) <- sample_pairs

compareCluster_enrichment(gene_list, comparisons, sample_pairs, directions, c="DMRcate", 8, 8)

# overlap with cgmaptools #### ----------------------------------------------------------
library(rtracklayer)
library(GenomicPlot)

export.bed(ir_filtered_dmr[ir_filtered_dmr$meandiff>0], "ir_dmr_hyper.bed")
export.bed(ir_filtered_dmr[ir_filtered_dmr$meandiff<0], "ir_dmr_hypo.bed")
export.bed(mut_filtered_dmr[mut_filtered_dmr$meandiff>0], "mut_dmr_hyper.bed")
export.bed(mut_filtered_dmr[mut_filtered_dmr$meandiff<0], "mut_dmr_hypo.bed")

export.bed(ir_filtered_dmr[irWT_filtered_dmr$meandiff>0], "irWT_dmr_hyper.bed")
export.bed(ir_filtered_dmr[irWT_filtered_dmr$meandiff<0], "irWT_dmr_hypo.bed")
export.bed(mut_filtered_dmr[mutNIR_filtered_dmr$meandiff>0], "mutNIR_dmr_hyper.bed")
export.bed(mut_filtered_dmr[mutNIR_filtered_dmr$meandiff<0], "mutNIR_dmr_hypo.bed")

cgmaptools_ir_hyper <- "C:/PROJECTS/Shane/Harding_combined/h4h/cgmaptools/cgmaptools_dmr/WT-NIR_vs_WT-IR_dmr.CG/HOMER/DMRs_hyper.bed"
cgmaptools_ir_hypo <- "C:/PROJECTS/Shane/Harding_combined/h4h/cgmaptools/cgmaptools_dmr/WT-NIR_vs_WT-IR_dmr.CG/HOMER/DMRs_hypo.bed"
cgmaptools_mut_hyper <- "C:/PROJECTS/Shane/Harding_combined/h4h/cgmaptools/cgmaptools_dmr/WT-NIR_vs_R172K-NIR_dmr.CG/HOMER/DMRs_hyper.bed"
cgmaptools_mut_hypo <- "C:/PROJECTS/Shane/Harding_combined/h4h/cgmaptools/cgmaptools_dmr/WT-NIR_vs_R172K-NIR_dmr.CG/HOMER/DMRs_hypo.bed"

dmrseq_ir_hyper <- "C:/PROJECTS/Shane/Harding_combined/h4h/DMRichR_treat/HOMER/DMRs_hyper.bed"
dmrseq_ir_hypo <- "C:/PROJECTS/Shane/Harding_combined/h4h/DMRichR_treat/HOMER/DMRs_hypo.bed"
dmrseq_mut_hyper <- "C:/PROJECTS/Shane/Harding_combined/h4h/DMRichR_mutant/HOMER/DMRs_hyper.bed"
dmrseq_mut_hypo <- "C:/PROJECTS/Shane/Harding_combined/h4h/DMRichR_mutant/HOMER/DMRs_hypo.bed"

bedimportParams <- setImportParams(
  offset = 0, fix_width = 0, fix_point = "center", norm = FALSE,
  useScore = FALSE, outRle = FALSE, useSizeFactor = FALSE, genome = "hg38")

bedfiles <- c(#"ir_dmr_hyper.bed", "ir_dmr_hypo.bed", "mut_dmr_hyper.bed", "mut_dmr_hypo.bed",
              "irWT_dmr_hyper.bed", "irWT_dmr_hypo.bed", "mutNIR_dmr_hyper.bed", "mutNIR_dmr_hypo.bed",
              cgmaptools_ir_hyper, cgmaptools_ir_hypo, cgmaptools_mut_hyper, cgmaptools_mut_hypo,
              dmrseq_ir_hyper, dmrseq_ir_hypo, dmrseq_mut_hyper, dmrseq_mut_hypo)
names(bedfiles) <- c("dmrcate_ir_hyper", "dmrcate_ir_hypo", "dmrcate_mut_hyper", "dmrcate_mut_hypo",
                     #"dmrcate_irWT_hyper", "dmrcate_irWT_hypo", "dmrcate_mutNIR_hyper", "dmrcate_mutNIR_hypo",
                     "cgmaptools_ir_hyper", "cgmaptools_ir_hypo", "cgmaptools_mut_hyper", "cgmaptools_mut_hypo",
                     "dmrseq_ir_hyper", "dmrseq_ir_hypo", "dmrseq_mut_hyper", "dmrseq_mut_hypo")


for(type in c("ir", "mut")){
  for(direction in c("hyper", "hypo")){
    combo <- paste(type, direction, sep="_")
    plot_overlap_bed(
      bedList = bedfiles[paste(c("dmrcate", "cgmaptools", "dmrseq"), combo, sep="_")], 
      importParams = bedimportParams, pairOnly = FALSE,
      stranded = FALSE, outPrefix = paste("overlap_for", combo, sep="_")
    )
    
  }
}


## ----sessionInfo--------------------------------------------------------------
sessionInfo()
