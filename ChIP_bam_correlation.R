library(GenomicPlot)

setwd("C:\\PROJECTS\\Shane\\ChIPseq_240424\\02.Bam")

bamfiles <- list.files(pattern="\\.bam$")
names(bamfiles) <- gsub(".bam", "", bamfiles, fixed = TRUE)

grouping <- c("ChIP", "ChIP", "Input", "Input")
names(grouping) <- names(bamfiles)

plot_bam_correlation(
  bamFiles = bamfiles,
  binSize = 1e+06,
  outPrefix = NULL,
  importParams = setImportParams(genome="hg38", outRle = FALSE),
  grouping = grouping,
  verbose = TRUE,
  hw = c(8, 8),
  nc = 5
)

tmp <- handle_bam(bamfiles[1], importParams = setImportParams(genome="hg38", outRle = FALSE), verbose = TRUE)

library(GenomicAlignments)
flag <- scanBamFlag(isPaired = TRUE, isProperPair=T, isUnmappedQuery = F, hasUnmappedMate = F, isSecondMateRead = T)
param <- ScanBamParam(mapqFilter = 10, flag = flag)
ga <- readGAlignments(bamfiles[1], use.names = TRUE, param = param)
queryRegions <- unlist(grglist(ga))
score(queryRegions) <- 1
head(queryRegions)


if("NCBI" %in% seqlevelsStyle(queryRegions)){
  seqlevelsStyle(queryRegions) <- "UCSC"
}


seqlevels(queryRegions)
