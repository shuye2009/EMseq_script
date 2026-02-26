rm(list = ls())

library(GenomicPlot)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(ggplot2)

source("C:/PROJECTS/Repositories/EMseq_script/create_chrominfo_from_fasta.R")

feature_dir <- "C:/PROJECTS/resource/T2T_CHM13"
bedgraph_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/bedGraph"
fasta_file <- "C:/PROJECTS/resource/T2T_CHM13/chm13v2.0.fa.gz"

# create chromInfo for T2T_CHM13 using fasta file
if(!file.exists(file.path(feature_dir, "chromInfo.tsv"))) {
    chromInfo <- create_chrominfo_from_fasta(fasta_file)
    chromInfo <- chromInfo %>%
        mutate(chr = sub(" .*$", "", chrom)) %>%  # Remove everything after first space
        mutate(start = 0, end = length) %>%
        select(chr, start, end) %>%
        filter(chr != "chrM")
    write.table(chromInfo, file.path(feature_dir, "chromInfo.tsv"), sep = "\t", row.names = FALSE)
} else {
    chromInfo <- read.table(file.path(feature_dir, "chromInfo.tsv"), sep = "\t", header = TRUE)
}

bgImportParams <- GenomicPlot::setImportParams(genome = "hs1", chromInfo = chromInfo, val=4, skip=1)

center_file <- file.path(feature_dir, "chm13v2.0_centromere.bed")
names(center_file) <- "centromere"
bedfiles <- list.files(bedgraph_dir, pattern = "*-1_CpG.bedGraph.gz$", full.names = TRUE)
names(bedfiles) <- sub("_CpG.bedGraph.gz$", "", basename(bedfiles))

# Debug: Print file information
cat("Found", length(bedfiles), "bedGraph files:\n")
print(names(bedfiles))
cat("Center file:", center_file, "\n")
cat("ChromInfo structure:\n")
str(chromInfo)

# Try with a single file first to isolate the issue
if (length(bedfiles) > 0) {
  test_file <- bedfiles[1]
  cat("Testing with single file:", basename(test_file), "\n")
  
  # Plot methylation profiles around centromeres
  GenomicPlot::plot_region(
    queryFiles = bedfiles,
    centerFiles = center_file,
    importParams = bgImportParams,
    nbins = 100,
    fiveP = -1000,
    threeP = 1000,
    regionName = "Centromere",
    heatmap = TRUE,
    hw = c(8, 8),
    rmOutlier = 0,
    stranded = FALSE,
    transform = NA,
    smooth = TRUE,
    scale = FALSE,
    Ylab = "Methylation level",
    outPrefix = file.path(bedgraph_dir, "centromere_methylation_profile"),
    nc = 5,
    verbose = TRUE
  )
} else {
  cat("No bedGraph files found in", bedgraph_dir, "\n")
}

# Plot methylation profiles of enhancers
enhancer_file <- file.path(feature_dir, "ENCFF912FUA_MCF10A_element_gene_links_thresholded_Engreitz_T2T_6col.bed")
names(enhancer_file) <- "enhancer"

GenomicPlot::plot_locus(
  queryFiles = bedfiles,
  centerFiles = enhancer_file,
  txdb = NULL,
  ext = c(-500, 500),
  hl = c(0, 0),
  shade = TRUE,
  smooth = FALSE,
  importParams = bgImportParams,
  verbose = FALSE,
  binSize = 10,
  refPoint = "center",
  Xlab = "Distance from enhancer center (bp)",
  Ylab = "Coverage/base/gene",
  inputFiles = NULL,
  stranded = FALSE,
  heatmap = TRUE,
  scale = FALSE,
  outPrefix = file.path(feature_dir, "enhancer_methylation_profile"),
  rmOutlier = 0,
  transform = NA,
  statsMethod = "wilcox.test",
  heatRange = NULL,
  hw = c(8, 8),
  nc = 5
)

