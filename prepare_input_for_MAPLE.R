#!/usr/bin/env Rscript

# Script: plot_regions_bsseq.R
# Purpose: Plot methylation profiles using bsseq plotRegions function
# Input: Bismark RData file and BED regions file
# Output: Methylation profile plots

# Clear workspace
rm(list = ls())

# Load required libraries
library(bsseq)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(IRanges)
library(S4Vectors)
library(SummarizedExperiment)
library(GenomeInfoDb)
library(methylKit)
library(parallel)
library(R.utils)

# Direct liftOver: T2T-CHM13 -> hg19
message("Setting up liftOver: T2T-CHM13 -> hg19...")
chain_dir <- "C:/PROJECTS/resource/T2T_CHM13"

# T2T-CHM13 to hg19 chain from UCSC
chain_gz <- file.path(chain_dir, "Hs1ToHg19.over.chain.gz")
chain_file <- file.path(chain_dir, "Hs1ToHg19.over.chain")

if (!file.exists(chain_file)) {
    if (!file.exists(chain_gz)) {
        message("Downloading T2T-CHM13 to hg19 chain file...")
        download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/Hs1ToHg19.over.chain.gz",
                      chain_gz, mode = "wb")
    }
    message("Decompressing chain file...")
    gunzip(chain_gz, destname = chain_file, remove = FALSE)
}

# Load chain file
message("Loading chain file...")
chain <- import.chain(chain_file)
message("Loaded T2T-CHM13->hg19 chain: ", length(chain), " chains")

bismark_file <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/Enhancer_targeted/IR10Gy6d_vs_NIR/RData/bismark.RData"
output_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/MAPLE_input"
comparison_tag <- "IR10Gy6d_vs_NIR"


# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# Check if input files exist
if (!file.exists(bismark_file)) {
  stop("Bismark RData file not found: ", bismark_file)
}

cat("Loading Bismark data from:", bismark_file, "\n")
# Load the Bismark RData object
load(bismark_file)

# Check what objects are loaded
cat("Objects in RData file:\n")
print(ls())

# Use the specific bs.filtered object
if (exists("bs.filtered")) {
  bsseq_obj <- bs.filtered
  cat("Using bs.filtered object\n")
} else {
  stop("bs.filtered object not found in RData file")
}

cat("BSseq object created successfully\n")
print(bsseq_obj)

# output format for MAPLE
# "",GSM6425155,GSM6425157,GSM6425173,GSM6425633,GSM6425989,GSM6426061
# cg00050873,,,,,,
# cg00212031,,0.0365512792947753,0.0493803253292022,,0.160498220640569,
# cg00213748,,0.813299232736573,0.911537167843929,,0.95327868852459,

beta_values <- getMeth(bsseq_obj, type = "raw")
head(beta_values)

pos_data <- GenomicRanges::granges(bsseq_obj)
print(head(pos_data))
  
# liftOver: T2T-CHM13 -> hg19
pos_data_hg19_list <- liftOver(pos_data, chain)

# Track which positions successfully lifted (some may fail)
keep_idx <- lengths(pos_data_hg19_list) > 0

# Unlist and keep only successful liftOvers
pos_data_hg19 <- unlist(pos_data_hg19_list)

# Fix seqlevels ordering
seqlevels(pos_data_hg19) <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

# CRITICAL: Subset beta_values to match positions that lifted successfully
beta_values <- beta_values[keep_idx, ]

message("Positions before liftOver: ", length(pos_data))
message("Positions after liftOver: ", length(pos_data_hg19))
message("Positions dropped: ", sum(!keep_idx))

print(head(pos_data_hg19))
  
# use pos_data to create the rownames
rownames(beta_values) <- paste0(seqnames(pos_data_hg19), "_", start(pos_data_hg19))
print(head(rownames(beta_values)))
  
# Get sample information
meta_data <- data.frame(sample_id = colnames(beta_values))
rownames(meta_data) <- colnames(beta_values)
print(head(meta_data))

write.csv(beta_values, file.path(output_dir, "beta_values.csv"), row.names = TRUE)
write.csv(meta_data, file.path(output_dir, "sample_metadata.csv"), row.names = TRUE)

