rm(list = ls())

library(dmrseq)
library(bsseq)
library(tidyr)
library(readr)
library(tibble)
library(BiocParallel)

# get home directory from command line argument
home_dir <- commandArgs(trailingOnly = TRUE)[1]
# get em_dir from command line argument
em_dir <- commandArgs(trailingOnly = TRUE)[2]

working_dir <- file.path(home_dir, em_dir, "RData")
dmr_dir <- file.path(home_dir, em_dir, "DMR")
if (!dir.exists(dmr_dir)) {
    dir.create(dmr_dir)
}
# load filtered bismark object: bs.filtered which is saved as bismark.RData
load(file.path(working_dir, "bismark.RData"))

print(pData(bs.filtered))
if(!file.exists(file.path(working_dir, "dmrseq_dmr.RData"))) {
    # call dmr using dmrseq
    dmr <- dmrseq::dmrseq(
      bs = bs.filtered,
      cutoff = 0.05,
      minNumRegion = 5,
      maxPerms = 100,
      testCovariate = "Group",
      adjustCovariate = NULL,
      matchCovariate = NULL,
      BPPARAM = BiocParallel::MulticoreParam(workers = 10)
    )

    # save dmr
    save(dmr, file = file.path(working_dir, "dmrseq_dmr.RData"))
}else{
    load(file.path(working_dir, "dmrseq_dmr.RData"))
}

dmr_table <- as.data.frame(dmr)
write_tsv(dmr_table, file.path(dmr_dir, "DMR_table.tsv"))
