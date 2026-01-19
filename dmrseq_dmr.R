rm(list = ls())

library(dmrseq)
library(bsseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)


home_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/TE_targeted"

em_dirs <- c("IR2Gy24h_vs_NIR", "IR10Gy24h_vs_NIR", "IR10Gy24h_vs_IR2Gy24h", "IR2Gy6d_vs_NIR", "IR10Gy6d_vs_NIR", "IR10Gy6d_vs_IR2Gy6d", "IR2Gy6d_vs_IR2Gy24h", "IR10Gy6d_vs_IR10Gy24h")

for (em_dir in em_dirs) {
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
      BPPARAM = BiocParallel::SnowParam(workers = 2)
    )

    # save dmr
    save(dmr, file = file.path(working_dir, "dmrseq_dmr.RData"))
  }else{
    load(file.path(working_dir, "dmrseq_dmr.RData"))
  }

  dmr_table <- as.data.frame(dmr)
  write_tsv(dmr_table, file.path(dmr_dir, "DMR_table.tsv"))
}
