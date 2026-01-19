#!/usr/bin/bash

script=C:/PROJECTS/Shane/Harding_240124/script/GenomicPlot_dmr.R

for COMPARISON in WT_NIR_vs_R172K_NIR WT_NIR_vs_WT_IR WT_IR_vs_R172K_IR R172K_NIR_vs_R172K_IR; do
  for CONTEXT in CG CHG CHH; do
    Rscript $script TRUE ERV $CONTEXT $COMPARISON
  done
done