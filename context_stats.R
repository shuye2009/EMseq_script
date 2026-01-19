args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 1) {
  message("Usage: Rscript context_stats.R FALSE")
}
LOCAL = args[1] # TRUE or FALSE

wd <- ifelse(LOCAL, "C:/PROJECTS/Shane/Harding_240124/local",
             "/cluster/projects/hardinggroup/Shuye/EMseq/Harding_240124/analysis")
setwd(wd)

library(dplyr)
library(ggplot2)

datapath <- ifelse(LOCAL, "C:/PROJECTS/Shane/Harding_240124/h4h/data",
                   "/cluster/projects/hardinggroup/Shuye/EMseq/Harding_240124/alignment/MethylDackel/report")


for(CONTEXT in c("CpG", "CHG", "CHH")){
  
  queryfiles <- list.files(path = datapath, pattern = paste0(CONTEXT, "\\.bedGraph\\.gz$"))
  shortName <- gsub("\\.bedGraph\\.gz|_dox|_S\\d", "", queryfiles, fixed = FALSE)
  queryfiles <- file.path(datapath, queryfiles)
  names(queryfiles) <- shortName
  
  print(queryfiles)
  sink(paste0(CONTEXT, "_methylation_stats.txt"))
  for (n in names(queryfiles)){
    message("reading ", n)
    df <- read.delim(queryfiles[n], header = FALSE, skip = 1)
    colnames(df) <- c("chr", "start", "end", "methylation_level", "methylatedC", "unmethylatedC")
    df_limit <- df %>%
      dplyr::filter(methylatedC > 0)
    
    print(n)
    print(summary(df))
    print(summary(df_limit))
    
  }
  sink()
}



