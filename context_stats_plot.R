args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 1) {
  message("Usage: Rscript context_stats.R FALSE")
}
LOCAL = args[1] # TRUE or FALSE

wd <- ifelse(LOCAL, "C:/PROJECTS/Shane",
             "/cluster/projects/hardinggroup/Shuye/EMseq/Harding_240124/analysis")
setwd(wd)

library(dplyr)
library(ggplot2)

datapath <- ifelse(LOCAL, wd, 
                   "/cluster/projects/hardinggroup/Shuye/EMseq/Harding_240124/alignment/MethylDackel/report")

pdf("Methylation_level_distribution_by_CONTEXT.pdf")
for(CONTEXT in c("CpG", "CHG", "CHH")){
  
  queryfiles <- list.files(path = datapath, pattern = paste0(CONTEXT, "\\.bedGraph\\.gz$"))
  shortName <- gsub("\\.bedGraph\\.gz|_dox|_S\\d", "", queryfiles, fixed = FALSE)
  queryfiles <- file.path(datapath, queryfiles)
  names(queryfiles) <- shortName
  
  print(queryfiles)
  
  bg_dfs <- lapply(names(queryfiles), function(n){
    message("reading ", n)
    df <- read.delim(queryfiles[n], header = FALSE, skip = 1)
    colnames(df) <- c("chr", "start", "end", "methylation_level", "methylatedC", "unmethylatedC")
    df$treatment <- rep(n, nrow(df))
    return(df)
  })
  
  bg_df <- as.data.frame(do.call("rbind", bg_dfs))
  colnames(bg_df) <- c("chr", "start", "end", "methylation_level", "methylatedC", "unmethylatedC", "treatment")
  bg_df_limit <- bg_df %>%
    dplyr::filter(methylation_level > 0)
  
  p1 <- ggplot2::ggplot(bg_df) +
    geom_density(aes(x=methylation_level, color=treatment), linewidth = 1.2) +
    ggtitle(CONTEXT) +
    theme_classic()
  print(p1)
  
  p2 <- ggplot2::ggplot(bg_df_limit) +
    geom_density(aes(x=methylation_level, color=treatment), linewidth = 1.2) +
    ggtitle(paste0(CONTEXT, ">0")) +
    theme_classic()
  print(p2)
  
}
dev.off()

pdf("Methylation_level_distribution_by_TREATMENT.pdf")
for(t in c("R172K_dox_IR_S3", "R172K_dox_NIR_S4", "WT_dox_IR_S1", "WT_dox_NIR_S2")){
  queryfiles <- list.files(path = datapath, pattern = paste0("^", t, ".*\\.bedGraph\\.gz$"))
  shortName <- gsub(paste0("\\.bedGraph\\.gz|_|", t), "", queryfiles, fixed = FALSE)
  queryfiles <- file.path(datapath, queryfiles)
  names(queryfiles) <- shortName
  
  print(queryfiles)
  bg_dfs <- bg_df <- NULL
  
  bg_dfs <- lapply(names(queryfiles), function(n){
    message("reading ", n)
    df <- read.delim(queryfiles[n], header = FALSE, skip = 1)
    colnames(df) <- c("chr", "start", "end", "methylation_level", "methylatedC", "unmethylatedC")
    df$context <- rep(n, nrow(df))
    return(df)
  })
  
  bg_df <- as.data.frame(do.call("rbind", bg_dfs))
  colnames(bg_df) <- c("chr", "start", "end", "methylation_level", "methylatedC", "unmethylatedC", "context")
  bg_df_limit <- bg_df %>%
    dplyr::filter(methylation_level > 0)
  
  p1 <- ggplot2::ggplot(bg_df) +
    geom_density(aes(x=methylation_level, color=context), linewidth = 1.2) +
    ggtitle(t) +
    theme_classic()
  print(p1)
  
  p2 <- ggplot2::ggplot(bg_df_limit) +
    geom_density(aes(x=methylation_level, color=treatment), linewidth = 1.2) +
    ggtitle(paste0(t, ">0")) +
    theme_classic()
  print(p2)
}
dev.off()
