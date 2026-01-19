rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(readr)
library(tibble)
library(ggExtra)

centromere_file <- "C:\\PROJECTS\\resource\\T2T_CHM13\\chm13v2.0_censat_v2.1.bed"
centromere_data <- read_tsv(centromere_file, col_names = FALSE, col_types = cols(), skip = 1)
head(centromere_data)

# only keep the first 6 columns
centromere_data <- centromere_data[, 1:6]
colnames(centromere_data) <- c("chr", "start", "end", "name", "score", "strand")

# replace strand with * for all rows
centromere_data$strand <- "*"   

# extract satellite names and chromosomes from column 4 like "hor_5_6(S1C1/5/19H1L)", where "hor" is the satellite name, "5" is the chromosome, and "6" is the bin number
sat_names <- sapply(centromere_data[[4]], function(x) strsplit(x, "_")[[1]][1])
sat_chrs <- sapply(centromere_data[[4]], function(x) strsplit(x, "_")[[1]][2])

centromere_data$sat_name <- sat_names
centromere_data$length <- centromere_data$end - centromere_data$start
centromere_data$chr <- factor(centromere_data$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY"))

# order by chr names using numeric values followed by X and Y and sort by start
centromere_data <- centromere_data[order(centromere_data$chr, centromere_data$start), ]

centromere_range <- read_tsv("C:\\PROJECTS\\resource\\T2T_CHM13\\chm13v2.0_centromere.bed", col_names = FALSE)
colnames(centromere_range) <- c("chr", "cen_start", "cen_end")

# Create custom color palette for sat_names with 'hor' as red
all_sat_names <- unique(centromere_data$sat_name)
# Define a palette of 12 highly distinguishable colors
distinct_colors <- c(
  "#E41A1C",  # Red
  "#377EB8",  # Blue
  "#4DAF4A",  # Green
  "#984EA3",  # Purple
  "#FF7F00",  # Orange
  "#FFFF33",  # Yellow
  "#A65628",  # Brown
  "#F781BF",  # Pink
  "#00CED1",  # Turquoise
  "#8B4513",  # Saddle Brown
  "#32CD32",  # Lime Green
  "#700b59"   # Crimson
)

# Assign colors: 'hor' gets red, others get remaining colors
sat_colors <- rep(distinct_colors, length.out = length(all_sat_names))
names(sat_colors) <- all_sat_names

# Ensure 'hor' gets red color if it exists
if("hor" %in% all_sat_names) {
  sat_colors["hor"] <- "#E41A1C"  # Red
  # Redistribute other colors to remaining satellites
  other_sats <- setdiff(all_sat_names, "hor")
  sat_colors[other_sats] <- distinct_colors[2:(length(other_sats) + 1)]
}

# load data
wd <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/RNAseq/deseq2/single_factor_analysis_transcript_level"
setwd(wd)
EM_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/TE_targeted"


# list directory with pattern "vs"
#dirs <- list.files(path = wd, pattern = "_vs_", full.names = FALSE)
dirs <- c("2Gy24h_vs_NIR", "10Gy24h_vs_NIR", "10Gy24h_vs_2Gy24h", "2Gy6d_vs_NIR", "10Gy6d_vs_NIR", "10Gy6d_vs_2Gy6d", "2Gy6d_vs_2Gy24h", "10Gy6d_vs_10Gy24h")
print(dirs)
em_dirs <- c("IR2Gy24h_vs_NIR", "IR10Gy24h_vs_NIR", "IR10Gy24h_vs_IR2Gy24h", "IR2Gy6d_vs_NIR", "IR10Gy6d_vs_NIR", "IR10Gy6d_vs_IR2Gy6d", "IR2Gy6d_vs_IR2Gy24h", "IR10Gy6d_vs_IR10Gy24h")
names(em_dirs) <- dirs
# Initialize master dataframes
master_sat_df <- data.frame()
master_chr_df <- data.frame()
master_chr_sat_up_df <- data.frame()
master_chr_sat_down_df <- data.frame()

# loop through each directory and load data
for(dir in dirs){
   # read in DESeq2 results for RNAseq
   sat_file <- file.path(wd, dir, paste0(dir,  "_DESeq2_results_SAT.csv"))
   sat_data <- read_csv(sat_file)
   
   # read in targeted differential analysis results for EMseq
   em_file <- file.path(EM_dir, em_dirs[[dir]], "Targeted",  "all_diff_chm13v2.0_censat_v2.1_simplified.tab")
   em_data <- read_tsv(em_file) %>%
     mutate(empvalue = pvalue) %>%
     select(-pvalue)
   
   # extract up and down regulated transcripts
   up_genes <- sat_data %>% 
     filter(log2FoldChange > 1 & padj < 0.05) %>%
     pull(transcript_id)
   
   down_genes <- sat_data %>% 
     filter(log2FoldChange < -1 & padj < 0.05) %>%
     pull(transcript_id)

  # extract satellite names and chromosomes from transcript_id like "hor_5_6(S1C1/5/19H1L)", where "hor" is the satellite name, "5" is the chromosome, and "6" is the bin number
  
  # Extract satellite names for up and down regulated genes
  up_sat_names <- character(0)
  up_sat_chrs <- character(0)
  down_sat_names <- character(0)
  down_sat_chrs <- character(0)
  
  # extract satellite names and chromosomes from transcript_id, using strsplit because chromosome names may not be digits
  if(length(up_genes) > 0){
    up_sat_names <- sapply(up_genes, function(x) strsplit(x, "_")[[1]][1])
    up_sat_chrs <- paste0("chr", sapply(up_genes, function(x) strsplit(x, "_")[[1]][2]))
  }
  
  if(length(down_genes) > 0){
    down_sat_names <- sapply(down_genes, function(x) strsplit(x, "_")[[1]][1])
    down_sat_chrs <- paste0("chr", sapply(down_genes, function(x) strsplit(x, "_")[[1]][2]))
  }
  
  # Create dataframes for chromosome-satellite combinations
  if(length(up_sat_chrs) > 0){
    chr_sat_up_df <- data.frame(
      sat_chr = up_sat_chrs,
      sat_name = up_sat_names,
      comparison = dir
    ) %>%
      count(sat_chr, sat_name, comparison, name = "count") %>%
      mutate(sat_chr = factor(sat_chr, levels = c(paste0("chr", 1:22), "chrX", "chrY")))
    
    master_chr_sat_up_df <- bind_rows(master_chr_sat_up_df, chr_sat_up_df)
  }
  
  if(length(down_sat_chrs) > 0){
    chr_sat_down_df <- data.frame(
      sat_chr = down_sat_chrs,
      sat_name = down_sat_names,
      comparison = dir
    ) %>%
      count(sat_chr, sat_name, comparison, name = "count") %>%
      mutate(sat_chr = factor(sat_chr, levels = c(paste0("chr", 1:22), "chrX", "chrY")))
    
    master_chr_sat_down_df <- bind_rows(master_chr_sat_down_df, chr_sat_down_df)
  }
  
  # Create consolidated dataframe for satellite names
  all_sat_names <- sat_names
  
  if(length(all_sat_names) > 0){
    # Count up-regulated satellites
    up_counts <- data.frame(sat_name = up_sat_names) %>%
      count(sat_name, name = "up_count")
    
    # Count down-regulated satellites
    down_counts <- data.frame(sat_name = down_sat_names) %>%
      count(sat_name, name = "down_count")
    
    # Create complete dataframe with all satellite names for this directory
    consolidated_sat_df <- data.frame(sat_name = all_sat_names) %>%
      left_join(up_counts, by = "sat_name") %>%
      left_join(down_counts, by = "sat_name") %>%
      replace_na(list(up_count = 0, down_count = 0)) %>%
      pivot_longer(cols = c(up_count, down_count), 
                   names_to = "regulation", 
                   values_to = "count") %>%
      mutate(regulation = factor(regulation, levels = c("up_count", "down_count"),
                                  labels = c("Up-regulated", "Down-regulated")),
             comparison = dir)
    
    # Add to master dataframe
    master_sat_df <- bind_rows(master_sat_df, consolidated_sat_df)
  }
  
  # Create consolidated dataframe for chromosomes
  all_sat_chrs <- sat_chrs
  
  if(length(all_sat_chrs) > 0){
    # Count up-regulated chromosomes
    up_chr_counts <- data.frame(sat_chr = up_sat_chrs) %>%
      count(sat_chr, name = "up_count")
    
    # Count down-regulated chromosomes
    down_chr_counts <- data.frame(sat_chr = down_sat_chrs) %>%
      count(sat_chr, name = "down_count")
    
    # Create complete dataframe with all chromosomes for this directory
    consolidated_chr_df <- data.frame(sat_chr = c("chrX", "chrY", paste0("chr", 1:22))) %>%
      left_join(up_chr_counts, by = "sat_chr") %>%
      left_join(down_chr_counts, by = "sat_chr") %>%
      replace_na(list(up_count = 0, down_count = 0)) %>%
      pivot_longer(cols = c(up_count, down_count), 
                   names_to = "regulation", 
                   values_to = "count") %>%
      mutate(regulation = factor(regulation, levels = c("up_count", "down_count"),
                                  labels = c("Up-regulated", "Down-regulated")),
             sat_chr = factor(sat_chr, levels = c(paste0("chr", 1:22), "chrX", "chrY")),
             comparison = dir)
    
    # Add to master dataframe
    master_chr_df <- bind_rows(master_chr_df, consolidated_chr_df)
  }

  joined_sat_data <- right_join(sat_data, centromere_data, by = c("transcript_id"="name")) %>%
    mutate(logP = -log10(pvalue)*sign(stat))
  joined_sat_data <- left_join(joined_sat_data, em_data) %>%
    mutate(logPvalue = -log10(empvalue)*sign(meth.diff)) %>%
    mutate(chr = factor(chr, levels = c(paste0("chr", 1:22), "chrX", "chrY")))

  joined_sat_data <- joined_sat_data[order(joined_sat_data$chr, joined_sat_data$start), ]

  # Set name as factor to preserve order as they appear in sorted data
  joined_sat_data$transcript_id <- factor(joined_sat_data$transcript_id, levels = unique(joined_sat_data$transcript_id))

  # Determine stat cutoff for padj < 0.05
  stat_cutoff <- -log10(0.05)

  
  # Plot DESeq2 results
  p_expr <- joined_sat_data %>%
    ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = logP, fill = sat_name)) +
    geom_rect(alpha = 0.8) +
    scale_fill_manual(values = sat_colors) +
    facet_wrap(~chr, scales = "free_x", ncol = 4) +
    coord_cartesian(ylim = c(-4, 4)) +
    labs(title = "Centromere Satellite Expression by Chromosome",
        x = "Position on chromosome (bp)",
        y = "-log10(pvalue) * sign(stat)",
        fill = "Satellite Name") +
    theme_classic() +
    theme(strip.text = element_text(size = 9, face = "bold"),
          legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 5))
  
  # Add horizontal reference lines for significance cutoff
  if(!is.na(stat_cutoff)) {
    p_expr <- p_expr + 
      geom_hline(yintercept = stat_cutoff, linetype = "dashed", color = "grey", linewidth = 0.5) +
      geom_hline(yintercept = -stat_cutoff, linetype = "dashed", color = "grey", linewidth = 0.5)
  }

  print(p_expr)
  ggsave(file.path(wd, dir, paste0(dir, "_satellite_expression.png")), p_expr, width = 16, height = 14)
  ggsave(file.path(wd, dir, paste0(dir, "_satellite_expression.pdf")), p_expr, width = 16, height = 14)

  # Plot EMseq results
  p_em <- joined_sat_data %>%
    ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = logPvalue, fill = sat_name)) +
    geom_rect(alpha = 0.8) +
    scale_fill_manual(values = sat_colors) +
    facet_wrap(~chr, scales = "free_x", ncol = 4) +
    coord_cartesian(ylim = c(-4, 4)) +
    labs(title = "Centromere Satellite Methylation by Chromosome",
        x = "Position on chromosome (bp)",
        y = "-log10(pvalue) * sign(meth.diff)",
        fill = "Satellite Name") +
    theme_classic() +
    theme(strip.text = element_text(size = 9, face = "bold"),
          legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 5))
  
  # Add horizontal reference lines for significance cutoff
  if(!is.na(stat_cutoff)) {
    p_em <- p_em + 
      geom_hline(yintercept = stat_cutoff, linetype = "dashed", color = "grey", linewidth = 0.5) +
      geom_hline(yintercept = -stat_cutoff, linetype = "dashed", color = "grey", linewidth = 0.5)
  }

  print(p_em)
  ggsave(file.path(wd, dir, paste0(dir, "_satellite_methylation.png")), p_em, width = 16, height = 14)
  ggsave(file.path(wd, dir, paste0(dir, "_satellite_methylation.pdf")), p_em, width = 16, height = 14)

  # save joined_sat_data
  write_tsv(joined_sat_data, file.path(wd, dir, paste0(dir, "_satellite_data.tsv")))

  # plot correlation between expression and methylation
  cor_data <- joined_sat_data %>% 
    filter(!is.na(log2FoldChange) & !is.na(meth.diff) & is.finite(log2FoldChange) & is.finite(meth.diff)) %>%
    mutate(expression_status = case_when(
      log2FoldChange > 1 ~ "upregulated",
      log2FoldChange < -1 ~ "downregulated",
      TRUE ~ "noChange"
    )) %>%
    mutate(expression_status = factor(expression_status, levels = c("upregulated","downregulated", "noChange"))) %>%
    mutate(methylation_status = ifelse(meth.diff > 0, "hypermethylated", "hypomethylated")) %>%
    mutate(methylation_status = factor(methylation_status, levels = c("hypomethylated", "hypermethylated")) )

  # Test association between expression status and methylation status
  fisher_or <- NA
  fisher_p <- NA
  fisher_label <- ""
  
  if(nrow(cor_data) > 0) {
    cat(sprintf("\n=== Association Test for %s ===\n", dir))
    
    # Create contingency table
    cont_table <- table(cor_data$expression_status, cor_data$methylation_status)
    print(cont_table)
    
    # Perform chi-square test
    chi_test <- chisq.test(cont_table)
    cat(sprintf("Chi-square test: X² = %.3f, df = %d, p-value = %.2e\n", 
                chi_test$statistic, chi_test$parameter, chi_test$p.value))
    
    # Also perform Fisher's exact test (more accurate for small samples)
    fisher_test <- fisher.test(cont_table)
    fisher_p <- fisher_test$p.value
    
    # For 2x2 tables, report OR; for larger tables, just report p-value
    if(!is.null(fisher_test$estimate)) {
      fisher_or <- fisher_test$estimate
      fisher_label <- sprintf("Fisher OR = %.3f\np = %.2e", fisher_or, fisher_p)
      cat(sprintf("Fisher's exact test: OR = %.3f, p-value = %.2e\n\n", 
                  fisher_or, fisher_p))
    } else {
      fisher_label <- sprintf("Fisher p = %.2e", fisher_p)
      cat(sprintf("Fisher's exact test: p-value = %.2e\n\n", fisher_p))
    }
  }

  if(nrow(cor_data) > 0) {
    # Filter data to match the plot limits for aligned marginal histograms
    cor_data_filtered <- cor_data %>%
      filter(meth.diff >= -10 & meth.diff <= 10 & log2FoldChange >= -5 & log2FoldChange <= 5)
    
    cor_test <- cor.test(cor_data_filtered$meth.diff, cor_data_filtered$log2FoldChange, method = "pearson")
    cor_r <- cor_test$estimate
    cor_p <- cor_test$p.value
    
    p_cor <- cor_data_filtered %>%
      ggplot(aes(y = log2FoldChange, x = meth.diff)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      geom_point(aes(color = sat_name), alpha = 0.8, size = 2) +
      geom_smooth(method = "lm", color = "black", linetype = "solid", linewidth = 1, se = TRUE) +
      scale_color_manual(values = sat_colors) +
      scale_x_continuous(limits = c(-10, 10)) +
      scale_y_continuous(limits = c(-5, 5)) +
      annotate("text", x = 10, y = 5, 
               label = sprintf("r = %.3f\np = %.2e\nn = %d", cor_r, cor_p, nrow(cor_data_filtered)),
               hjust = 1.1, vjust = 1.1, size = 5, fontface = "bold") +
      annotate("text", x = 10, y = -5,
               label = fisher_label,
               hjust = 1.1, vjust = -0.1, size = 5, fontface = "bold", color = "blue") +
      labs(title = sprintf("Correlation: Expression vs Methylation (%s)", dir),
           x = "DNA Methylation (meth.diff)",
           y = "RNA Expression (log2FoldChange)",
           color = "Satellite Name") +
      theme_classic() +
      theme(legend.position = "right")
    
    # Add marginal histograms
    p_cor_with_hist <- ggExtra::ggMarginal(p_cor, type = "histogram", fill = "lightgray", color = "black", size = 4)
    
    print(p_cor_with_hist)
    ggsave(file.path(wd, dir, paste0(dir, "_correlation_expr_meth.png")), p_cor_with_hist, width = 10, height = 8)
    ggsave(file.path(wd, dir, paste0(dir, "_correlation_expr_meth.pdf")), p_cor_with_hist, width = 10, height = 8)
  }
}

# Expand master_chr_sat dataframes to include all chromosomes with 0 counts
all_chrs <- c(paste0("chr", 1:22), "chrX", "chrY")

if(nrow(master_chr_sat_up_df) > 0){
  # Get all unique combinations of comparison and sat_name
  all_combinations_up <- expand.grid(
    comparison = dirs,
    sat_name = unique(master_chr_sat_up_df$sat_name),
    sat_chr = all_chrs,
    stringsAsFactors = FALSE
  )
  
  # Merge with actual data and fill missing with 0
  master_chr_sat_up_df <- all_combinations_up %>%
    left_join(master_chr_sat_up_df, by = c("comparison", "sat_name", "sat_chr")) %>%
    replace_na(list(count = 0)) %>%
    mutate(sat_chr = factor(sat_chr, levels = all_chrs))
}

if(nrow(master_chr_sat_down_df) > 0){
  # Get all unique combinations of comparison and sat_name
  all_combinations_down <- expand.grid(
    comparison = dirs,
    sat_name = unique(master_chr_sat_down_df$sat_name),
    sat_chr = all_chrs,
    stringsAsFactors = FALSE
  )
  
  # Merge with actual data and fill missing with 0
  master_chr_sat_down_df <- all_combinations_down %>%
    left_join(master_chr_sat_down_df, by = c("comparison", "sat_name", "sat_chr")) %>%
    replace_na(list(count = 0)) %>%
    mutate(sat_chr = factor(sat_chr, levels = all_chrs))
}

# Create faceted plots from master dataframes
if(nrow(master_sat_df) > 0){
  # Calculate ordering based on total counts across all comparisons
  sat_order_master <- master_sat_df %>%
    group_by(sat_name) %>%
    summarise(total = sum(count)) %>%
    arrange(desc(total)) %>%
    pull(sat_name)
  
  master_sat_df$sat_name <- factor(master_sat_df$sat_name, levels = sat_order_master)
  master_sat_df$comparison <- factor(master_sat_df$comparison, levels = dirs)
  
  # Create faceted plot for satellite names
  p_master_sat <- ggplot(master_sat_df, aes(x = sat_name, y = count, fill = regulation)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Down-regulated" = "steelblue", "Up-regulated" = "coral")) +
    scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, max(x) + 1))))) +
    facet_wrap(~comparison, ncol = 3) +
    labs(title = "Satellite Names - Up vs Down Regulated Across All Comparisons",
         x = "Satellite Name", 
         y = "Count",
         fill = "Regulation") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          strip.text = element_text(size = 10, face = "bold"))
  
  ggsave(file.path(wd, "master_consolidated_sat_names.png"), 
         p_master_sat, width = 16, height = 12)
  ggsave(file.path(wd, "master_consolidated_sat_names.pdf"), 
         p_master_sat, width = 16, height = 12)
}

if(nrow(master_chr_df) > 0){
  master_chr_df$comparison <- factor(master_chr_df$comparison, levels = dirs)
  
  # Create faceted plot for chromosomes
  p_master_chr <- ggplot(master_chr_df, aes(x = sat_chr, y = count, fill = regulation)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Up-regulated" = "coral", "Down-regulated" = "steelblue")) +
    scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, max(x) + 1))))) +
    facet_wrap(~comparison, ncol = 3) +
    labs(title = "Chromosome Distribution - Up vs Down Regulated Across All Comparisons",
         x = "Chromosome", 
         y = "Count",
         fill = "Regulation") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          strip.text = element_text(size = 10, face = "bold"))
  
  ggsave(file.path(wd, "master_consolidated_sat_chrs.png"), 
         p_master_chr, width = 16, height = 10)
  ggsave(file.path(wd, "master_consolidated_sat_chrs.pdf"), 
         p_master_chr, width = 16, height = 10)
}

# Create stacked bar plots for chromosome-satellite combinations
# Create consistent color palette for all satellite names
all_sat_names_combined <- unique(c(
  if(nrow(master_chr_sat_up_df) > 0) unique(master_chr_sat_up_df$sat_name) else character(0),
  if(nrow(master_chr_sat_down_df) > 0) unique(master_chr_sat_down_df$sat_name) else character(0)
))
sat_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(all_sat_names_combined)),
  all_sat_names_combined
)

if(nrow(master_chr_sat_up_df) > 0){
  master_chr_sat_up_df$comparison <- factor(master_chr_sat_up_df$comparison, levels = dirs)
  
  # Create stacked plot for up-regulated genes
  p_chr_sat_up <- ggplot(master_chr_sat_up_df, aes(x = sat_chr, y = count, fill = sat_name)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = sat_colors) +
    scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, max(x) + 1))))) +
    facet_wrap(~comparison, ncol = 3) +
    labs(title = "Up-regulated: Satellite Name Composition by Chromosome Across All Comparisons",
         x = "Chromosome", 
         y = "Count",
         fill = "Satellite Name") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          strip.text = element_text(size = 10, face = "bold")) +
    guides(fill = guide_legend(ncol = 5))
  
  ggsave(file.path(wd, "master_chr_sat_up_stacked.png"), 
         p_chr_sat_up, width = 16, height = 10)
  ggsave(file.path(wd, "master_chr_sat_up_stacked.pdf"), 
         p_chr_sat_up, width = 16, height = 10)
}

if(nrow(master_chr_sat_down_df) > 0){
  master_chr_sat_down_df$comparison <- factor(master_chr_sat_down_df$comparison, levels = dirs)
  
  # Create stacked plot for down-regulated genes
  p_chr_sat_down <- ggplot(master_chr_sat_down_df, aes(x = sat_chr, y = count, fill = sat_name)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = sat_colors) +
    scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, max(x) + 1))))) +
    facet_wrap(~comparison, ncol = 3) +
    labs(title = "Down-regulated: Satellite Name Composition by Chromosome Across All Comparisons",
         x = "Chromosome", 
         y = "Count",
         fill = "Satellite Name") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          strip.text = element_text(size = 10, face = "bold")) +
    guides(fill = guide_legend(ncol = 5))
  
  ggsave(file.path(wd, "master_chr_sat_down_stacked.png"), 
         p_chr_sat_down, width = 16, height = 10)
  ggsave(file.path(wd, "master_chr_sat_down_stacked.pdf"), 
         p_chr_sat_down, width = 16, height = 10)
}

# Save master dataframes
write_csv(master_sat_df, file.path(wd, "master_sat_names_data.csv"))
write_csv(master_chr_df, file.path(wd, "master_chr_data.csv"))
write_csv(master_chr_sat_up_df, file.path(wd, "master_chr_sat_up_data.csv"))
write_csv(master_chr_sat_down_df, file.path(wd, "master_chr_sat_down_data.csv"))
