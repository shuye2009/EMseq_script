
library(data.table)
library(ggplot2)
library(dplyr)

# ==============================================================================
# Settings and Paths
# ==============================================================================

# Define analysis types and their base directories
analyses <- list(
  "Enhancer" = "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/Enhancer_targeted",
  "Promoter" = "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/Promoter_targeted",
  "Centromere" = "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/Centromere_targeted"
)

# Thresholds for coloring
pval_cutoff <- 0.05
meth_diff_cutoff <- 10  # 10% methylation difference

# ==============================================================================
# Processing Loop
# ==============================================================================

for (analysis_name in names(analyses)) {
  base_dir <- analyses[[analysis_name]]
  output_dir <- file.path(base_dir, "Volcano_Plots")
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  message(paste0("\nProcessing ", analysis_name, " analysis..."))
  
  # Get comparisons
  comparisons <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
  comparisons <- comparisons[grepl("_vs_", comparisons)]
  
  for (comp in comparisons) {
    message("  Comparison: ", comp)
    
    # Path to Targeted folder
    targeted_dir <- file.path(base_dir, comp, "Targeted")
    
    if (dir.exists(targeted_dir)) {
      # Find the 'all_diff' file
      all_diff_files <- list.files(targeted_dir, pattern = "^all_diff_.*\\.tab$", full.names = TRUE)
      
      if (length(all_diff_files) > 0) {
        # Use the first match if multiple (usually only one)
        diff_file <- all_diff_files[1]
        
        # Read Data
        data <- fread(diff_file)
        
        # Check columns
        req_cols <- c("meth.diff", "pvalue")
        if (!all(req_cols %in% colnames(data))) {
          warning("    Missing required columns in ", basename(diff_file))
          next
        }
        
        # Prepare Plot Data
        # Handle p=0 for log transformation
        min_nonzero_p <- min(data$pvalue[data$pvalue > 0], na.rm = TRUE)
        if (is.infinite(min_nonzero_p)) min_nonzero_p <- 1e-300 # Fallback if all are 0 or NA
        
        plot_data <- data %>%
          mutate(
            plot_pval = ifelse(pvalue == 0, min_nonzero_p * 0.1, pvalue),
            logP = -log10(plot_pval),
            Direction = case_when(
              pvalue < pval_cutoff & meth.diff > meth_diff_cutoff ~ "Hypermethylated",
              pvalue < pval_cutoff & meth.diff < -meth_diff_cutoff ~ "Hypomethylated",
              TRUE ~ "Not Significant"
            )
          )
        
        # Counts
        n_hyper <- sum(plot_data$Direction == "Hypermethylated")
        n_hypo  <- sum(plot_data$Direction == "Hypomethylated")
        
        # Plot
        p <- ggplot(plot_data, aes(x = meth.diff, y = logP, color = Direction)) +
          geom_point(alpha = 0.5, size = 1) +
          scale_color_manual(values = c(
            "Hypermethylated" = "#E74C3C", # Red
            "Hypomethylated"  = "#3498DB", # Blue
            "Not Significant" = "grey80"
          )) +
          labs(title = paste(analysis_name, "Volcano Plot:", comp),
               subtitle = paste0("Hyper: ", n_hyper, " | Hypo: ", n_hypo, 
                               " (p < ", pval_cutoff, ", |diff| > ", meth_diff_cutoff, "%)"),
               x = "Methylation Difference (%)",
               y = "-log10(p-value)") +
          theme_bw() +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "top"
          ) +
          geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "black", alpha = 0.3) +
          geom_vline(xintercept = c(-meth_diff_cutoff, meth_diff_cutoff), linetype = "dashed", color = "black", alpha = 0.3)
        
        # Save
        out_prefix <- file.path(output_dir, paste0("Volcano_methDiff_pval_", comp))
        ggsave(paste0(out_prefix, ".pdf"), p, width = 8, height = 6)
        ggsave(paste0(out_prefix, ".png"), p, width = 8, height = 6, dpi = 300)
        
        message("    Saved plots.")
        
      } else {
        warning("    No 'all_diff_*.tab' file found in ", targeted_dir)
      }
    }
  }
}

message("\nAll Done.")
