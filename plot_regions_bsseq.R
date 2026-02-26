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

# Define input and output paths
bismark_file <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/Enhancer_targeted/IR10Gy6d_vs_NIR/RData/bismark.RData"
hyper_regions_file <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/Enhancer_methylation_summary/IR10Gy6d_vs_NIR_Hypermethylated_Enhancers_Union.bed"
hypo_regions_file <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/Enhancer_methylation_summary/IR10Gy6d_vs_NIR_Hypomethylated_Enhancers_Union.bed"
output_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/Enhancer_methylation_summary"
comparison_tag <- "IR10Gy6d_vs_NIR"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Check if input files exist
if (!file.exists(bismark_file)) {
  stop("Bismark RData file not found: ", bismark_file)
}

if (!file.exists(hyper_regions_file)) {
  stop("Hypermethylated regions BED file not found: ", hyper_regions_file)
}

if (!file.exists(hypo_regions_file)) {
  stop("Hypomethylated regions BED file not found: ", hypo_regions_file)
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

# Define sample groups (adjust as needed if sample order changes)
sample_groups <- c("IR10Gy6d", "IR10Gy6d", "IR10Gy6d", "NIR", "NIR", "NIR")
if (length(sample_groups) != ncol(bsseq_obj)) {
  stop("Length of sample_groups (", length(sample_groups), ") does not match number of samples (", 
       ncol(bsseq_obj), "). Please update the sample_groups vector.")
}
names(sample_groups) <- colnames(bsseq_obj)

# Helper to format p-values consistently
format_p_value <- function(p_value) {
  if (is.na(p_value)) return("NA")
  if (p_value < 1e-4) return("<1e-4")
  if (p_value < 0.001) return("<0.001")
  sprintf("%.3f", p_value)
}

# Helper to retrieve region annotation (BED column 4) with fallback to name
get_region_label <- function(region, fallback_name) {
  label <- fallback_name
  meta <- mcols(region)
  if (!is.null(meta)) {
    if ("name" %in% names(meta)) {
      candidate <- as.character(meta$name)
      if (!is.null(candidate) && length(candidate) > 0 && candidate[1] != "") {
        label <- candidate[1]
      }
    } else if (ncol(meta) >= 1) {
      candidate <- as.character(meta[[1]])
      if (!is.null(candidate) && length(candidate) > 0 && candidate[1] != "") {
        label <- candidate[1]
      }
    }
  }
  label
}

# Compute methylKit-style p-values (logistic z-test and F-test) using methylated/unmethylated counts
compute_methylkit_pvalues <- function(region_data, sample_groups) {
  tryCatch({
    # Extract methylated and coverage counts from bsseq object
    m_counts <- getCoverage(region_data, type = "M")
    cov_counts <- getCoverage(region_data, type = "Cov")
    if (is.null(m_counts) || is.null(cov_counts)) {
      return(list(logistic = NA_real_, f = NA_real_))
    }
    u_counts <- cov_counts - m_counts
    # Flatten matrices into long data frame (row per CpG per sample)
    df <- data.frame(
      m = as.vector(m_counts),
      u = as.vector(u_counts),
      coverage = as.vector(cov_counts),
      group = rep(sample_groups, each = nrow(m_counts)),
      stringsAsFactors = FALSE
    )
    # Keep entries with valid coverage
    df <- df[!is.na(df$m) & !is.na(df$u) & df$coverage > 0, ]
    if (nrow(df) < 6) {
      return(list(logistic = NA_real_, f = NA_real_))
    }
    df$group <- factor(df$group, levels = unique(sample_groups))
    if (length(levels(df$group)) < 2) {
      return(list(logistic = NA_real_, f = NA_real_))
    }
    # Fit logistic regression mimicking methylKit's calculateDiffMeth
    model <- glm(cbind(m, u) ~ group, family = binomial(), data = df)
    coef_summary <- summary(model)$coefficients
    target_level <- levels(df$group)[2]
    term_name <- paste0("group", target_level)
    if (!term_name %in% rownames(coef_summary)) {
      return(list(logistic = NA_real_, f = NA_real_))
    }
    logistic_p <- coef_summary[term_name, "Pr(>|z|)"]
    z_value <- coef_summary[term_name, "Estimate"] / coef_summary[term_name, "Std. Error"]
    f_stat <- z_value^2
    df_resid <- max(1, model$df.residual)
    f_p <- pf(f_stat, 1, df_resid, lower.tail = FALSE)
    list(logistic = logistic_p, f = f_p)
  }, error = function(e) {
    list(logistic = NA_real_, f = NA_real_)
  })
}

# Load hypermethylated regions
cat("Loading hypermethylated regions from:", hyper_regions_file, "\n")
hyper_regions <- import.bed(hyper_regions_file)
cat("Loaded", length(hyper_regions), "hypermethylated regions\n")

# Load hypomethylated regions
cat("Loading hypomethylated regions from:", hypo_regions_file, "\n")
hypo_regions <- import.bed(hypo_regions_file)
cat("Loaded", length(hypo_regions), "hypomethylated regions\n")

# Check if regions have proper names
if (is.null(names(hyper_regions))) {
  names(hyper_regions) <- paste0("Hyper_", seq_along(hyper_regions))
}

if (is.null(names(hypo_regions))) {
  names(hypo_regions) <- paste0("Hypo_", seq_along(hypo_regions))
}

cat("Creating methylation profile plots for both hyper and hypomethylated enhancers...\n")

# Function to plot regions
plot_enhancer_regions <- function(regions, region_type, bsseq_obj, output_dir, sample_groups) {
  cat("Plotting", length(regions), region_type, "enhancers...\n")
  
  # Create a comprehensive plot similar to plotManyRegions
  pdf(file.path(output_dir, paste0(comparison_tag, "_custom_", region_type, "_methylation_profiles.pdf")),
      width = 15, height = 10)
  
  # Set up layout for multiple regions
  num_regions <- length(regions)
  if (num_regions <= 6) {
    par(mfrow = c(2, 3))
  } else if (num_regions <= 9) {
    par(mfrow = c(3, 3))
  } else {
    par(mfrow = c(4, 3))
  }
  
  # Plot each region
  for (i in seq_along(regions)) {
    region <- regions[i]
    region_name <- names(regions)[i]
    region_label <- get_region_label(region, region_name)
    
    # Extract methylation data for this region with 2kb extension
    region_gr <- GRanges(seqnames = seqnames(region), 
                         ranges = IRanges(start = start(region) - 2000, 
                                        end = end(region) + 2000))
    
    # Get methylation data overlapping this region
    overlaps <- findOverlaps(region_gr, rowRanges(bsseq_obj))
    
    if (length(overlaps) > 0) {
      region_data <- bsseq_obj[subjectHits(overlaps)]
      
      # Get methylation data for all samples
      meth_matrix <- getMeth(region_data, type = "raw")
      positions <- start(region_data)
      
      # Compute methylKit-style p-values
      p_values <- compute_methylkit_pvalues(region_data, sample_groups)
      logistic_text <- format_p_value(p_values$logistic)
      f_text <- format_p_value(p_values$f)
      
      # Create plot with all samples
      plot(positions, meth_matrix[,1], type = "n", 
           ylim = c(0, 1), xlab = "Genomic Position", ylab = "Methylation (%)",
           main = paste0(region_label, "\n", length(overlaps), " sites, p_log=", logistic_text,
                        ", p_F=", f_text))
      
      # Define colors for experimental groups
      sample_colors <- ifelse(sample_groups == "NIR", "red", "blue")
      
      # Plot each sample as points with group-specific colors
      for (j in seq_len(ncol(meth_matrix))) {
        points(positions, meth_matrix[,j], col = sample_colors[j], pch = 16, cex = 0.5)
      }
      
      # Add smooth lines for each sample with better handling
      for (j in seq_len(ncol(meth_matrix))) {
        # Remove NA values and get valid data
        valid_idx <- !is.na(meth_matrix[,j])
        if (sum(valid_idx) > 3) {  # Need at least 3 points for smoothing
          valid_positions <- positions[valid_idx]
          valid_meth <- meth_matrix[valid_idx, j]
          
          # Sort by position for proper line drawing
          sort_order <- order(valid_positions)
          sorted_positions <- valid_positions[sort_order]
          sorted_meth <- valid_meth[sort_order]
          
          # Try Gaussian kernel smoothing first, fallback to loess if needed
          smooth_line <- tryCatch({
            # Use ksmooth for Gaussian kernel smoothing
            ksmooth_result <- ksmooth(sorted_positions, sorted_meth, 
                                     bandwidth = 1000,  # Adjust based on region size
                                     kernel = "normal")
            list(x = ksmooth_result$x, y = ksmooth_result$y)
          }, error = function(e) {
            # Fallback to loess if ksmooth fails
            lowess(sorted_positions, sorted_meth, f = 0.4)
          })
          
          # Draw smooth line
          lines(smooth_line$x, smooth_line$y, col = sample_colors[j], lwd = 1.5, lty = 1)
        }
      }
      
      # Add CpG position markers on x-axis
      # Get unique positions (avoid overplotting if multiple samples have same positions)
      unique_positions <- unique(positions)
      
      # Add small tick marks for CpG positions
      axis(1, at = unique_positions, labels = FALSE, tcl = -0.2, col = "gray50")
      
      # Add region boundaries
      abline(v = c(start(region), end(region)), col = "black", lty = 2, lwd = 2)
      
      # Add legend for experimental groups
      legend("topright", c("IR10Gy6d", "NIR"), col = c("blue", "red"), 
             pch = 16, cex = 0.8)
    } else {
      # No data for this region
      plot(1, 1, type = "n", xlab = "", ylab = "",
           main = paste(region_label, "- No data"))
    }
  }
  
  # Add overall title
  mtext(paste0(region_type, "Enhancers - Methylation Profiles (All Samples)"), 
        outer = TRUE, cex = 1.5, line = -1)
  
  dev.off()
  cat("Custom", region_type, "multi-sample plot saved\n")
  
  # Create individual plots for each region (more detailed)
  pdf(file.path(output_dir, paste0(comparison_tag, "_detailed_", region_type, "_region_plots.pdf")),
      width = 12, height = 8)
  
  for (i in seq_along(regions)) {
    region <- regions[i]
    region_name <- names(regions)[i]
    
    # Extract methylation data for this region with 2kb extension
    region_gr <- GRanges(seqnames = seqnames(region), 
                         ranges = IRanges(start = start(region) - 2000, 
                                        end = end(region) + 2000))
    
    # Get methylation data overlapping this region
    overlaps <- findOverlaps(region_gr, rowRanges(bsseq_obj))
    
    if (length(overlaps) > 0) {
      region_data <- bsseq_obj[subjectHits(overlaps)]
      
      # Get methylation data for all samples
      meth_matrix <- getMeth(region_data, type = "raw")
      positions <- start(region_data)
      
      # Compute methylKit-style p-values
      p_values <- compute_methylkit_pvalues(region_data, sample_groups)
      logistic_text <- format_p_value(p_values$logistic)
      f_text <- format_p_value(p_values$f)
      
      # Create detailed plot
      plot(positions, meth_matrix[,1], type = "n", 
           ylim = c(0, 1), xlab = "Genomic Position", ylab = "Methylation (%)",
           main = paste("Region:", region_label, "(", length(overlaps), " sites, p_log=",
                        logistic_text, ", p_F=", f_text, ")"))
      
      # Define colors for experimental groups
      sample_colors <- ifelse(sample_groups == "NIR", "red", "blue")
      
      # Plot each sample as points with group-specific colors
      for (j in seq_len(ncol(meth_matrix))) {
        points(positions, meth_matrix[,j], col = sample_colors[j], pch = 16, cex = 0.8)
      }
      
      # Add smooth lines for each sample with better handling
      for (j in seq_len(ncol(meth_matrix))) {
        # Remove NA values and get valid data
        valid_idx <- !is.na(meth_matrix[,j])
        if (sum(valid_idx) > 3) {  # Need at least 3 points for smoothing
          valid_positions <- positions[valid_idx]
          valid_meth <- meth_matrix[valid_idx, j]
          
          # Sort by position for proper line drawing
          sort_order <- order(valid_positions)
          sorted_positions <- valid_positions[sort_order]
          sorted_meth <- valid_meth[sort_order]
          
          # Try Gaussian kernel smoothing first, fallback to loess if needed
          smooth_line <- tryCatch({
            # Use ksmooth for Gaussian kernel smoothing
            ksmooth_result <- ksmooth(sorted_positions, sorted_meth, 
                                     bandwidth = 1000,  # Adjust based on region size
                                     kernel = "normal")
            list(x = ksmooth_result$x, y = ksmooth_result$y)
          }, error = function(e) {
            # Fallback to loess if ksmooth fails
            lowess(sorted_positions, sorted_meth, f = 0.4)
          })
          
          # Draw smooth line
          lines(smooth_line$x, smooth_line$y, col = sample_colors[j], lwd = 1.8, lty = 1)
        }
      }
      
      # Add CpG position markers on x-axis
      # Get unique positions (avoid overplotting if multiple samples have same positions)
      unique_positions <- unique(positions)
      
      # Add small tick marks for CpG positions
      axis(1, at = unique_positions, labels = FALSE, tcl = -0.3, col = "gray50")
      
      # Add region boundaries
      abline(v = c(start(region), end(region)), col = "black", lty = 2, lwd = 2)
      
      # Add legend for experimental groups
      legend("topright", c("IR10Gy6d", "NIR"), col = c("blue", "red"), 
             pch = 16, cex = 0.9)
      
      # Add grid
      grid()
    } else {
      # No data for this region
      plot(1, 1, type = "n", xlab = "", ylab = "",
           main = paste("Region:", region_label, "- No data"))
    }
  }
  
  dev.off()
  cat("Detailed", region_type, "individual plots saved\n")
  
  # Create a summary statistics plot with p-values
  pdf(file.path(output_dir, paste0(comparison_tag, "_", region_type, "_methylation_summary_stats.pdf")),
      width = 12, height = 8)

  # Calculate average methylation and p-values for each region across all samples
  region_means <- c()
  region_names <- c()
  region_labels <- rep(NA_character_, length(regions))
  region_pvalues_log <- c()
  region_pvalues_f <- c()

  for (i in seq_along(regions)) {
    region <- regions[i]
    region_name <- names(regions)[i]
    region_label <- get_region_label(region, region_name)
    region_labels[i] <- region_label
    
    # Extract methylation data for this region
    region_gr <- GRanges(seqnames = seqnames(region), 
                         ranges = IRanges(start = start(region) - 2000, 
                                        end = end(region) + 2000))
    
    overlaps <- findOverlaps(region_gr, rowRanges(bsseq_obj))
    
    if (length(overlaps) > 0) {
      region_data <- bsseq_obj[subjectHits(overlaps)]
      meth_matrix <- getMeth(region_data, type = "raw")
      region_means[i] <- mean(meth_matrix, na.rm = TRUE)
      region_names[i] <- region_name
      
      # Perform methylKit-style tests between groups
      p_vals <- compute_methylkit_pvalues(region_data, sample_groups)
      region_pvalues_log[i] <- p_vals$logistic
      region_pvalues_f[i] <- p_vals$f
    }
  }
  
  # Create barplot of average methylation with significance indicators
  if (length(region_means) > 0) {
    # Create colors based on significance
    colors <- ifelse(is.na(region_pvalues_log), "gray80",
                    ifelse(region_pvalues_log < 0.05, "red", 
                          ifelse(region_pvalues_log < 0.1, "orange", "steelblue")))
    
    barplot(region_means, names.arg = region_names, 
            ylab = "Average Methylation (%)", 
            main = paste("Average Methylation Across All", region_type, "Regions"),
            col = colors, las = 2, cex.names = 0.8)
    abline(h = 0.5, col = "red", lty = 2)
    
    # Add legend for significance
    legend("topright", c("p < 0.05", "p < 0.1", "p ≥ 0.1", "NA"), 
           fill = c("red", "orange", "steelblue", "gray80"), cex = 0.8)
  } else {
    plot(1, 1, type = "n", main = "No data available for summary")
  }
  
  dev.off()
  cat("Summary statistics plot saved for", region_type, "regions\n")
  
  # Create a summary table with p-values
  summary_table <- data.frame(
    Region = region_names,
    Region_Annotation = region_labels,
    Num_Sites = sapply(seq_along(regions), function(i) {
      region <- regions[i]
      region_gr <- GRanges(seqnames = seqnames(region), 
                           ranges = IRanges(start = start(region) - 2000, 
                                          end = end(region) + 2000))
      overlaps <- findOverlaps(region_gr, rowRanges(bsseq_obj))
      length(overlaps)
    }),
    Avg_Methylation = region_means,
    P_Value_Logistic = region_pvalues_log,
    P_Value_F = region_pvalues_f,
    stringsAsFactors = FALSE
  )
  
  # Format p-values in the table
  summary_table$Significance <- ifelse(is.na(summary_table$P_Value_Logistic), "NA",
                                       ifelse(summary_table$P_Value_Logistic < 0.001, "***",
                                              ifelse(summary_table$P_Value_Logistic < 0.01, "**",
                                                     ifelse(summary_table$P_Value_Logistic < 0.05, "*",
                                                            ifelse(summary_table$P_Value_Logistic < 0.1, ".", "ns")))))
  
  # Save summary table
  write.table(summary_table, 
              file = file.path(output_dir, paste0(comparison_tag, "_", region_type, "_regions_summary_table.tsv")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("Summary table with p-values saved for", region_type, "regions\n")
}

# Plot hypermethylated regions
plot_enhancer_regions(hyper_regions, "Hypermethylated", bsseq_obj, output_dir, sample_groups)

# Plot hypomethylated regions
plot_enhancer_regions(hypo_regions, "Hypomethylated", bsseq_obj, output_dir, sample_groups)

cat("Script completed successfully!\n")
cat("Output files saved in:", output_dir, "\n")
