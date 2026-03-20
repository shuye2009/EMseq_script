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

# Define input and output paths
genome <- "hs1"
minCpGs <- 5
cutoff <- 0.05
cores <- 1
testCovariate <- "Group"

region_bed_file <- "C:/PROJECTS/resource/T2T_CHM13/ENCFF912FUA_MCF10A_element_gene_links_thresholded_Engreitz_T2T_6col.bed"
bismark_file <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/Enhancer_targeted/IR10Gy6d_vs_NIR/RData/bismark.RData"
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

# MethylKit-style filtering/normalization parameters
min_coverage_threshold <- 1L
min_cpg_threshold <- 5L
min_per_group_threshold <- 2L

# Helper to format p-values consistently
format_p_value <- function(p_value) {
  if (is.na(p_value)) return("NA")
  if (p_value < 0.001) return("<0.001")
  sprintf("%.3f", p_value)
}

# Normalize coverage across samples similar to methylKit's normalizeCoverage
normalize_coverage_counts <- function(m_counts, cov_counts) {
  if (is.null(m_counts) || is.null(cov_counts)) {
    return(list(m = NULL, cov = NULL))
  }
  totals <- colSums(cov_counts, na.rm = TRUE)
  valid_totals <- totals > 0
  if (!any(valid_totals)) {
    return(list(m = m_counts, cov = cov_counts))
  }
  target <- median(totals[valid_totals])
  scale_factors <- totals / target
  scale_factors[!valid_totals] <- 1
  cov_norm <- sweep(cov_counts, 2, scale_factors, "/", check.margin = FALSE)
  m_norm <- sweep(m_counts, 2, scale_factors, "/", check.margin = FALSE)
  cov_norm <- round(cov_norm)
  m_norm <- round(m_norm)
  m_norm <- pmin(m_norm, cov_norm)
  list(m = m_norm, cov = cov_norm)
}

# Apply methylKit-style coverage/CpG filters
filter_counts_by_thresholds <- function(m_counts, cov_counts, sample_groups,
                                        min_coverage, min_cpg, min_per_group) {
  if (is.null(m_counts) || is.null(cov_counts)) {
    return(list(m = NULL, cov = NULL))
  }
  coverage_mask <- cov_counts >= min_coverage
  keep_rows <- rep(TRUE, nrow(cov_counts))
  for (grp in unique(sample_groups)) {
    grp_cols <- which(sample_groups == grp)
    grp_support <- rowSums(coverage_mask[, grp_cols, drop = FALSE], na.rm = TRUE)
    keep_rows <- keep_rows & (grp_support >= min_per_group)
  }
  keep_rows[is.na(keep_rows)] <- FALSE
  keep_rows <- keep_rows & (rowSums(cov_counts, na.rm = TRUE) > 0)
  if (sum(keep_rows) < min_cpg) {
    return(list(m = NULL, cov = NULL))
  }
  list(m = m_counts[keep_rows, , drop = FALSE],
       cov = cov_counts[keep_rows, , drop = FALSE])
}

# Internal helper to fit GLM given coverage matrices
fit_glm_from_counts <- function(m_counts, cov_counts, sample_groups, aggregate = FALSE) {
  if (is.null(m_counts) || is.null(cov_counts) || length(m_counts) == 0) {
    return(list(logistic = NA_real_, f = NA_real_))
  }
  if (aggregate) {
    agg_m <- colSums(m_counts, na.rm = TRUE)
    agg_cov <- colSums(cov_counts, na.rm = TRUE)
    df <- data.frame(
      m = as.numeric(agg_m),
      u = as.numeric(agg_cov - agg_m),
      coverage = as.numeric(agg_cov),
      group = sample_groups,
      stringsAsFactors = FALSE
    )
  } else {
    df <- data.frame(
      m = as.vector(m_counts),
      u = as.vector(cov_counts - m_counts),
      coverage = as.vector(cov_counts),
      group = rep(sample_groups, each = nrow(m_counts)),
      stringsAsFactors = FALSE
    )
  }
  df <- df[!is.na(df$m) & !is.na(df$u) & df$coverage > 0, , drop = FALSE]
  if (nrow(df) < 6) {
    return(list(logistic = NA_real_, f = NA_real_))
  }
  df$group <- factor(df$group, levels = levels(as.factor(sample_groups)))
  if (length(levels(df$group)) < 2) {
    return(list(logistic = NA_real_, f = NA_real_))
  }
  model <- tryCatch(glm(cbind(m, u) ~ group, family = binomial(), data = df),
                    error = function(e) NULL)
  if (is.null(model)) {
    return(list(logistic = NA_real_, f = NA_real_))
  }
  model_summary <- summary(model)
  coef_summary <- model_summary$coefficients
  target_level <- levels(df$group)[2]
  term_name <- paste0("group", target_level)
  if (!term_name %in% rownames(coef_summary)) {
    return(list(logistic = NA_real_, f = NA_real_))
  }
  dispersion <- model_summary$dispersion
  if (is.null(dispersion) || is.na(dispersion)) {
    dispersion <- 1
  }
  estimate <- coef_summary[term_name, "Estimate"]
  std_error <- coef_summary[term_name, "Std. Error"] * sqrt(dispersion)
  if (is.na(std_error) || std_error == 0) {
    return(list(logistic = NA_real_, f = NA_real_))
  }
  z_value <- estimate / std_error
  logistic_p <- 2 * pnorm(-abs(z_value))
  f_stat <- z_value^2
  df_resid <- max(1, model$df.residual)
  f_p <- pf(f_stat, 1, df_resid, lower.tail = FALSE)
  list(logistic = logistic_p, f = f_p)
}

# Compute GLM-style p-values (logistic from per-CpG counts, F from region-level counts)
compute_GLM_pvalues <- function(region_data, sample_groups) {
  group_levels <- levels(sample_groups)
  default_group_means <- setNames(rep(NA_real_, length(group_levels)), group_levels)
  raw_m_counts <- getCoverage(region_data, type = "M")
  raw_cov_counts <- getCoverage(region_data, type = "Cov")
  normalized <- normalize_coverage_counts(raw_m_counts, raw_cov_counts)
  filtered <- filter_counts_by_thresholds(
    normalized$m,
    normalized$cov,
    sample_groups,
    min_coverage_threshold,
    min_cpg_threshold,
    min_per_group_threshold
  )
  if (is.null(filtered$m)) {
    return(list(
      logistic = NA_real_,
      f = NA_real_,
      rlogistic = NA_real_,
      rf = NA_real_,
      group_means = default_group_means,
      meth_diff = NA_real_,
      avg_meth = NA_real_,
      num_sites = 0L
    ))
  }


  per_cpg <- fit_glm_from_counts(filtered$m, filtered$cov, sample_groups, aggregate = FALSE)
  per_region <- fit_glm_from_counts(filtered$m, filtered$cov, sample_groups, aggregate = TRUE)
  
  # Ensure scalar values for p-values
  cpg_log <- if(length(per_cpg$logistic) == 1) per_cpg$logistic else NA_real_
  cpg_f <- if(length(per_cpg$f) == 1) per_cpg$f else NA_real_
  reg_log <- if(length(per_region$logistic) == 1) per_region$logistic else NA_real_
  reg_f <- if(length(per_region$f) == 1) per_region$f else NA_real_
  
  group_means <- default_group_means
  for (grp in group_levels) {
    grp_cols <- which(sample_groups == grp)
    grp_m <- sum(filtered$m[, grp_cols, drop = FALSE], na.rm = TRUE)
    grp_cov <- sum(filtered$cov[, grp_cols, drop = FALSE], na.rm = TRUE)
    group_means[grp] <- if (grp_cov > 0) (grp_m / grp_cov) * 100 else NA_real_
  }
  meth_diff <- if (length(group_levels) >= 2) {
    group_means[group_levels[2]] - group_means[group_levels[1]]
  } else {
    NA_real_
  }
  total_m <- sum(filtered$m, na.rm = TRUE)
  total_cov <- sum(filtered$cov, na.rm = TRUE)
  avg_meth <- if (total_cov > 0) (total_m / total_cov) * 100 else NA_real_
  list(
    logistic = cpg_log,
    f = cpg_f,
    rlogistic = reg_log,
    rf = reg_f,
    group_means = group_means,
    meth_diff = meth_diff,
    avg_meth = avg_meth,
    num_sites = nrow(filtered$m)
  )
}

# Run targeted methylKit analysis mimicking DM.R pipeline
run_methylKit_style_region_test <- function(target_region_file, bsseq_obj, testCovariate,
                                      output_dir, comparison_tag, genome) {
  if (!file.exists(target_region_file)) {
    cat("Target region file", target_region_file, "does not exist\n")
    return(NULL)
  }

  targetRegions <- import.bed(target_region_file)
  print(head(targetRegions))
  print(length(targetRegions))
  
  print(glue::glue("Loaded {length(targetRegions)} target regions"))
  
  # Convert BSseq to methylKit format for region-based testing
  cat("Converting BSseq to methylKit format...\n")
      
  print(head(bs.filtered@assays@data$M))
  # Extract methylation data from BSseq object
  meth_data <- bsseq::getCoverage(bs.filtered, type = "M")
  print(head(meth_data))
  cov_data <- bsseq::getCoverage(bs.filtered, type = "Cov")
  print(head(cov_data))
  pos_data <- GenomicRanges::granges(bs.filtered)
  print(head(pos_data))
  
  # Get sample information
  sample_info <- bsseq::pData(bs.filtered)
  
  # Create methylKit objects for each sample
  cat("Building methylRaw list...\n")
  methyRaw_list <- list()
  for(i in 1:ncol(meth_data)) {
    sample_name <- colnames(meth_data)[i]
    
    # Create data frame for this sample
    sample_df <- data.frame(
      chr = as.character(seqnames(pos_data)),
      start = start(pos_data),
      end = end(pos_data),
      strand = strand(pos_data),
      coverage = cov_data[,i],
      numC = meth_data[,i],
      numT = cov_data[,i] - meth_data[,i]
    )
    
    # Remove rows with NA values
    #sample_df <- sample_df[complete.cases(sample_df),]
    print(dim(sample_df))

    # Convert dataframe to methylRaw
    
    methylraw_obj <- new("methylRaw",
                  sample_df,
                  sample.id = sample_name,
                  assembly = genome,
                  context = "CpG",
                  resolution = "base")
    
    methyRaw_list[[i]] <- methylraw_obj
  }
      
  # Create treatment vector based on your design
  cat("Creating methylRawList object...\n")
  treatment_vector <- as.numeric(as.factor(sample_info[[testCovariate]])) - 1
  
  # Create methylRawList
  myMethylListobj <- new("methylRawList", methyRaw_list, treatment = treatment_vector)
  print(head(myMethylListobj))

  # Unite samples (keep sites covered in at least 2 samples per group)
  cat("Uniting methylRawList objects...\n")
  meth_united <- methylKit::unite(myMethylListobj, destrand = FALSE, min.per.group = 2L)
  print(head(meth_united))
  # Calculate regional methylation for target regions
  cat("Calculating regional methylation for target regions...\n")
  
  # Get regional methylation
  regional_meth <- methylKit::regionCounts(meth_united, targetRegions, 
                                            cov.bases = minCpGs, 
                                            strand.aware = FALSE)
  print(head(regional_meth))
  # Perform differential methylation testing on regions
  cat("Performing differential methylation testing on target regions...\n")
      
  # Calculate differential methylation for regions
  target_diff <- methylKit::calculateDiffMeth(regional_meth, 
                                              overdispersion = "MN",
                                              test = "F",
                                              mc.cores = cores)
  print(head(target_diff))
  target_diff_df <- getData(target_diff)

  if(nrow(target_diff_df) > 0) {
    
    target_diff_df <- merge(target_diff_df, targetBed, by = c("chr", "start", "end", "strand"), all.x = TRUE)
    # Export results
    print(glue::glue("Exporting targeted region results..."))
    write.table(target_diff_df, file.path(output_dir, paste0("all_diff_", comparison_tag, ".tab")), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
    # Filter significant results
    cat("Filtering significant results...\n")
    target_diff_sig <- methylKit::getMethylDiff(target_diff, 
                                                difference = cutoff * 100,  # Convert to percentage
                                                qvalue = 0.05)
    print(head(target_diff_sig))
    target_diff_sig_df <- getData(target_diff_sig)
    if(nrow(target_diff_sig_df) > 0) {
      
      target_diff_sig_df <- merge(target_diff_sig_df, targetBed, by = c("chr", "start", "end", "strand"), all.x = TRUE)

      print(glue::glue("Found {nrow(target_diff_sig_df)} significant targeted regions"))
      write.table(target_diff_sig_df, file.path(output_dir, paste0("significant_diff_", comparison_tag, ".tab")), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

    
    }else{
      print(glue::glue("No significant differentially methylated targeted regions found"))
    }
  }
  return(target_diff_df)
}



# Function to plot regions
run_cpg_level_region_test <- function(region_bed_file, bsseq_obj, testCovariate, output_dir, methylkit_results = NULL) {

  regions <- import.bed(region_bed_file)
  sample_info <- colData(bsseq_obj)
  print(sample_info)

  # Define sample groups (adjust as needed if sample order changes)
  sample_groups <- sample_info[[testCovariate]]
  if (length(sample_groups) != ncol(bsseq_obj)) {
    stop("Length of sample_groups (", length(sample_groups), ") does not match number of samples (", 
        ncol(bsseq_obj), "). Please update the sample_groups vector.")
  }
  names(sample_groups) <- rownames(sample_info)
  
  # Parallelize stats computation
  ncores <- min(cores, detectCores() - 1, length(regions))
  if (ncores < 1) ncores <- 1
  
  cat("Computing statistics using", ncores, "cores...\n")
  
  # Function to process each region
  process_region <- function(i) {
    region <- regions[i]
    region_label <- regions$name[i]
    
    overlaps_stats <- findOverlaps(region, rowRanges(bsseq_obj))
    
    if (length(overlaps_stats) > 0) {
      region_data_stats <- bsseq_obj[subjectHits(overlaps_stats)]
      p_values <- compute_GLM_pvalues(region_data_stats, sample_groups)
      return(list(label = region_label, values = p_values))
    } else {
      return(NULL)
    }
  }
  
  if (ncores > 1) {
    # Use different parallelization strategy based on OS
    if (.Platform$OS.type == "unix") {
      # Unix/Linux/Mac: use mclapply (forking)
      cat("Using forking-based parallelization (mclapply)...\n")
      results <- mclapply(seq_along(regions), process_region, 
                          mc.cores = ncores)
    } else {
      # Windows: use socket cluster
      cat("Using socket-based parallelization (parLapply)...\n")
      cl <- makeCluster(ncores)
      on.exit(stopCluster(cl), add = TRUE)
      
      # Export necessary objects to cluster
      clusterExport(cl, c("bsseq_obj", "sample_groups", "regions",
                          "compute_GLM_pvalues", 
                          "normalize_coverage_counts",
                          "filter_counts_by_thresholds", 
                          "fit_glm_from_counts",
                          "min_coverage_threshold", "min_cpg_threshold",
                          "min_per_group_threshold"),
                    envir = environment())
      
      # Load required libraries on each worker
      clusterEvalQ(cl, {
        library(bsseq)
        library(GenomicRanges)
      })
      
      results <- parLapply(cl, seq_along(regions), process_region)
    }
    
    # Convert list to named list, removing NULL entries
    stats_list <- setNames(
      lapply(results[!sapply(results, is.null)], `[[`, "values"),
      sapply(results[!sapply(results, is.null)], `[[`, "label")
    )
  } else {
    stats_list <- list()
    for (i in seq_along(regions)) {
      region <- regions[i]
      region_label <- regions$name[i]
      
      overlaps_stats <- findOverlaps(region, rowRanges(bsseq_obj))
      
      if (length(overlaps_stats) > 0) {
        region_data_stats <- bsseq_obj[subjectHits(overlaps_stats)]
        p_values <- compute_GLM_pvalues(region_data_stats, sample_groups)
        stats_list[[region_label]] <- p_values
      }else{
        stats_list[[region_label]] <- NULL
      }
    }
  }

  # Create a summary statistics with multiple p-values

  Num_NAs <- sum(sapply(stats_list, is.null))
  num_regions <- length(regions)
  group_levels <- levels(sample_groups)
  group_mean_cols <- paste0("Mean_", group_levels)
  
  group_means_mat <- matrix(NA_real_, nrow = num_regions, ncol = length(group_levels))
  colnames(group_means_mat) <- group_mean_cols
  region_labels <- regions$name

  stats_list <- stats_list[region_labels]
  
  default_group_means <- setNames(rep(NA_real_, length(group_levels)), 
                                  group_levels)
  group_means_mat <- lapply(stats_list, function(x) {
    if(is.null(x)) default_group_means else x$group_means
  })
  group_means_mat <- do.call(rbind, group_means_mat)
  rownames(group_means_mat) <- NULL
  
  summary_table <- data.frame(
    chr = as.character(seqnames(regions)),
    start = as.integer(start(regions)),
    end = as.integer(end(regions)),
    strand = as.character(strand(regions)),
    Region_Annotation = region_labels,
    Num_Sites = sapply(stats_list, function(x) if(is.null(x)) 0L else x$num_sites),
    Avg_Methylation = sapply(stats_list, function(x) if(is.null(x)) NA_real_ else x$avg_meth),
    Meth_Diff = sapply(stats_list, function(x) if(is.null(x)) NA_real_ else x$meth_diff),
    P_Value_Logistic = sapply(stats_list, function(x) if(is.null(x)) NA_real_ else x$rlogistic),
    P_Value_F = sapply(stats_list, function(x) if(is.null(x)) NA_real_ else x$rf),
    cpg_P_Value_Logistic = sapply(stats_list, function(x) if(is.null(x)) NA_real_ else x$logistic),
    cpg_P_Value_F = sapply(stats_list, function(x) if(is.null(x)) NA_real_ else x$f),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  summary_table <- cbind(summary_table, group_means_mat)
  
  # Remove rows with no CpG sites (regions without overlaps)
  summary_table <- summary_table[summary_table$Num_Sites > 0, ]

  if (!is.null(methylkit_results)) {
    methylkit_results_subset <- methylkit_results[, c("chr", "start", "end", "strand", "pvalue", "qvalue", "meth.diff", "name")]
    colnames(methylkit_results_subset) <- c("chr", "start", "end", "strand", "MethylKit_pvalue", "MethylKit_qvalue", "MethylKit_meth_diff", "name")
    
    summary_table <- merge(summary_table, methylkit_results_subset,
                           by.x = c("chr", "start", "end", "strand"),
                           by.y = c("chr", "start", "end", "strand"),
                           all.x = TRUE, sort = FALSE)
  }
  
  # Format p-values in the table
  summary_table$cpg_qvalue_F <- p.adjust(as.numeric(summary_table$cpg_P_Value_F), method = "fdr")
  summary_table$Significance <- ifelse(is.na(summary_table$cpg_qvalue_F), "NA",
                                       ifelse(summary_table$cpg_qvalue_F < 0.001, "***",
                                              ifelse(summary_table$cpg_qvalue_F < 0.01, "**",
                                                     ifelse(summary_table$cpg_qvalue_F < 0.05, "*",
                                                            ifelse(summary_table$cpg_qvalue_F < 0.1, ".", "ns")))))
  sig_regions_table <- summary_table[summary_table$cpg_qvalue_F < 0.05,]

  
  # Save summary table
  write.table(summary_table, 
              file = file.path(output_dir, paste0(comparison_tag, "_regions_summary_table.tsv")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(sig_regions_table, 
              file = file.path(output_dir, paste0(comparison_tag, "_regions_sig_regions_table.tsv")), 
              sep = "\t", quote = FALSE, row.names = FALSE)  

  cat("Summary table with p-values saved for", region_type, "regions\n")
  
  sig_regions <- regions[regions$name %in% sig_regions_table$Region_Annotation,]

  pdf(file.path(output_dir, paste0(comparison_tag, "_sig_regions.pdf")))
  par(mfrow = c(2, 2))
  for (i in seq_along(sig_regions)) {
    region <- sig_regions[i]
    
    region_label <- sig_regions$name[i]
    
    region_gr_plot <- GRanges(seqnames = seqnames(region), 
                              ranges = IRanges(start = start(region) - 2000, 
                                               end = end(region) + 2000))
    region_gr_stats <- GRanges(seqnames = seqnames(region), 
                               ranges = IRanges(start = start(region), 
                                                end = end(region)))
    
    overlaps_plot <- findOverlaps(region_gr_plot, rowRanges(bsseq_obj))
    overlaps_stats <- findOverlaps(region_gr_stats, rowRanges(bsseq_obj))
    
    if (length(overlaps_plot) > 0 && length(overlaps_stats) > 0) {
      region_data_plot <- bsseq_obj[subjectHits(overlaps_plot)]
      meth_matrix <- getMeth(region_data_plot, type = "raw")
      positions <- start(region_data_plot)
      
      region_data_stats <- bsseq_obj[subjectHits(overlaps_stats)]
      p_values <- stats_list[[region_label]][["p_"]]
      stats_list[[region_label]] <- p_values
      f_text <- format_p_value(p_values$f)
      
      plot(positions, meth_matrix[,1], type = "n", 
           ylim = c(0, 1), xlab = "Genomic Position", ylab = "Methylation (%)",
           main = paste("Region:", region_label, "(", p_values$num_sites, " sites, p_F=", f_text, ")"))
      
      sample_colors <- ifelse(sample_groups == "cNIR", "blue", "red")
      for (j in seq_len(ncol(meth_matrix))) {
        points(positions, meth_matrix[,j], col = sample_colors[j], pch = 16, cex = 0.5)
      }
      for (j in seq_len(ncol(meth_matrix))) {
        valid_idx <- !is.na(meth_matrix[,j])
        if (sum(valid_idx) > 3) {
          valid_positions <- positions[valid_idx]
          valid_meth <- meth_matrix[valid_idx, j]
          sort_order <- order(valid_positions)
          sorted_positions <- valid_positions[sort_order]
          sorted_meth <- valid_meth[sort_order]
          smooth_line <- tryCatch({
            ksmooth_result <- ksmooth(sorted_positions, sorted_meth, 
                                      bandwidth = 1000, kernel = "normal")
            list(x = ksmooth_result$x, y = ksmooth_result$y)
          }, error = function(e) {
            lowess(sorted_positions, sorted_meth, f = 0.4)
          })
          lines(smooth_line$x, smooth_line$y, col = sample_colors[j], lwd = 1.5, lty = 1)
        }
      }
      unique_positions <- unique(positions)
      axis(1, at = unique_positions, labels = FALSE, tcl = -0.2, col = "gray50")
      abline(v = c(start(region), end(region)), col = "black", lty = 2, lwd = 2)
      legend("topright", levels(sample_groups), col = c("blue", "red"), 
             pch = 16, cex = 0.8)
    } else {
      plot(1, 1, type = "n", xlab = "", ylab = "",
           main = paste(region_label, "- No data"))
    }
  }
  
  dev.off()
  cat("Detailed", region_type, "individual plots saved\n")
  
  
  
}


# MAIN starts here #########################

#Load hypermethylated regions
cat("Loading hypermethylated regions from:", hyper_regions_file, "\n")
hyper_regions <- import.bed(hyper_regions_file)
cat("Loaded", length(hyper_regions), "hypermethylated regions\n")

# Load hypomethylated regions
cat("Loading hypomethylated regions from:", hypo_regions_file, "\n")
hypo_regions <- import.bed(hypo_regions_file)
cat("Loaded", length(hypo_regions), "hypomethylated regions\n")

# DM.R targeted output not loaded (set to NULL)
dmr_results <- NULL

# Check if regions have proper names
if (is.null(names(hyper_regions))) {
  names(hyper_regions) <- paste0("Hyper_", seq_along(hyper_regions))
}

if (is.null(names(hypo_regions))) {
  names(hypo_regions) <- paste0("Hypo_", seq_along(hypo_regions))
}

# Run methyKit style region test
results <-run_methylKit_style_region_test(enhancer_bed_file, bsseq_obj, output_dir, comparison_tag, genome)

cat("Creating methylation profile plots for both hyper and hypomethylated enhancers...\n")
# Plot hypermethylated regions
plot_enhancer_regions(hyper_regions, "Hypermethylated", bsseq_obj, output_dir, sample_groups, results)

# Plot hypomethylated regions
plot_enhancer_regions(hypo_regions, "Hypomethylated", bsseq_obj, output_dir, sample_groups, results)



cat("Script completed successfully!\n")
cat("Output files saved in:", output_dir, "\n")
