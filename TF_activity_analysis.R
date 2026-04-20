rm(list = ls())
# --constraint=avx2 for using of newer CPU instructions on H4H
# ==============================================================================
# Package Installation and Loading
# ==============================================================================

# Function to check and install CRAN packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message("Installing package: ", pkg)
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Function to check and install Bioconductor packages
install_bioc_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message("Installing Bioconductor package: ", pkg)
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Install and load CRAN packages
cran_packages <- c("data.table", "openxlsx", "dplyr", "tidyr", "stringr", "ggplot2", "pheatmap", "ggrepel", "tibble", "nnls", "limSolve")
for (pkg in cran_packages) {
  install_if_missing(pkg)
}

# Install and load Bioconductor packages
bioc_packages <- c("viper", "org.Hs.eg.db", "mygene", "decoupleR", "OmnipathR", "limma", "netZooR", "pcaMethods", "singscore", "SummarizedExperiment")
for (pkg in bioc_packages) {
  install_bioc_if_missing(pkg)
}

# Install DeconRNASeq from local zip file
rnaseq_base_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/RNAseq/deseq2"
if (!require("DeconRNASeq", character.only = TRUE, quietly = TRUE)) {
  message("Installing DeconRNASeq from local zip file")
  deconrnaseq_zip <- file.path(rnaseq_base_dir, "DeconRNASeq_1.52.0.zip")
  if (file.exists(deconrnaseq_zip)) {
    install.packages(deconrnaseq_zip, repos = NULL, type = "source")
    library(DeconRNASeq)
    cat("DeconRNASeq installed successfully from:", deconrnaseq_zip, "\n")
  } else {
    warning("DeconRNASeq zip file not found at:", deconrnaseq_zip)
  }
}

# Load regulatory network
net <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
netCorr <- decoupleR::check_corr(net)


# Load gene expression data
gene_expression <- read.csv(file.path(rnaseq_base_dir, "rsem_gene_level_abundance.tsv"), 
                            row.names = 1, 
                            sep = "\t")

# Extract gene IDs (Ensembl) and log2-transform
gene_symbols <- gene_expression[, 1]
sum(grepl("ZBTB48", gene_symbols))
log2_expression <- log2(gene_expression[, -1] + 1)

# Convert to matrix and set rownames to gene symbols
expression_matrix <- as.matrix(log2_expression)
rownames(expression_matrix) <- gene_symbols

# Verify the data
cat("Gene expression matrix dimensions:", dim(expression_matrix), "\n")
cat("First few rows:\n")
print(head(expression_matrix[grepl("ZNF", rownames(expression_matrix)), ]))

sample_info <- data.frame(sample = colnames(expression_matrix),
                          group = gsub("_.*", "", colnames(expression_matrix)))

# run univariate linear model
res_ulm <- decoupleR::run_ulm(expression_matrix, 
                   net, 
                   .source='source', 
                   .target='target',
                   .mor='mor', minsize = 5)
res_ulm
ntop <- 20
# Add group information
res_ulm <- res_ulm %>%
  mutate(group = gsub("_.*", "", condition))  # Extract group from sample name
 
# Aggregate scores by group (mean across replicates)
res_ulm_grouped <- res_ulm %>%
  group_by(source, group) %>%
  summarise(
    score = mean(score),
    p_value = mean(p_value),  # or use median, or meta-analysis
    n_replicates = n(),
    .groups = "drop"
  )
 
 # select to tops based on min p_value
tops <- res_ulm_grouped %>%
  group_by(source) %>%
  summarise(minp = min(p_value)) %>%
  arrange(minp) %>%
  head(ntop) %>%
  pull(source) 
# Create heatmap matrix
mat_ulm_grouped <- res_ulm_grouped %>%
  pivot_wider_profile(id_cols = source, 
              names_from = group, 
              values_from = score) %>%
  as.matrix()

pdf(file.path(rnaseq_base_dir, "decoupleR_ulm_analysis_TFs.pdf"), width=8, height=10)
pheatmap::pheatmap(mat_ulm_grouped[tops, ], cluster_rows = TRUE, cluster_cols = TRUE, cellwidth = 15, 
         border_color = NA, cellheight = 10, main = "Transcription factor : ulm")
dev.off()

# ==============================================================================
# VIPER Analysis
# ==============================================================================

# Run VIPER (Virtual Inference of Protein-activity by Enriched Regulon)
res_viper <- decoupleR::run_viper(expression_matrix, 
                                   net, 
                                   .source='source', 
                                   .target='target',
                                   .mor='mor', 
                                   minsize = 5)
res_viper

# Add group information
res_viper <- res_viper %>%
  mutate(group = gsub("_.*", "", condition))

# Aggregate scores by group (mean across replicates)
res_viper_grouped <- res_viper %>%
  group_by(source, group) %>%
  summarise(
    score = mean(score),
    p_value = mean(p_value),
    n_replicates = n(),
    .groups = "drop"
  )

# Select top TFs based on min p_value
tops_viper <- res_viper_grouped %>%
  group_by(source) %>%
  summarise(minp = min(p_value)) %>%
  arrange(minp) %>%
  head(ntop) %>%
  pull(source)

# Create heatmap matrix
mat_viper_grouped <- res_viper_grouped %>%
  pivot_wider_profile(id_cols = source, 
              names_from = group, 
              values_from = score) %>%
  as.matrix()

# Generate VIPER heatmap
pdf(file.path(rnaseq_base_dir, "decoupleR_viper_analysis_TFs.pdf"), width=8, height=10)
pheatmap::pheatmap(mat_viper_grouped[tops_viper, ], 
                   cluster_rows = TRUE, 
                   cluster_cols = TRUE, 
                   cellwidth = 15, 
                   border_color = NA, 
                   cellheight = 10, 
                   main = "Transcription factor : VIPER")
dev.off()


# ==================================================================================
# tiger analysis through netZooR
# ==================================================================================


# Create adjacency matrix: TFs (rows) x Target genes (columns) decoupleR::get_collectri
# Values represent mode of regulation: 1 (activation), -1 (repression), 0 (no interaction)
 
build_adjacency_matrix <- function(network) {
  
  # Get unique TFs and target genes
  tfs <- unique(network$source)
  targets <- unique(network$target)
  
  cat("Building adjacency matrix:\n")
  cat("  Number of TFs:", length(tfs), "\n")
  cat("  Number of target genes:", length(targets), "\n")
  
  # Initialize matrix with zeros
  adj_matrix <- matrix(0, 
                       nrow = length(tfs), 
                       ncol = length(targets),
                       dimnames = list(tfs, targets))
  
  # Fill in the interactions
  for (i in seq_len(nrow(network))) {
    tf <- network$source[i]
    target <- network$target[i]
    mor <- network$mor[i]  # mode of regulation
    
    adj_matrix[tf, target] <- mor
  }
  
  cat("  Matrix dimensions:", dim(adj_matrix), "\n")
  cat("  Non-zero entries:", sum(adj_matrix != 0), "\n")
  cat("  Activating interactions (+1):", sum(adj_matrix == 1), "\n")
  cat("  Repressing interactions (-1):", sum(adj_matrix == -1), "\n")
  
  return(adj_matrix)
}
 
# Build the adjacency matrix
adjacency_matrix <- build_adjacency_matrix(net)
 
# Optional: Save to file
write.csv(adjacency_matrix, 
          file.path(rnaseq_base_dir, "TF_target_adjacency_matrix.csv"),
          row.names = TRUE)
 
cat("Adjacency matrix saved to:", file.path(rnaseq_base_dir, "TF_target_adjacency_matrix.csv"), "\n")

# run tiger
if(!file.exists(file.path(rnaseq_base_dir, "tiger_result.rds"))) {
  tiger_result <- netZooR::tiger(expression_matrix, adjacency_matrix)
  saveRDS(tiger_result, file.path(rnaseq_base_dir, "tiger_result.rds"))
}
tiger_result <- readRDS(file.path(rnaseq_base_dir, "tiger_result.rds"))
names(tiger_result)

# ==============================================================================
# Limma Differential Analysis for TF Activities
# ==============================================================================

library(limma)

# Function to perform limma differential analysis
run_limma_tf_analysis <- function(res_data, method_name) {
  
  cat("\n=== Running limma differential analysis for", method_name, "===\n")
  
  # Create matrix: TFs (rows) x samples (columns)
  if(method_name == "TIGER") {
    tf_matrix <- res_data
  } else {
    tf_matrix <- res_data %>%
      dplyr::select(source, condition, score) %>%
      tidyr::pivot_wider(names_from = condition, values_from = score) %>%
      tibble::column_to_rownames("source") %>%
      as.matrix()
  }
  
  # Ensure columns match sample_info order
  tf_matrix <- tf_matrix[, sample_info$sample]
  
  # Create design matrix
  groups <- factor(sample_info$group)
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  
  cat("Design matrix:\n")
  print(design)
  
  # Fit linear model
  fit <- lmFit(tf_matrix, design)
  
  # Define contrasts - comparing each treatment to control (NIR)
  # Adjust contrast names based on your actual groups
  group_levels <- levels(groups)
  contrasts_list <- paste0(group_levels[group_levels != "NIR"], "-NIR")
  
  contrast_matrix <- makeContrasts(contrasts = contrasts_list, levels = design)
  
  cat("Contrast matrix:\n")
  print(contrast_matrix)
  
  # Fit contrasts
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  # Extract results for all contrasts
  all_results <- list()
  
  for (i in 1:ncol(contrast_matrix)) {
    contrast_name <- colnames(contrast_matrix)[i]
    cat("\n--- Results for contrast:", contrast_name, "---\n")
    
    # Get top table
    tt <- topTable(fit2, coef = i, number = Inf, adjust.method = "BH")
    tt$TF <- rownames(tt)
    tt$contrast <- contrast_name
    
    # Count significant TFs
    sig_up <- sum(tt$adj.P.Val < 0.05 & tt$logFC >= 1)
    sig_down <- sum(tt$adj.P.Val < 0.05 & tt$logFC <= -1)
    
    cat("Significant TFs (adj.P.Val < 0.05, |logFC| >= 1):\n")
    cat("  Up-regulated:", sig_up, "\n")
    cat("  Down-regulated:", sig_down, "\n")
    
    all_results[[contrast_name]] <- tt
  }
  
  # Combine all results
  combined_results <- bind_rows(all_results) 

  out_results <- combined_results %>%
    dplyr::filter(adj.P.Val < 0.05 & abs(logFC) >= 1)
  
  
  # order tf_matrix by max row score and get the top 20 rows
  top_rows <- order(apply(tf_matrix, 1, max), decreasing = TRUE)[1:20]
  top_tf_matrix <- tf_matrix[top_rows, ]
  # Create heatmap of top 20 TFs
  pdf_file <- file.path(rnaseq_base_dir, paste0("limma_TF_differential_", method_name, "_heatmap.pdf"))
  pdf(pdf_file, width = 10, height = 10)
  pheatmap::pheatmap(top_tf_matrix, 
                     cluster_rows = TRUE, 
                     cluster_cols = TRUE,
                     main = paste0("TFA for top 20 TFs - ", method_name))

  if (!is.null(out_results) && nrow(out_results) > 0) {
    
    # get significant TFs
    sub_tf_matrix <- tf_matrix[unique(out_results$TF), ]
    # scale by row
    sub_tf_matrix_scaled <- t(scale(t(sub_tf_matrix), center = TRUE, scale = TRUE))
    
    pheatmap::pheatmap(sub_tf_matrix, 
                       cluster_rows = TRUE, 
                       cluster_cols = TRUE,
                       main = paste0("TFA for significant TFs - ", method_name))
    pheatmap::pheatmap(sub_tf_matrix_scaled, 
                       cluster_rows = TRUE, 
                       cluster_cols = TRUE,
                       main = paste0("TFA for significant TFs (scaled) - ", method_name))
    # Save results
    output_file <- file.path(rnaseq_base_dir, paste0("limma_TF_differential_", method_name, ".csv"))
    write.csv(out_results, output_file, row.names = FALSE)
    cat("\nResults saved to:", output_file, "\n")
  } else {
    cat("\nNo significant TFs found - skipping heatmaps for significant TFs\n")
  }
  cat("\nHeatmap saved to:", pdf_file, "\n")
  dev.off()

  
  
  return(list(
    fit = fit2,
    results = combined_results,
    tf_matrix = tf_matrix,
    design = design
  ))
}


# Run limma analysis for ULM
ulm_limma <- run_limma_tf_analysis(res_ulm, "ULM")

# Run limma analysis for VIPER
viper_limma <- run_limma_tf_analysis(res_viper, "VIPER")

# Run limma analysis for TIGER
tiger_limma <- run_limma_tf_analysis(tiger_result$Z, "TIGER")

# ==============================================================================
# Visualize Significant TFs
# ==============================================================================

# Function to create volcano plot
create_volcano_plot <- function(results_df, contrast_name, method_name) {
  
  contrast_data <- results_df %>%
    filter(contrast == contrast_name) %>%
    mutate(
      significance = case_when(
        adj.P.Val < 0.05 & logFC > 1 ~ "Up",
        adj.P.Val < 0.05 & logFC < -1 ~ "Down",
        adj.P.Val < 0.05 ~ "Significant",
        TRUE ~ "NS"
      )
    )
  
  p <- ggplot(contrast_data, aes(x = logFC, y = -log10(adj.P.Val), color = significance)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", 
                                   "Significant" = "orange", "NS" = "grey")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    labs(title = paste0(method_name, ": ", contrast_name),
         x = "log2 Fold Change",
         y = "-log10(adj.P.Val)") +
    theme_bw() +
    theme(legend.position = "right")
  
  # Add labels for top TFs
  top_tfs <- contrast_data %>%
    filter(adj.P.Val < 0.05) %>%
    arrange(adj.P.Val) %>%
    head(10)
  
  if (nrow(top_tfs) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = top_tfs,
      aes(label = TF),
      size = 3,
      max.overlaps = 20
    )
  }
  
  return(p)
}

# Generate volcano plots for each contrast
ulm_contrasts <- unique(ulm_limma$results$contrast)
viper_contrasts <- unique(viper_limma$results$contrast)
tiger_contrasts <- unique(tiger_limma$results$contrast)

# ULM volcano plots
pdf(file.path(rnaseq_base_dir, "limma_ULM_volcano_plots.pdf"), width = 10, height = 8)
for (contrast in ulm_contrasts) {
  p <- create_volcano_plot(ulm_limma$results, contrast, "ULM")
  print(p)
}
dev.off()

# VIPER volcano plots
pdf(file.path(rnaseq_base_dir, "limma_VIPER_volcano_plots.pdf"), width = 10, height = 8)
for (contrast in viper_contrasts) {
  p <- create_volcano_plot(viper_limma$results, contrast, "VIPER")
  print(p)
}
dev.off()

# TIGER volcano plots
pdf(file.path(rnaseq_base_dir, "limma_TIGER_volcano_plots.pdf"), width = 10, height = 8)
for (contrast in tiger_contrasts) {
  p <- create_volcano_plot(tiger_limma$results, contrast, "TIGER")
  print(p)
}
dev.off()

cat("\n=== Limma differential analysis completed ===\n")


# ==============================================================================
# TF activity from DESeq2 results
# ==============================================================================

comp_dir <- file.path(rnaseq_base_dir, "single_factor_analysis_gene_level")
if (!dir.exists(comp_dir)) {
  stop("Comparison directory not found: ", comp_dir)
}

comparisons <- list.dirs(comp_dir, full.names = FALSE, recursive = FALSE)
cat("\nFound", length(comparisons), "comparisons:\n")
print(comparisons)

# Initialize list to store TF activity results
deseq2_tf_activities <- list()

# Loop through each comparison
for (comp in comparisons) {
  
  cat("\n=== Processing comparison:", comp, "===\n")
  
  # Construct path to DESeq2 results file
  deseq2_file <- file.path(comp_dir, comp, paste0(comp, "_DESeq2_results_regular.csv"))
  
  if (!file.exists(deseq2_file)) {
    cat("WARNING: DESeq2 results file not found:", deseq2_file, "\n")
    next
  }
  
  # Load DESeq2 results
  deseq2_res <- read.csv(deseq2_file)
  cat("Loaded DESeq2 results with", nrow(deseq2_res), "genes\n")
  
  # Check if stat column exists
  if (!"stat" %in% colnames(deseq2_res)) {
    cat("WARNING: 'stat' column not found in", deseq2_file, "\n")
    next
  }
  
  # Extract gene symbols and stat values
  # Assume row names are gene symbols or need to extract from first column
  gene_names <- deseq2_res$gene_name
  stat_values <- deseq2_res$stat
  
  # Remove NAs
  valid_idx <- !is.na(stat_values)
  gene_names <- gene_names[valid_idx]
  stat_values <- stat_values[valid_idx]
  
  cat("Valid genes after removing NAs:", length(gene_names), "\n")
  
  # Create a matrix for run_ulm (genes x samples, here just one "sample" = the stat column)
  stat_matrix <- matrix(stat_values, ncol = 1, dimnames = list(gene_names, comp))
  
  # Run ULM to infer TF activity
  cat("Running ULM to infer TF activity...\n")
  tf_activity <- decoupleR::run_ulm(
    mat = stat_matrix,
    network = net,
    .source = 'source',
    .target = 'target',
    .mor = 'mor',
    minsize = 5
  )
  
  cat("TF activity computed for", nrow(tf_activity), "TFs\n")
  
  # Add comparison name to results
  tf_activity$comparison <- comp
  
  # Store results
  deseq2_tf_activities[[comp]] <- tf_activity
  
  # Save individual comparison results
  output_file <- file.path(comp_dir, comp, paste0(comp, "_TF_activity_ULM.csv"))
  write.csv(tf_activity, output_file, row.names = FALSE)
  cat("TF activity saved to:", output_file, "\n")
}

# Combine all TF activities
if (length(deseq2_tf_activities) > 0) {
  combined_tf_activities <- dplyr::bind_rows(deseq2_tf_activities)
  
  # Save combined results
  combined_output <- file.path(rnaseq_base_dir, "all_comparisons_TF_activity_ULM.csv")
  write.csv(combined_tf_activities, combined_output, row.names = FALSE)
  cat("\nCombined TF activities saved to:", combined_output, "\n")
  
  # Create summary: TF activity matrix (TFs x comparisons)
  tf_activity_wide <- combined_tf_activities %>%
    dplyr::select(source, comparison, score) %>%
    tidyr::pivot_wider(names_from = comparison, values_from = score)
  
  # Save wide format
  wide_output <- file.path(rnaseq_base_dir, "TF_activity_matrix_comparisons.csv")
  write.csv(tf_activity_wide, wide_output, row.names = FALSE)
  cat("TF activity matrix saved to:", wide_output, "\n")
  
  # Generate heatmap of TF activities across comparisons
  tf_matrix_comp <- tf_activity_wide %>%
    tibble::column_to_rownames("source") %>%
    as.matrix()
  
  # Remove rows with all NAs
  tf_matrix_comp <- tf_matrix_comp[rowSums(is.na(tf_matrix_comp)) < ncol(tf_matrix_comp), ]
  
  # Select top TFs by variance across comparisons
  top_n_tfs <- 50
  if (nrow(tf_matrix_comp) > top_n_tfs) {
    tf_var <- apply(tf_matrix_comp, 1, var, na.rm = TRUE)
    top_tf_idx <- order(tf_var, decreasing = TRUE)[1:top_n_tfs]
    tf_matrix_comp_top <- tf_matrix_comp[top_tf_idx, ]
  } else {
    tf_matrix_comp_top <- tf_matrix_comp
  }
  
  # Generate heatmap
  pdf(file.path(rnaseq_base_dir, "TF_activity_heatmap_comparisons.pdf"), width = 12, height = 14)
  pheatmap::pheatmap(
    tf_matrix_comp_top,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "row",
    main = paste0("Top ", nrow(tf_matrix_comp_top), " TF Activities Across DESeq2 Comparisons"),
    fontsize_row = 8,
    cellwidth = 15,
    cellheight = 10
  )
  dev.off()
  
  cat("\nTF activity heatmap saved to:", file.path(rnaseq_base_dir, "TF_activity_heatmap_comparisons.pdf"), "\n")
  
} else {
  cat("\nNo TF activities computed - no valid comparisons found\n")
}


# ==============================================================================
# Luminal to Basal Epithelial Transition Analysis
# ==============================================================================

cat("\n=== Analyzing Luminal to Basal Epithelial Transition ===\n")

# Define marker genes
basal_markers <- c("KRT5", "KRT14", "KRT15", "TP63", "ITGA6", "ITGB4", "COL17A1", "LAMA3")
luminal_markers <- c("KRT8", "KRT18", "KRT19", "EPCAM", "CDH1", "MUC1")
basal_tfs <- c("TP63", "SOX2", "MYC", "TP53", "RUNX1", "GRHL3")

# Extract TPM expression data (not log2 transformed)
tpm_expression <- expression_matrix

# Get sample groups
sample_groups <- gsub("_.*", "", colnames(tpm_expression))

# Function to extract marker expression
extract_marker_expression <- function(markers, expr_data, gene_col = NULL) {
  available_markers <- markers[markers %in% rownames(expr_data)]
  if (length(available_markers) == 0) {
    cat("WARNING: No markers found in expression data\n")
    return(NULL)
  }
  cat("Found", length(available_markers), "out of", length(markers), "markers\n")
  cat("Available markers:", paste(available_markers, collapse = ", "), "\n")
  return(expr_data[available_markers, , drop = FALSE])
}

# Extract basal and luminal marker expression
cat("\n--- Basal markers ---\n")
basal_expr <- extract_marker_expression(basal_markers, tpm_expression)

cat("\n--- Luminal markers ---\n")
luminal_expr <- extract_marker_expression(luminal_markers, tpm_expression)

# Calculate mean expression for each marker type per sample
if (!is.null(basal_expr) && !is.null(luminal_expr)) {
  
  # Mean expression across markers
  basal_score <- colMeans(basal_expr, na.rm = TRUE)
  luminal_score <- colMeans(luminal_expr, na.rm = TRUE)
  
  # Calculate basal/luminal ratio (transition score)
  transition_score <- log2((basal_score + 1) / (luminal_score + 1))
  
  # Create summary data frame
  transition_df <- data.frame(
    sample = names(transition_score),
    group = sample_groups,
    basal_score = basal_score,
    luminal_score = luminal_score,
    transition_score = transition_score,
    stringsAsFactors = FALSE
  )
  
  # Save transition scores
  write.csv(transition_df, 
            file.path(rnaseq_base_dir, "epithelial_transition_scores.csv"),
            row.names = FALSE)
  
  cat("\nTransition scores saved\n")
  
  # Visualize marker expression heatmap
  pdf(file.path(rnaseq_base_dir, "epithelial_markers_heatmap.pdf"), width = 12, height = 10)
  
  # Combine basal and luminal markers
  all_markers <- rbind(basal_expr, luminal_expr)
  
  # Add annotation for marker type
  marker_annotation <- data.frame(
    Type = c(rep("Basal", nrow(basal_expr)), rep("Luminal", nrow(luminal_expr))),
    row.names = rownames(all_markers)
  )
  
  # Sample annotation
  sample_annotation <- data.frame(
    Group = sample_groups,
    row.names = colnames(all_markers)
  )
  
  # Log2 transform for visualization
  all_markers_log <- log2(all_markers + 1)
  
  pheatmap::pheatmap(
    all_markers_log,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_row = marker_annotation,
    annotation_col = sample_annotation,
    main = "Basal vs Luminal Epithelial Markers Expression (TPM)",
    fontsize_row = 10,
    cellwidth = 15,
    cellheight = 12
  )
  
  dev.off()
  cat("Marker expression heatmap saved\n")
  
  # Plot transition scores by group
  pdf(file.path(rnaseq_base_dir, "epithelial_transition_scores.pdf"), width = 10, height = 8)
  
  # Box plot of transition scores
  p1 <- ggplot(transition_df, aes(x = group, y = transition_score, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    labs(title = "Epithelial Transition Score (Basal/Luminal Ratio)",
         x = "Group",
         y = "log2(Basal/Luminal)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")
  print(p1)
  
  # Basal vs Luminal scores
  p2 <- ggplot(transition_df, aes(x = luminal_score, y = basal_score, color = group)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(title = "Basal vs Luminal Marker Expression",
         x = "Mean Luminal Marker Expression (TPM)",
         y = "Mean Basal Marker Expression (TPM)") +
    theme_bw() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray")
  print(p2)
  
  # Individual marker expression by group
  basal_long <- as.data.frame(basal_expr) %>% 
    tibble::rownames_to_column("gene") %>% 
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "TPM") %>% 
    dplyr::mutate(group = gsub("_.*", "", sample),
                  type = "Basal")
  
  luminal_long <- as.data.frame(luminal_expr) %>% 
    tibble::rownames_to_column("gene") %>% 
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "TPM") %>% 
    dplyr::mutate(group = gsub("_.*", "", sample),
                  type = "Luminal")
  
  marker_long <- rbind(basal_long, luminal_long)
  
  p3 <- ggplot(marker_long, aes(x = group, y = log2(TPM + 1), fill = type)) +
    geom_boxplot() +
    facet_wrap(~ gene, scales = "free_y", ncol = 4) +
    labs(title = "Individual Marker Expression Across Groups",
         x = "Group",
         y = "log2(TPM + 1)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p3)
  
  dev.off()
  cat("Transition score plots saved\n")
}

# Extract TF activities for basal-related TFs
cat("\n--- Basal TF Activities (from sample-based ULM) ---\n")

if (exists("res_ulm")) {
  basal_tf_activity <- res_ulm %>% 
    dplyr::filter(source %in% basal_tfs) %>% 
    dplyr::mutate(group = gsub("_.*", "", condition))
  
  if (nrow(basal_tf_activity) > 0) {
    cat("Found", length(unique(basal_tf_activity$source)), "basal TFs in ULM results\n")
    
    # Save basal TF activities
    write.csv(basal_tf_activity, 
              file.path(rnaseq_base_dir, "basal_TF_activities_ULM.csv"),
              row.names = FALSE)
    
    # Plot basal TF activities
    pdf(file.path(rnaseq_base_dir, "basal_TF_activities.pdf"), width = 12, height = 8)
    
    p_tf <- ggplot(basal_tf_activity, aes(x = group, y = score, fill = group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      facet_wrap(~ source, scales = "free_y", ncol = 3) +
      labs(title = "Basal TF Activities Across Groups (ULM)",
           x = "Group",
           y = "TF Activity Score") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    print(p_tf)
    
    dev.off()
    cat("Basal TF activity plots saved\n")
  }
}

# Extract TF activities from DESeq2-based analysis
cat("\n--- Basal TF Activities (from DESeq2 comparisons) ---\n")

if (exists("combined_tf_activities")) {
  basal_tf_deseq2 <- combined_tf_activities %>% 
    dplyr::filter(source %in% basal_tfs)
  
  if (nrow(basal_tf_deseq2) > 0) {
    cat("Found", length(unique(basal_tf_deseq2$source)), "basal TFs in DESeq2-based results\n")
    
    # Save
    write.csv(basal_tf_deseq2,
              file.path(rnaseq_base_dir, "basal_TF_activities_DESeq2.csv"),
              row.names = FALSE)
    
    # Plot
    pdf(file.path(rnaseq_base_dir, "basal_TF_activities_DESeq2.pdf"), width = 14, height = 8)
    
    p_tf_deseq2 <- ggplot(basal_tf_deseq2, aes(x = comparison, y = score, fill = source)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(title = "Basal TF Activities Across DESeq2 Comparisons",
           x = "Comparison",
           y = "TF Activity Score",
           fill = "TF") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p_tf_deseq2)
    
    # Heatmap of basal TFs
    basal_tf_wide <- basal_tf_deseq2 %>% 
      dplyr::select(source, comparison, score) %>% 
      tidyr::pivot_wider(names_from = comparison, values_from = score) %>% 
      tibble::column_to_rownames("source")
    
    pheatmap::pheatmap(
      as.matrix(basal_tf_wide),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      main = "Basal TF Activities - DESeq2 Comparisons",
      fontsize_row = 10,
      cellwidth = 12,
      cellheight = 15
    )
    
    dev.off()
    cat("DESeq2-based basal TF activity plots saved\n")
  }
}

cat("\n=== Epithelial transition analysis completed ===\n")
cat("\nInterpretation guide:\n")
cat("- Positive transition score: Shift toward basal phenotype\n")
cat("- Negative transition score: Shift toward luminal phenotype\n")
cat("- High TP63/SOX2 activity: Indicative of basal epithelial state\n")
cat("- Look for radiation-induced increases in basal markers and TF activities\n\n")

# ==============================================================================
# Cell Type Deconvolution: Basal vs Luminal Epithelial Cells
# ==============================================================================

cat("\n=== Deconvoluting Basal and Luminal Cell Types ===\n")

# Use log2-transformed expression for deconvolution
# expression_matrix is already log2(TPM + 1)

if (!is.null(basal_expr) && !is.null(luminal_expr)) {
  
  # Create signature matrix from marker genes
  # Rows = genes, Columns = cell types
  
  # Get log2 expression of markers
  basal_markers_log <- expression_matrix[rownames(basal_expr), , drop = FALSE]
  luminal_markers_log <- expression_matrix[rownames(luminal_expr), , drop = FALSE]
  
  # Create reference signatures (mean expression of markers for each cell type)
  # For a simple deconvolution, we'll use relative expression patterns
  
  # ==============================================================================
  # Cell Type Deconvolution using DeconRNASeq Package
  # ==============================================================================
  
  cat("\n=== Running DeconRNASeq Cell Type Deconvolution ===\n")
  
  # DeconRNASeq requires:
  # 1. Signature matrix: genes x cell types (reference profiles)
  # 2. Mixture matrix: genes x samples (your bulk expression data)
  
  # Get all marker genes
  all_marker_genes <- unique(c(rownames(basal_markers_log), rownames(luminal_markers_log)))
  
  # Create signature matrix with realistic expression profiles
  # Use actual mean expression across all samples as reference
  
  # Calculate mean expression for each marker set (in log2 space)
  basal_ref <- rowMeans(basal_markers_log, na.rm = TRUE)
  luminal_ref <- rowMeans(luminal_markers_log, na.rm = TRUE)
  
  # Convert back from log2(x+1) to TPM: 2^value - 1
  basal_ref_tpm <- 2^basal_ref - 1
  luminal_ref_tpm <- 2^luminal_ref - 1
  
  # Build signature matrix: rows = genes, cols = cell types
  # Each column represents the expression profile of that cell type
  sig_matrix <- matrix(0.01, nrow = length(all_marker_genes), ncol = 2)  # Small baseline instead of 0
  rownames(sig_matrix) <- all_marker_genes
  colnames(sig_matrix) <- c("Basal", "Luminal")
  
  # Basal cell signature: high basal markers, low luminal markers
  sig_matrix[names(basal_ref_tpm), "Basal"] <- basal_ref_tpm
  sig_matrix[names(luminal_ref_tpm), "Basal"] <- pmin(luminal_ref_tpm * 0.1, 0.1)  # Low but not zero
  
  # Luminal cell signature: low basal markers, high luminal markers  
  sig_matrix[names(basal_ref_tpm), "Luminal"] <- pmin(basal_ref_tpm * 0.1, 0.1)  # Low but not zero
  sig_matrix[names(luminal_ref_tpm), "Luminal"] <- luminal_ref_tpm
  
  # Prepare mixture matrix (your samples)
  # DeconRNASeq expects non-log transformed data (TPM scale)
  mixture_matrix <- 2^expression_matrix[all_marker_genes, ] - 1  # Convert from log2(x+1) back to TPM
  mixture_matrix[mixture_matrix < 0] <- 0  # Ensure no negative values
  
  cat("Signature matrix dimensions:", dim(sig_matrix), "\n")
  cat("Mixture matrix dimensions:", dim(mixture_matrix), "\n")
  cat("Signature matrix preview:\n")
  print(head(sig_matrix))
  cat("Signature matrix range - Basal:", range(sig_matrix[, "Basal"]), "\n")
  cat("Signature matrix range - Luminal:", range(sig_matrix[, "Luminal"]), "\n")
  
  # Run DeconRNASeq
  tryCatch({
    cat("\nRunning DeconRNASeq...\n")
    
    decon_result <- DeconRNASeq::DeconRNASeq(
      datasets = as.data.frame(mixture_matrix),
      signatures = as.data.frame(sig_matrix),
      proportions = NULL,
      checksig = FALSE,
      known.prop = FALSE,
      use.scale = TRUE,
      fig = FALSE
    )
    
    cat("DeconRNASeq completed. Checking results...\n")
    cat("Result structure:\n")
    print(names(decon_result))
    
    # Extract proportions
    deconrnaseq_props <- decon_result$out.all
    
    cat("Proportions matrix dimensions:", dim(deconrnaseq_props), "\n")
    
    # Check if we have valid results
    if (is.null(deconrnaseq_props) || nrow(deconrnaseq_props) == 0 || ncol(deconrnaseq_props) == 0) {
      stop("DeconRNASeq returned empty results. Check signature and mixture matrices.")
    }
    
    cat("Proportions preview:\n")
    print(head(deconrnaseq_props))
    cat("Row names of proportions:\n")
    print(rownames(deconrnaseq_props))
    cat("Column names of proportions:\n")
    print(colnames(deconrnaseq_props))
    
    # DeconRNASeq already returns: samples x cell types (15 samples x 2 cell types)
    # No transpose needed!
    
    # Create data frame directly
    deconrnaseq_df <- as.data.frame(deconrnaseq_props)
    
    cat("Data frame dimensions:", dim(deconrnaseq_df), "\n")
    cat("Data frame column names:", colnames(deconrnaseq_df), "\n")
    
    # Add sample names from the mixture matrix column names
    deconrnaseq_df$sample <- colnames(mixture_matrix)
    deconrnaseq_df$group <- gsub("_.*", "", deconrnaseq_df$sample)
    
    # Check for valid column names
    if (!("Basal" %in% colnames(deconrnaseq_df)) || !("Luminal" %in% colnames(deconrnaseq_df))) {
      cat("ERROR: Actual columns are:", paste(colnames(deconrnaseq_df), collapse=", "), "\n")
      stop("DeconRNASeq results missing expected cell type columns")
    }
    
    deconrnaseq_df$basal_luminal_ratio <- deconrnaseq_df$Basal / deconrnaseq_df$Luminal
    
    # Save results
    write.csv(deconrnaseq_df,
              file.path(rnaseq_base_dir, "cell_type_deconvolution_DeconRNASeq.csv"),
              row.names = FALSE)
    cat("DeconRNASeq results saved\n")
    
    # Print summary
    cat("\n--- DeconRNASeq Summary by Group ---\n")
    deconrnaseq_summary <- deconrnaseq_df %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        mean_basal = mean(Basal),
        sd_basal = sd(Basal),
        mean_luminal = mean(Luminal),
        sd_luminal = sd(Luminal),
        mean_ratio = mean(basal_luminal_ratio),
        n = n(),
        .groups = "drop"
      )
    print(deconrnaseq_summary)
    
    # Visualize deconvolution results
    cat("\n--- Creating Visualizations ---\n")
    pdf(file.path(rnaseq_base_dir, "cell_type_deconvolution_DeconRNASeq.pdf"), width = 12, height = 10)
    
    # Stacked bar plot of cell type proportions
    deconv_long <- deconrnaseq_df %>%
      dplyr::select(sample, group, Basal, Luminal) %>%
      tidyr::pivot_longer(cols = c(Basal, Luminal),
                          names_to = "cell_type",
                          values_to = "proportion")
    
    p1 <- ggplot(deconv_long, aes(x = sample, y = proportion, fill = cell_type)) +
      geom_bar(stat = "identity") +
      facet_wrap(~ group, scales = "free_x") +
      scale_fill_manual(values = c("Basal" = "#E69F00", "Luminal" = "#56B4E9")) +
      labs(title = "Cell Type Deconvolution: Basal vs Luminal (DeconRNASeq)",
           x = "Sample",
           y = "Proportion",
           fill = "Cell Type") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
    print(p1)
    
    # Box plot of basal proportion by group
    p2 <- ggplot(deconrnaseq_df, aes(x = group, y = Basal, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      labs(title = "Basal Cell Proportion by Group",
           x = "Group",
           y = "Basal Cell Proportion") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")
    print(p2)
    
    # Box plot of basal/luminal ratio by group
    p3 <- ggplot(deconrnaseq_df, aes(x = group, y = log2(basal_luminal_ratio), fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      labs(title = "Basal/Luminal Ratio by Group",
           x = "Group",
           y = "log2(Basal/Luminal Ratio)") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red")
    print(p3)
    
    # Scatter plot: basal vs luminal proportions
    p4 <- ggplot(deconrnaseq_df, aes(x = Luminal, y = Basal, color = group)) +
      geom_point(size = 3, alpha = 0.7) +
      labs(title = "Basal vs Luminal Cell Proportions",
           x = "Luminal Proportion",
           y = "Basal Proportion",
           color = "Group") +
      theme_bw() +
      geom_abline(slope = -1, intercept = 1, linetype = "dashed", color = "gray") +
      coord_fixed()
    print(p4)
    
    # Heatmap of proportions
    deconrnaseq_matrix <- t(as.matrix(deconrnaseq_df[, c("Basal", "Luminal")]))
    colnames(deconrnaseq_matrix) <- deconrnaseq_df$sample
    
    sample_annot <- data.frame(
      Group = deconrnaseq_df$group,
      row.names = deconrnaseq_df$sample
    )
    
    pheatmap::pheatmap(
      deconrnaseq_matrix,
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      annotation_col = sample_annot,
      main = "Cell Type Proportions Across Samples",
      color = colorRampPalette(c("white", "darkblue"))(100),
      breaks = seq(0, 1, length.out = 101),
      cellwidth = 10,
      cellheight = 20
    )
    
    dev.off()
    cat("Deconvolution plots saved\n")
    
    # Statistical test: compare basal proportion across groups
    cat("\n--- Statistical Testing ---\n")
    
    if (length(unique(deconrnaseq_df$group)) > 2) {
      # Kruskal-Wallis test (non-parametric)
      kw_test <- kruskal.test(Basal ~ group, data = deconrnaseq_df)
      cat("Kruskal-Wallis test for Basal proportion across groups:\n")
      cat("  Chi-squared =", kw_test$statistic, "\n")
      cat("  p-value =", kw_test$p.value, "\n")
      
      # Post-hoc pairwise comparisons if significant
      if (kw_test$p.value < 0.05) {
        cat("\nPerforming pairwise Wilcoxon tests...\n")
        pairwise_results <- pairwise.wilcox.test(deconrnaseq_df$Basal, 
                                                  deconrnaseq_df$group, 
                                                  p.adjust.method = "BH")
        print(pairwise_results)
      }
    }
    
    cat("\n=== Cell Type Deconvolution Completed ===\n")
    cat("\nInterpretation:\n")
    cat("- Basal proportion > 0.5: Sample dominated by basal cells\n")
    cat("- Luminal proportion > 0.5: Sample dominated by luminal cells\n")
    cat("- Look for radiation-induced shifts in cell type composition\n")
    cat("- Higher basal/luminal ratio indicates enrichment of basal epithelial phenotype\n\n")
    
  }, error = function(e) {
    cat("ERROR running DeconRNASeq:", e$message, "\n")
    cat("Cell type deconvolution failed. Please check input data and marker genes.\n")
  })
  
} else {
  cat("Skipping deconvolution - marker expression data not available\n")
}

# ==============================================================================
# Oxygen Depletion and cGAS-STING Pathway Analysis
# ==============================================================================

cat("\n=== Analyzing Oxygen Depletion and cGAS-STING Pathway Activation ===\n")

# Define pathway genes based on literature
# Oxygen depletion / Hypoxia response pathway genes
oxygen_depletion_genes <- c(
  # HIF pathway (master regulators of hypoxia response)
  "HIF1A", "HIF2A", "EPAS1", "ARNT", "HIF3A",
  # HIF target genes
  "VEGFA", "VEGFB", "VEGFC", "LDHA", "LDHB", "PDK1", "PDK2", "PDK3", "PDK4",
  "SLC2A1", "SLC2A3", "GLUT1", "GLUT3",  # Glucose transporters
  "PGK1", "ENO1", "ENO2", "ALDOA", "ALDOC",  # Glycolytic enzymes
  "BNIP3", "BNIP3L", "NIX",  # Autophagy/mitophagy
  "EPO", "HMOX1", "NOS2", "NOS3",  # Oxygen sensing/metabolism
  "CA9", "CA12",  # Carbonic anhydrases (hypoxia markers)
  "ADM", "ANGPTL4", "LOX", "LOXL2",  # Hypoxia-induced genes
  "PHD1", "PHD2", "PHD3", "EGLN1", "EGLN2", "EGLN3",  # Prolyl hydroxylases
  "VHL",  # VHL tumor suppressor
  "CITED2", "DEC1", "DEC2", "BHLHE40", "BHLHE41"  # HIF modulators
)

# cGAS-STING pathway genes
cgas_sting_genes <- c(
  # Core pathway components
  "CGAS", "MB21D1",  # cGAS (also known as MB21D1)
  "STING1", "TMEM173",  # STING (also known as TMEM173)
  "TBK1",  # TANK-binding kinase 1
  "IRF3", "IRF7",  # Interferon regulatory factors
  # Type I interferons
  "IFNA1", "IFNA2", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8",
  "IFNA10", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA21",
  "IFNB1",  # IFN-beta
  "IFNL1", "IFNL2", "IFNL3",  # Type III interferons
  # Interferon-stimulated genes (ISGs)
  "ISG15", "ISG20", "IFIT1", "IFIT2", "IFIT3", "IFIT5",
  "OAS1", "OAS2", "OAS3", "OASL",
  "MX1", "MX2",
  "RSAD2",  # Viperin
  "IFI44", "IFI44L", "IFI27", "IFI6", "IFI35",
  "IFIH1", "DDX58",  # MDA5 and RIG-I
  "STAT1", "STAT2", "JAK1", "TYK2",  # JAK-STAT signaling
  # Inflammatory cytokines downstream of STING
  "IL6", "TNF", "CXCL10", "CCL5", "CXCL9", "CXCL11",
  # DNA damage sensors that feed into cGAS-STING
  "TREX1",  # DNA exonuclease (negative regulator)
  "ATM", "ATR", "CHEK1", "CHEK2"  # DNA damage response
)

# Load DESeq2 comparison results
comp_dir <- file.path(rnaseq_base_dir, "single_factor_analysis_gene_level")

if (dir.exists(comp_dir)) {
  
  comparisons <- list.dirs(comp_dir, full.names = FALSE, recursive = FALSE)
  # Filter for radiation vs NIR comparisons, matching "Gy" and "NIR"
  rad_comparisons <- comparisons[grepl("Gy", comparisons) & grepl("NIR", comparisons)]
  
  cat("Found", length(rad_comparisons), "radiation vs NIR comparisons:\n")
  print(rad_comparisons)
  
  # Initialize lists to store pathway gene stats
  oxygen_stats_list <- list()
  cgas_sting_stats_list <- list()
  
  # Loop through each comparison
  for (comp in rad_comparisons) {
    
    cat("\n--- Processing:", comp, "---\n")
    
    # Load DESeq2 results
    deseq2_file <- file.path(comp_dir, comp, paste0(comp, "_DESeq2_results_regular.csv"))
    
    if (!file.exists(deseq2_file)) {
      cat("WARNING: File not found:", deseq2_file, "\n")
      next
    }
    
    deseq2_res <- read.csv(deseq2_file)
    
    # Extract oxygen depletion pathway genes
    oxygen_data <- deseq2_res %>%
      dplyr::filter(gene_name %in% oxygen_depletion_genes) %>%
      dplyr::select(gene_name, stat, log2FoldChange, pvalue, padj) %>%
      dplyr::mutate(comparison = comp, pathway = "Oxygen_Depletion")
    
    cat("Oxygen depletion genes found:", nrow(oxygen_data), "\n")
    
    # Extract cGAS-STING pathway genes
    cgas_data <- deseq2_res %>%
      dplyr::filter(gene_name %in% cgas_sting_genes) %>%
      dplyr::select(gene_name, stat, log2FoldChange, pvalue, padj) %>%
      dplyr::mutate(comparison = comp, pathway = "cGAS_STING")
    
    cat("cGAS-STING genes found:", nrow(cgas_data), "\n")
    
    # Store results
    oxygen_stats_list[[comp]] <- oxygen_data
    cgas_sting_stats_list[[comp]] <- cgas_data
  }
  
  # Combine all results
  oxygen_stats <- dplyr::bind_rows(oxygen_stats_list)
  cgas_sting_stats <- dplyr::bind_rows(cgas_sting_stats_list)
  
  # Save pathway gene statistics
  write.csv(oxygen_stats, 
            file.path(rnaseq_base_dir, "oxygen_depletion_pathway_stats.csv"),
            row.names = FALSE)
  write.csv(cgas_sting_stats, 
            file.path(rnaseq_base_dir, "cGAS_STING_pathway_stats.csv"),
            row.names = FALSE)
  
  cat("\nPathway statistics saved\n")
  
  # ==============================================================================
  # Summarize Pathway Activation
  # ==============================================================================
  
  cat("\n=== Oxygen Depletion Pathway Summary ===\n")
  
  # Summary by comparison
  oxygen_summary <- oxygen_stats %>%
    dplyr::group_by(comparison) %>%
    dplyr::summarise(
      n_genes = n(),
      mean_stat = mean(stat, na.rm = TRUE),
      median_stat = median(stat, na.rm = TRUE),
      n_upregulated = sum(stat > 0 & padj < 0.05, na.rm = TRUE),
      n_downregulated = sum(stat < 0 & padj < 0.05, na.rm = TRUE),
      mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(mean_stat))
  
  print(oxygen_summary)
  
  cat("\n=== cGAS-STING Pathway Summary ===\n")
  
  cgas_summary <- cgas_sting_stats %>%
    dplyr::group_by(comparison) %>%
    dplyr::summarise(
      n_genes = n(),
      mean_stat = mean(stat, na.rm = TRUE),
      median_stat = median(stat, na.rm = TRUE),
      n_upregulated = sum(stat > 0 & padj < 0.05, na.rm = TRUE),
      n_downregulated = sum(stat < 0 & padj < 0.05, na.rm = TRUE),
      mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(mean_stat))
  
  print(cgas_summary)
  
  # ==============================================================================
  # Visualizations
  # ==============================================================================
  
  pdf(file.path(rnaseq_base_dir, "pathway_activation_analysis.pdf"), width = 14, height = 12)
  
  # 1. Heatmap of oxygen depletion pathway genes
  if (nrow(oxygen_stats) > 0) {
    
    oxygen_wide <- oxygen_stats %>%
      dplyr::select(gene_name, comparison, stat) %>%
      tidyr::pivot_wider(names_from = comparison, values_from = stat) %>%
      tibble::column_to_rownames("gene_name")
    
    # Remove rows with all NAs
    oxygen_wide <- oxygen_wide[rowSums(!is.na(oxygen_wide)) > 0, , drop = FALSE]
    
    if (nrow(oxygen_wide) > 0) {
      pheatmap::pheatmap(
        as.matrix(oxygen_wide),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        scale = "none",
        main = "Oxygen Depletion / Hypoxia Pathway Genes (DESeq2 stat)",
        fontsize_row = 8,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        breaks = seq(-10, 10, length.out = 101),
        na_col = "grey90"
      )
    }
  }
  
  # 2. Heatmap of cGAS-STING pathway genes
  if (nrow(cgas_sting_stats) > 0) {
    
    cgas_wide <- cgas_sting_stats %>%
      dplyr::select(gene_name, comparison, stat) %>%
      tidyr::pivot_wider(names_from = comparison, values_from = stat) %>%
      tibble::column_to_rownames("gene_name")
    
    # Remove rows with all NAs
    cgas_wide <- cgas_wide[rowSums(!is.na(cgas_wide)) > 0, , drop = FALSE]
    
    if (nrow(cgas_wide) > 0) {
      pheatmap::pheatmap(
        as.matrix(cgas_wide),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        scale = "none",
        main = "cGAS-STING Pathway Genes (DESeq2 stat)",
        fontsize_row = 8,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        breaks = seq(-10, 10, length.out = 101),
        na_col = "grey90"
      )
    }
  }
  
  # 3. Mean +/- SE bar plot with p-value significance for each pathway
  all_pathway_stats <- rbind(oxygen_stats, cgas_sting_stats)
  
  # Calculate mean, SE, and p-value for each pathway/comparison
  pathway_summary_stats <- all_pathway_stats %>%
    dplyr::group_by(pathway, comparison) %>%
    dplyr::summarise(
      mean_stat = mean(stat, na.rm = TRUE),
      se_stat = sd(stat, na.rm = TRUE) / sqrt(sum(!is.na(stat))),
      n = sum(!is.na(stat)),
      # One-sample t-test against 0
      p_value = ifelse(n >= 3, t.test(stat, mu = 0)$p.value, NA),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # Significance labels
      sig_label = case_when(
        is.na(p_value) ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      # Position for significance label (above or below bar)
      label_y = ifelse(mean_stat >= 0, mean_stat + se_stat + 0.3, mean_stat - se_stat - 0.3)
    )
  
  # Faceted bar plot with mean +/- SE and significance
  p1 <- ggplot(pathway_summary_stats, aes(x = comparison, y = mean_stat, fill = pathway)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = mean_stat - se_stat, ymax = mean_stat + se_stat),
                  width = 0.25) +
    geom_text(aes(y = label_y, label = sig_label), size = 4, vjust = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ pathway, ncol = 1, scales = "free_y") +
    labs(title = "Pathway Activation (Mean ± SE of DESeq2 stat)",
         subtitle = "* p<0.05, ** p<0.01, *** p<0.001 (one-sample t-test vs 0)",
         x = "Comparison",
         y = "Mean DESeq2 stat") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    scale_fill_manual(values = c("Oxygen_Depletion" = "#E69F00", "cGAS_STING" = "#56B4E9"))
  print(p1)
  
  # 5. Highlight key genes in each pathway
  # Oxygen depletion key genes
  oxygen_key_genes <- c("HIF1A", "EPAS1", "VEGFA", "LDHA", "PDK1", "CA9", "BNIP3", "EGLN1")
  oxygen_key <- oxygen_stats %>% dplyr::filter(gene_name %in% oxygen_key_genes)
  
  if (nrow(oxygen_key) > 0) {
    p3 <- ggplot(oxygen_key, aes(x = comparison, y = stat, fill = gene_name)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Key Oxygen Depletion / Hypoxia Genes",
           x = "Comparison",
           y = "DESeq2 stat",
           fill = "Gene") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p3)
  }
  
  # cGAS-STING key genes
  cgas_key_genes <- c("CGAS", "MB21D1", "STING1", "TMEM173", "TBK1", "IRF3", "IRF7", 
                      "IFNB1", "ISG15", "CXCL10", "STAT1", "MX1")
  cgas_key <- cgas_sting_stats %>% dplyr::filter(gene_name %in% cgas_key_genes)
  
  if (nrow(cgas_key) > 0) {
    p4 <- ggplot(cgas_key, aes(x = comparison, y = stat, fill = gene_name)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Key cGAS-STING Pathway Genes",
           x = "Comparison",
           y = "DESeq2 stat",
           fill = "Gene") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p4)
  }
  
  # 6. Volcano-style plot for significant pathway genes
  p5 <- ggplot(all_pathway_stats, aes(x = log2FoldChange, y = -log10(padj), color = pathway)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    facet_wrap(~ comparison, ncol = 2) +
    labs(title = "Pathway Genes: log2FC vs Significance",
         x = "log2 Fold Change",
         y = "-log10(padj)",
         color = "Pathway") +
    theme_bw() +
    scale_color_manual(values = c("Oxygen_Depletion" = "#E69F00", "cGAS_STING" = "#56B4E9"))
  print(p5)
  
  dev.off()
  
  cat("\nPathway activation plots saved to:", file.path(rnaseq_base_dir, "pathway_activation_analysis.pdf"), "\n")
  
  # ==============================================================================
  # Statistical Testing for Pathway Activation
  # ==============================================================================
  
  cat("\n=== Statistical Testing for Pathway Activation ===\n")
  
  # Test if mean stat is significantly different from 0 for each pathway/comparison
  cat("\n--- Oxygen Depletion Pathway ---\n")
  for (comp in unique(oxygen_stats$comparison)) {
    comp_data <- oxygen_stats %>% dplyr::filter(comparison == comp)
    if (nrow(comp_data) >= 3) {
      t_result <- t.test(comp_data$stat, mu = 0)
      cat(sprintf("%s: mean stat = %.3f, p-value = %.4f %s\n", 
                  comp, 
                  mean(comp_data$stat, na.rm = TRUE),
                  t_result$p.value,
                  ifelse(t_result$p.value < 0.05, "*", "")))
    }
  }
  
  cat("\n--- cGAS-STING Pathway ---\n")
  for (comp in unique(cgas_sting_stats$comparison)) {
    comp_data <- cgas_sting_stats %>% dplyr::filter(comparison == comp)
    if (nrow(comp_data) >= 3) {
      t_result <- t.test(comp_data$stat, mu = 0)
      cat(sprintf("%s: mean stat = %.3f, p-value = %.4f %s\n", 
                  comp, 
                  mean(comp_data$stat, na.rm = TRUE),
                  t_result$p.value,
                  ifelse(t_result$p.value < 0.05, "*", "")))
    }
  }
  
  cat("\n=== DESeq2 stat-based Pathway Analysis Completed ===\n")
  
  # ==============================================================================
  # Singscore Analysis for Pathway Activation
  # ==============================================================================
  
  cat("\n=== Computing Singscore for Pathway Activation ===\n")
  
  # Singscore is a rank-based single-sample gene set scoring method
  # It doesn't require a reference and works on individual samples
  
  library(singscore)
  library(SummarizedExperiment)
  
  # Create gene sets as named lists (GSEABase format not required for simpleScore)
  oxygen_geneset <- oxygen_depletion_genes[oxygen_depletion_genes %in% rownames(expression_matrix)]
  cgas_geneset <- cgas_sting_genes[cgas_sting_genes %in% rownames(expression_matrix)]
  
  cat("Oxygen depletion genes in expression matrix:", length(oxygen_geneset), "\n")
  cat("cGAS-STING genes in expression matrix:", length(cgas_geneset), "\n")
  
  if (length(oxygen_geneset) >= 5 && length(cgas_geneset) >= 5) {
    
    # Rank the expression matrix (singscore requires ranked data)
    # Higher expression = higher rank
    ranked_expr <- singscore::rankGenes(expression_matrix)
    
    cat("Expression matrix ranked for singscore\n")
    
    # Compute singscores for each pathway
    # simpleScore returns scores for each sample
    oxygen_scores <- singscore::simpleScore(
      rankData = ranked_expr,
      upSet = oxygen_geneset,
      knownDirection = TRUE
    )
    
    cgas_scores <- singscore::simpleScore(
      rankData = ranked_expr,
      upSet = cgas_geneset,
      knownDirection = TRUE
    )
    
    cat("Singscores computed\n")
    
    # Extract scores and create data frame
    singscore_df <- data.frame(
      sample = colnames(expression_matrix),
      group = gsub("_.*", "", colnames(expression_matrix)),
      Oxygen_Depletion_Score = oxygen_scores$TotalScore,
      Oxygen_Depletion_Dispersion = oxygen_scores$TotalDispersion,
      cGAS_STING_Score = cgas_scores$TotalScore,
      cGAS_STING_Dispersion = cgas_scores$TotalDispersion,
      stringsAsFactors = FALSE
    )
    
    # Save singscore results
    write.csv(singscore_df,
              file.path(rnaseq_base_dir, "pathway_singscores.csv"),
              row.names = FALSE)
    cat("Singscores saved to pathway_singscores.csv\n")
    
    # Summary by group
    cat("\n--- Singscore Summary by Group ---\n")
    singscore_summary <- singscore_df %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        mean_oxygen = mean(Oxygen_Depletion_Score),
        sd_oxygen = sd(Oxygen_Depletion_Score),
        mean_cgas = mean(cGAS_STING_Score),
        sd_cgas = sd(cGAS_STING_Score),
        n = n(),
        .groups = "drop"
      )
    print(singscore_summary)
    
    # Statistical testing: compare each radiation group to NIR
    cat("\n--- Statistical Testing (vs NIR) ---\n")
    nir_oxygen <- singscore_df$Oxygen_Depletion_Score[singscore_df$group == "NIR"]
    nir_cgas <- singscore_df$cGAS_STING_Score[singscore_df$group == "NIR"]
    
    rad_groups <- unique(singscore_df$group[singscore_df$group != "NIR"])
    
    singscore_stats <- data.frame()
    
    for (grp in rad_groups) {
      grp_oxygen <- singscore_df$Oxygen_Depletion_Score[singscore_df$group == grp]
      grp_cgas <- singscore_df$cGAS_STING_Score[singscore_df$group == grp]
      
      # t-test vs NIR
      oxygen_test <- t.test(grp_oxygen, nir_oxygen)
      cgas_test <- t.test(grp_cgas, nir_cgas)
      
      cat(sprintf("%s vs NIR:\n", grp))
      cat(sprintf("  Oxygen Depletion: diff = %.4f, p = %.4f %s\n",
                  mean(grp_oxygen) - mean(nir_oxygen),
                  oxygen_test$p.value,
                  ifelse(oxygen_test$p.value < 0.05, "*", "")))
      cat(sprintf("  cGAS-STING: diff = %.4f, p = %.4f %s\n",
                  mean(grp_cgas) - mean(nir_cgas),
                  cgas_test$p.value,
                  ifelse(cgas_test$p.value < 0.05, "*", "")))
      
      singscore_stats <- rbind(singscore_stats, data.frame(
        comparison = paste0(grp, "_vs_NIR"),
        pathway = c("Oxygen_Depletion", "cGAS_STING"),
        mean_diff = c(mean(grp_oxygen) - mean(nir_oxygen), mean(grp_cgas) - mean(nir_cgas)),
        p_value = c(oxygen_test$p.value, cgas_test$p.value),
        stringsAsFactors = FALSE
      ))
    }
    
    # Add significance labels
    singscore_stats <- singscore_stats %>%
      dplyr::mutate(
        sig_label = case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01 ~ "**",
          p_value < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      )
    
    # Save statistics
    write.csv(singscore_stats,
              file.path(rnaseq_base_dir, "pathway_singscore_statistics.csv"),
              row.names = FALSE)
    
    # Visualizations
    pdf(file.path(rnaseq_base_dir, "pathway_singscore_analysis.pdf"), width = 12, height = 10)
    
    # 1. Bar plot of singscores by group with significance
    singscore_long <- singscore_df %>%
      tidyr::pivot_longer(
        cols = c(Oxygen_Depletion_Score, cGAS_STING_Score),
        names_to = "pathway",
        values_to = "score"
      ) %>%
      dplyr::mutate(pathway = gsub("_Score", "", pathway))
    
    # Calculate mean and SE for plotting
    singscore_plot_data <- singscore_long %>%
      dplyr::group_by(group, pathway) %>%
      dplyr::summarise(
        mean_score = mean(score),
        se_score = sd(score) / sqrt(n()),
        .groups = "drop"
      )
    
    # Add significance from comparison to NIR
    singscore_plot_data <- singscore_plot_data %>%
      dplyr::left_join(
        singscore_stats %>%
          dplyr::mutate(group = gsub("_vs_NIR", "", comparison)) %>%
          dplyr::select(group, pathway, sig_label),
        by = c("group", "pathway")
      ) %>%
      dplyr::mutate(
        sig_label = ifelse(is.na(sig_label), "", sig_label),
        label_y = mean_score + se_score + 0.01
      )
    
    p_sing1 <- ggplot(singscore_plot_data, aes(x = group, y = mean_score, fill = pathway)) +
      geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.8) +
      geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score),
                    position = position_dodge(0.9), width = 0.25) +
      geom_text(aes(y = label_y, label = sig_label),
                position = position_dodge(0.9), size = 4, vjust = 0) +
      labs(title = "Pathway Activation Singscores (Mean ± SE)",
           subtitle = "Significance vs NIR: * p<0.05, ** p<0.01, *** p<0.001",
           x = "Group",
           y = "Singscore",
           fill = "Pathway") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_manual(values = c("Oxygen_Depletion" = "#E69F00", "cGAS_STING" = "#56B4E9"))
    print(p_sing1)
    
    # 2. Faceted bar plot
    p_sing2 <- ggplot(singscore_plot_data, aes(x = group, y = mean_score, fill = pathway)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score),
                    width = 0.25) +
      geom_text(aes(y = label_y, label = sig_label), size = 4, vjust = 0) +
      facet_wrap(~ pathway, ncol = 1, scales = "free_y") +
      labs(title = "Pathway Singscores by Group (Mean ± SE)",
           subtitle = "Significance vs NIR: * p<0.05, ** p<0.01, *** p<0.001",
           x = "Group",
           y = "Singscore") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      scale_fill_manual(values = c("Oxygen_Depletion" = "#E69F00", "cGAS_STING" = "#56B4E9"))
    print(p_sing2)
    
    # 3. Individual sample scores with jitter
    p_sing3 <- ggplot(singscore_long, aes(x = group, y = score, fill = pathway)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(0.8)) +
      geom_jitter(aes(color = pathway), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                  alpha = 0.6, size = 2) +
      facet_wrap(~ pathway, ncol = 1, scales = "free_y") +
      labs(title = "Individual Sample Singscores",
           x = "Group",
           y = "Singscore") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      scale_fill_manual(values = c("Oxygen_Depletion" = "#E69F00", "cGAS_STING" = "#56B4E9")) +
      scale_color_manual(values = c("Oxygen_Depletion" = "#E69F00", "cGAS_STING" = "#56B4E9"))
    print(p_sing3)
    
    # 4. Scatter plot: Oxygen vs cGAS-STING scores
    p_sing4 <- ggplot(singscore_df, aes(x = Oxygen_Depletion_Score, y = cGAS_STING_Score, color = group)) +
      geom_point(size = 3, alpha = 0.7) +
      labs(title = "Oxygen Depletion vs cGAS-STING Pathway Scores",
           x = "Oxygen Depletion Singscore",
           y = "cGAS-STING Singscore",
           color = "Group") +
      theme_bw() +
      geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "gray")
    print(p_sing4)
    
    # 5. Heatmap of singscores
    singscore_matrix <- as.matrix(singscore_df[, c("Oxygen_Depletion_Score", "cGAS_STING_Score")])
    rownames(singscore_matrix) <- singscore_df$sample
    colnames(singscore_matrix) <- c("Oxygen Depletion", "cGAS-STING")
    
    sample_annot <- data.frame(
      Group = singscore_df$group,
      row.names = singscore_df$sample
    )
    
    pheatmap::pheatmap(
      t(singscore_matrix),
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      annotation_col = sample_annot,
      main = "Pathway Singscores Across Samples",
      color = colorRampPalette(c("blue", "white", "red"))(100),
      cellwidth = 12,
      cellheight = 25
    )
    
    dev.off()
    
    cat("\nSingscore plots saved to pathway_singscore_analysis.pdf\n")
    
  } else {
    cat("WARNING: Not enough pathway genes found in expression matrix for singscore analysis\n")
  }
  
  # ==============================================================================
  # EMT (Epithelial-Mesenchymal Transition) Analysis
  # ==============================================================================
  
  cat("\n=== Computing EMT Score for Each Sample ===\n")
  
  # EMT gene signatures based on literature (Tan et al., Groger et al., Hallmarks)
  # Epithelial markers - genes downregulated during EMT
  epithelial_genes <- c(
    # Tight junction proteins
    "CDH1", "EPCAM", "OCLN", "TJP1", "TJP2", "TJP3", "CLDN1", "CLDN3", "CLDN4", "CLDN7",
    # Cytokeratins
    "KRT8", "KRT18", "KRT19", "KRT7", "KRT5", "KRT14",
    # Desmosomal proteins
    "DSP", "DSG2", "DSC2", "PKP1", "PKP3", "JUP",
    # Other epithelial markers
    "MUC1", "GRHL2", "ESRP1", "ESRP2", "OVOL1", "OVOL2",
    # miR-200 family targets (epithelial)
    "CRB3", "ST14", "RAB25"
  )
  
  # Mesenchymal markers - genes upregulated during EMT
  mesenchymal_genes <- c(
    # Core EMT transcription factors
    "SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2",
    # Mesenchymal markers
    "VIM", "CDH2", "FN1", "ACTA2", "S100A4", "FSP1",
    # ECM components
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2",
    "MMP2", "MMP3", "MMP9", "MMP14",
    # EMT-associated signaling
    "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2",
    "WNT5A", "WNT5B", "CTNNB1",
    # Other mesenchymal markers
    "SPARC", "SERPINE1", "ITGB1", "ITGA5", "ITGAV",
    "THY1", "PDGFRB", "FAP", "POSTN", "LOX", "LOXL2",
    # Stemness markers associated with EMT
    "SOX2", "NANOG", "POU5F1", "CD44", "ALDH1A1"
  )
  
  cat("Epithelial signature genes:", length(epithelial_genes), "\n")
  cat("Mesenchymal signature genes:", length(mesenchymal_genes), "\n")
  
  # Filter to genes present in expression matrix
  epi_geneset <- epithelial_genes[epithelial_genes %in% rownames(expression_matrix)]
  mes_geneset <- mesenchymal_genes[mesenchymal_genes %in% rownames(expression_matrix)]
  
  cat("Epithelial genes in expression matrix:", length(epi_geneset), "\n")
  cat("Mesenchymal genes in expression matrix:", length(mes_geneset), "\n")
  
  if (length(epi_geneset) >= 5 && length(mes_geneset) >= 5) {
    
    # Compute EMT score using singscore
    # For EMT: upSet = mesenchymal genes, downSet = epithelial genes
    # Higher score = more mesenchymal (EMT activated)
    # Lower score = more epithelial
    
    emt_scores <- singscore::simpleScore(
      rankData = ranked_expr,
      upSet = mes_geneset,
      downSet = epi_geneset,
      knownDirection = TRUE
    )
    
    # Also compute individual pathway scores
    epi_scores <- singscore::simpleScore(
      rankData = ranked_expr,
      upSet = epi_geneset,
      knownDirection = TRUE
    )
    
    mes_scores <- singscore::simpleScore(
      rankData = ranked_expr,
      upSet = mes_geneset,
      knownDirection = TRUE
    )
    
    cat("EMT singscores computed\n")
    
    # Create EMT data frame
    emt_df <- data.frame(
      sample = colnames(expression_matrix),
      group = gsub("_.*", "", colnames(expression_matrix)),
      EMT_Score = emt_scores$TotalScore,
      EMT_Dispersion = emt_scores$TotalDispersion,
      Epithelial_Score = epi_scores$TotalScore,
      Mesenchymal_Score = mes_scores$TotalScore,
      stringsAsFactors = FALSE
    )
    
    # Save EMT scores
    write.csv(emt_df,
              file.path(rnaseq_base_dir, "EMT_singscores.csv"),
              row.names = FALSE)
    cat("EMT scores saved to EMT_singscores.csv\n")
    
    # Summary by group
    cat("\n--- EMT Score Summary by Group ---\n")
    emt_summary <- emt_df %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        mean_EMT = mean(EMT_Score),
        sd_EMT = sd(EMT_Score),
        mean_Epi = mean(Epithelial_Score),
        mean_Mes = mean(Mesenchymal_Score),
        n = n(),
        .groups = "drop"
      )
    print(emt_summary)
    
    # Statistical testing: compare each radiation group to NIR
    cat("\n--- EMT Statistical Testing (vs NIR) ---\n")
    nir_emt <- emt_df$EMT_Score[emt_df$group == "NIR"]
    nir_epi <- emt_df$Epithelial_Score[emt_df$group == "NIR"]
    nir_mes <- emt_df$Mesenchymal_Score[emt_df$group == "NIR"]
    
    rad_groups <- unique(emt_df$group[emt_df$group != "NIR"])
    
    emt_stats <- data.frame()
    
    for (grp in rad_groups) {
      grp_emt <- emt_df$EMT_Score[emt_df$group == grp]
      grp_epi <- emt_df$Epithelial_Score[emt_df$group == grp]
      grp_mes <- emt_df$Mesenchymal_Score[emt_df$group == grp]
      
      # t-test vs NIR
      emt_test <- t.test(grp_emt, nir_emt)
      epi_test <- t.test(grp_epi, nir_epi)
      mes_test <- t.test(grp_mes, nir_mes)
      
      cat(sprintf("%s vs NIR:\n", grp))
      cat(sprintf("  EMT Score: diff = %.4f, p = %.4f %s\n",
                  mean(grp_emt) - mean(nir_emt),
                  emt_test$p.value,
                  ifelse(emt_test$p.value < 0.05, "*", "")))
      cat(sprintf("  Epithelial: diff = %.4f, p = %.4f %s\n",
                  mean(grp_epi) - mean(nir_epi),
                  epi_test$p.value,
                  ifelse(epi_test$p.value < 0.05, "*", "")))
      cat(sprintf("  Mesenchymal: diff = %.4f, p = %.4f %s\n",
                  mean(grp_mes) - mean(nir_mes),
                  mes_test$p.value,
                  ifelse(mes_test$p.value < 0.05, "*", "")))
      
      emt_stats <- rbind(emt_stats, data.frame(
        comparison = paste0(grp, "_vs_NIR"),
        score_type = c("EMT", "Epithelial", "Mesenchymal"),
        mean_diff = c(mean(grp_emt) - mean(nir_emt), 
                      mean(grp_epi) - mean(nir_epi),
                      mean(grp_mes) - mean(nir_mes)),
        p_value = c(emt_test$p.value, epi_test$p.value, mes_test$p.value),
        stringsAsFactors = FALSE
      ))
    }
    
    # Add significance labels
    emt_stats <- emt_stats %>%
      dplyr::mutate(
        sig_label = case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01 ~ "**",
          p_value < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      )
    
    # Save statistics
    write.csv(emt_stats,
              file.path(rnaseq_base_dir, "EMT_singscore_statistics.csv"),
              row.names = FALSE)
    
    # Visualizations
    pdf(file.path(rnaseq_base_dir, "EMT_analysis.pdf"), width = 12, height = 10)
    
    # 1. EMT Score bar plot by group
    emt_plot_data <- emt_df %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        mean_score = mean(EMT_Score),
        se_score = sd(EMT_Score) / sqrt(n()),
        .groups = "drop"
      ) %>%
      dplyr::left_join(
        emt_stats %>%
          dplyr::filter(score_type == "EMT") %>%
          dplyr::mutate(group = gsub("_vs_NIR", "", comparison)) %>%
          dplyr::select(group, sig_label),
        by = "group"
      ) %>%
      dplyr::mutate(
        sig_label = ifelse(is.na(sig_label), "", sig_label),
        label_y = mean_score + se_score + 0.01
      )
    
    p_emt1 <- ggplot(emt_plot_data, aes(x = group, y = mean_score, fill = group)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score),
                    width = 0.25) +
      geom_text(aes(y = label_y, label = sig_label), size = 5, vjust = 0) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "EMT Score by Group (Mean ± SE)",
           subtitle = "Positive = Mesenchymal, Negative = Epithelial | * p<0.05 vs NIR",
           x = "Group",
           y = "EMT Singscore") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    print(p_emt1)
    
    # 2. Epithelial vs Mesenchymal scores
    emt_long <- emt_df %>%
      tidyr::pivot_longer(
        cols = c(Epithelial_Score, Mesenchymal_Score),
        names_to = "score_type",
        values_to = "score"
      ) %>%
      dplyr::mutate(score_type = gsub("_Score", "", score_type))
    
    emt_long_summary <- emt_long %>%
      dplyr::group_by(group, score_type) %>%
      dplyr::summarise(
        mean_score = mean(score),
        se_score = sd(score) / sqrt(n()),
        .groups = "drop"
      ) %>%
      dplyr::left_join(
        emt_stats %>%
          dplyr::mutate(group = gsub("_vs_NIR", "", comparison)) %>%
          dplyr::select(group, score_type, sig_label),
        by = c("group", "score_type")
      ) %>%
      dplyr::mutate(
        sig_label = ifelse(is.na(sig_label), "", sig_label),
        label_y = mean_score + se_score + 0.01
      )
    
    p_emt2 <- ggplot(emt_long_summary, aes(x = group, y = mean_score, fill = score_type)) +
      geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.8) +
      geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score),
                    position = position_dodge(0.9), width = 0.25) +
      geom_text(aes(y = label_y, label = sig_label),
                position = position_dodge(0.9), size = 4, vjust = 0) +
      labs(title = "Epithelial vs Mesenchymal Scores (Mean ± SE)",
           subtitle = "* p<0.05, ** p<0.01, *** p<0.001 vs NIR",
           x = "Group",
           y = "Singscore",
           fill = "Score Type") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_manual(values = c("Epithelial" = "#2E86AB", "Mesenchymal" = "#A23B72"))
    print(p_emt2)
    
    # 3. Faceted bar plot
    p_emt3 <- ggplot(emt_long_summary, aes(x = group, y = mean_score, fill = score_type)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score),
                    width = 0.25) +
      geom_text(aes(y = label_y, label = sig_label), size = 4, vjust = 0) +
      facet_wrap(~ score_type, ncol = 1, scales = "free_y") +
      labs(title = "EMT Component Scores by Group",
           subtitle = "* p<0.05, ** p<0.01, *** p<0.001 vs NIR",
           x = "Group",
           y = "Singscore") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      scale_fill_manual(values = c("Epithelial" = "#2E86AB", "Mesenchymal" = "#A23B72"))
    print(p_emt3)
    
    # 4. Individual sample scatter: Epithelial vs Mesenchymal
    p_emt4 <- ggplot(emt_df, aes(x = Epithelial_Score, y = Mesenchymal_Score, color = group)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
      labs(title = "Epithelial vs Mesenchymal Scores",
           subtitle = "Points above diagonal = more mesenchymal",
           x = "Epithelial Singscore",
           y = "Mesenchymal Singscore",
           color = "Group") +
      theme_bw()
    print(p_emt4)
    
    # 5. Boxplot with jitter for EMT score
    p_emt5 <- ggplot(emt_df, aes(x = group, y = EMT_Score, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "EMT Score Distribution by Group",
           subtitle = "Positive = Mesenchymal shift, Negative = Epithelial",
           x = "Group",
           y = "EMT Singscore") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    print(p_emt5)
    
    # 6. Heatmap of key EMT genes
    # Select key EMT markers present in data
    key_epi <- c("CDH1", "EPCAM", "KRT8", "KRT18", "KRT19", "CLDN4", "TJP1", "GRHL2", "ESRP1")
    key_mes <- c("VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "MMP2")
    key_emt_genes <- c(key_epi, key_mes)
    key_emt_present <- key_emt_genes[key_emt_genes %in% rownames(expression_matrix)]
    
    if (length(key_emt_present) >= 5) {
      emt_expr <- expression_matrix[key_emt_present, , drop = FALSE]
      
      # Z-score normalize
      emt_expr_z <- t(scale(t(emt_expr)))
      
      # Annotation
      gene_annot <- data.frame(
        Type = ifelse(rownames(emt_expr_z) %in% key_epi, "Epithelial", "Mesenchymal"),
        row.names = rownames(emt_expr_z)
      )
      
      sample_annot <- data.frame(
        Group = gsub("_.*", "", colnames(emt_expr_z)),
        row.names = colnames(emt_expr_z)
      )
      
      pheatmap::pheatmap(
        emt_expr_z,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        annotation_row = gene_annot,
        annotation_col = sample_annot,
        main = "Key EMT Marker Expression (Z-score)",
        color = colorRampPalette(c("blue", "white", "red"))(100),
        breaks = seq(-3, 3, length.out = 101),
        fontsize_row = 10
      )
    }
    
    dev.off()
    
    cat("\nEMT analysis plots saved to EMT_analysis.pdf\n")
    
  } else {
    cat("WARNING: Not enough EMT genes found in expression matrix\n")
  }
  
  cat("\n=== All Pathway Analyses Completed ===\n")
  cat("\nInterpretation:\n")
  cat("- Positive mean stat: Pathway genes are upregulated in irradiated vs NIR\n")
  cat("- Negative mean stat: Pathway genes are downregulated in irradiated vs NIR\n")
  cat("- Higher singscore: Greater pathway activation in that sample\n")
  cat("- EMT Score: Positive = mesenchymal shift, Negative = epithelial\n")
  cat("- * indicates statistically significant difference (p < 0.05)\n")
  cat("- Look for dose-dependent and time-dependent patterns\n\n")
  
} else {
  cat("ERROR: Comparison directory not found:", comp_dir, "\n")
}

