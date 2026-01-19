library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(corrplot)

wd <- "C:/PROJECTS/Shane/Harding_250611/cgmaptools"
output_path <- "C:/PROJECTS/Shane/Harding_250611/cgmaptools/Rplots"
if (!dir.exists(output_path)) {
  dir.create(output_path)
}
setwd(wd)


# Processing mstat files
standard_mstat_files <- list.files(path = wd, pattern = "mstat.tab", full.names = FALSE)
control_mstat_files <- list.files(path = wd, pattern = "mstat_control.tab", full.names = FALSE)

# Collect control stats by chromosome

control_stats <- lapply(control_mstat_files, function(afile){
    shortName <- gsub("\\_mstat_control.tab", "", afile, fixed = FALSE)
    df <- read.delim(afile, header = TRUE, skip = 0)

    chrom_stats <- df %>%
        dplyr::filter(MethStat == "mean_mC_byChr")

    # Convert to long format
    chrom_stats_long <- chrom_stats %>%
        tidyr::pivot_longer(cols = -c(MethStat, context), names_to = "chromosome", values_to = "methylation") %>%
        dplyr::mutate(MethStat = shortName) %>%
        dplyr::rename(Sample = MethStat, Chromosome = context, Context = chromosome, Methylation = methylation)
    
    return(chrom_stats_long)
  
})

control_df <- do.call("rbind", control_stats)

# Collect standard stats by chromosome

standard_stats <- lapply(standard_mstat_files, function(afile){
    shortName <- gsub("\\_mstat.tab", "", afile, fixed = FALSE)
    df <- read.delim(afile, header = TRUE, skip = 0)

    chrom_stats <- df %>%
        dplyr::filter(MethStat == "mean_mC_byChr", context != "chrY")

    # Convert to long format
    chrom_stats_long <- chrom_stats %>%
        tidyr::pivot_longer(cols = -c(MethStat, context), names_to = "chromosome", values_to = "methylation") %>%
        dplyr::mutate(MethStat = shortName) %>%
        dplyr::rename(Sample = MethStat, Chromosome = context, Context = chromosome, Methylation = methylation)
  
    return(chrom_stats_long)
  
})

standard_df <- do.call("rbind", standard_stats)

# Initialize storage for heatmap data
heatmap_data <- list()

for(context in c("C", "CG", "CHG", "CHH")){
    # transform control_df to wide format and plot heatmap  
    control_df_wide <- control_df %>%
        tidyr::pivot_wider(names_from = Sample, values_from = Methylation) %>%
        dplyr::filter(Context == context)
    control_matrix <- as.matrix(control_df_wide[, -c(1,2)])
    rownames(control_matrix) <- control_df_wide$Chromosome
    
    # transform standard_df to wide format and plot heatmap
    standard_df_wide <- standard_df %>%
        tidyr::pivot_wider(names_from = Sample, values_from = Methylation) %>%
        dplyr::filter(Context == context)
    standard_matrix <- as.matrix(standard_df_wide[, -c(1,2)])
    rownames(standard_matrix) <- standard_df_wide$Chromosome

    combined_matrix <- rbind(control_matrix, standard_matrix)

    # scale combined matrix
    scaled_matrix <- scale(t(combined_matrix))
    scaled_matrix <- t(scaled_matrix)

    # Matrix-based normalization using internal controls
    
    # Check available control rows
    available_controls <- rownames(control_matrix)
    cat("Available controls in", context, "context:", 
        paste(available_controls, collapse = ", "), "\n")
    
    # Find lambda and pUC19 controls (flexible naming)
    lambda_row <- grep("lambda|phage", available_controls, ignore.case = TRUE)
    puc19_row <- grep("plasmid|puc19", available_controls, ignore.case = TRUE)
    
    if(length(lambda_row) == 0 | length(puc19_row) == 0) {
        cat("Warning: Missing controls in", context, "context. Skipping normalization.\n")
        next
    }
    
    # Extract control values
    lambda_values <- control_matrix[lambda_row[1], ]
    puc19_values <- control_matrix[puc19_row[1], ]
    
    # Context-specific normalization strategy
    if(context %in% c("CG", "C")) {
        # For CG: pUC19 should be highly methylated (~95%), lambda unmethylated (~0%)
        cat(paste(context, " context: Using dual control normalization\n", sep=""))
        cat("Lambda (background):", round(mean(lambda_values, na.rm = TRUE), 3), "\n")
        cat("pUC19 (positive control):", round(mean(puc19_values, na.rm = TRUE), 3), "\n")
        
        # Background subtraction using lambda
        background_methylation <- lambda_values
        
        # Scaling factor using pUC19 efficiency
        expected_puc19 <- 0.24  # Expected 24% methylation for C context (pUC19 efficiency)
        if(context == "CG"){
            expected_puc19 <- 0.95 # Expected 95% methylation for CG context (pUC19 efficiency)
        } 
        observed_puc19 <- puc19_values - background_methylation  # Background-not-corrected pUC19
        scaling_factors <- (expected_puc19 / observed_puc19) # Square the scaling factor
        
        # Apply dual normalization: background subtraction + scaling
        normalized_combined_matrix <- sweep(combined_matrix, 2, background_methylation, "-")
        normalized_combined_matrix <- sweep(normalized_combined_matrix, 2, scaling_factors, "*")

        # Make sure all values of the plasmid row are equal to 0.95, to avoid unexpected zscores after scaling
        normalized_combined_matrix[puc19_row[1], ] <- mean(normalized_combined_matrix[puc19_row[1], ])

        names(scaling_factors) <- colnames(combined_matrix)
        sink(paste0(output_path, "/scaling_factors_", context, ".txt"))
        print(scaling_factors)
        sink()  
        
        
    } else {
        # For CHG/CHH: Both controls should be unmethylated (~0%)
        cat("Non-CG context (", context, "): Using background subtraction\n")
        cat("Lambda (background):", round(mean(lambda_values, na.rm = TRUE), 3), "\n")
        cat("pUC19 (also background):", round(mean(puc19_values, na.rm = TRUE), 3), "\n")
        
        # Use average of both controls as background (both should be ~0%)
        background_methylation <- (lambda_values + puc19_values) / 2
        
        # Simple background subtraction
        normalized_combined_matrix <- sweep(combined_matrix, 2, background_methylation, "-")
    }
    
    # Quality control: Print normalization summary
    cat("Normalization summary for", context, ":\n")
    cat("  Original range:", round(range(combined_matrix, na.rm = TRUE), 3), "\n")
    cat("  Normalized range:", round(range(normalized_combined_matrix, na.rm = TRUE), 3), "\n")
    cat("  Background correction:", round(mean(background_methylation, na.rm = TRUE), 3), "\n")

    # Prepare scaled normalized matrix
    scaled_normalized_combined_matrix <- scale(t(normalized_combined_matrix))
    scaled_normalized_combined_matrix <- t(scaled_normalized_combined_matrix)
    if(context == "CG") {
        print(normalized_combined_matrix[1:5, ])
        print(scaled_normalized_combined_matrix[1:5, ])
    }
    
    # Store matrices for combined plotting
    heatmap_data[[context]] <- list(
        original = combined_matrix,
        scaled = scaled_matrix,
        normalized = normalized_combined_matrix,
        normalized_scaled = scaled_normalized_combined_matrix
    )
}

# Load required library for grid arrangement
if (!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    library(gridExtra)
}

# Generate combined heatmaps for each context using grid graphics
for(context in names(heatmap_data)) {
    
    cat("Generating aggregated heatmap for", context, "context...\n")
    
    # Create pheatmap objects (without displaying them)
    p1 <- pheatmap(heatmap_data[[context]]$original, 
                  cluster_rows = FALSE, cluster_cols = TRUE, 
                  show_rownames = TRUE, show_colnames = TRUE, 
                  annotation_row = NULL, annotation_col = NULL, 
                  display_numbers = TRUE, number_format = "%.3f",
                  main = paste("Original Combined (", context, ")", sep = ""),
                  silent = TRUE)

    p2 <- pheatmap(heatmap_data[[context]]$scaled, 
                  cluster_rows = FALSE, cluster_cols = TRUE, 
                  show_rownames = TRUE, show_colnames = TRUE, 
                  annotation_row = NULL, annotation_col = NULL, 
                  display_numbers = TRUE, number_format = "%.3f",
                  main = paste("Scaled Combined (", context, ")", sep = ""),
                  silent = TRUE)
    
    p3 <- pheatmap(heatmap_data[[context]]$normalized, 
                  cluster_rows = FALSE, cluster_cols = TRUE, 
                  show_rownames = TRUE, show_colnames = TRUE, 
                  annotation_row = NULL, annotation_col = NULL, 
                  display_numbers = TRUE, number_format = "%.3f",
                  main = paste("Normalized Combined (", context, ")", sep = ""),
                  silent = TRUE)
    
    p4 <- pheatmap(heatmap_data[[context]]$normalized_scaled, 
                  cluster_rows = FALSE, cluster_cols = TRUE, 
                  show_rownames = TRUE, show_colnames = TRUE, 
                  annotation_row = NULL, annotation_col = NULL, 
                  display_numbers = TRUE, number_format = "%.3f",
                  main = paste("Normalized + Scaled (", context, ")", sep = ""),
                  silent = TRUE)
    
    # Create a single PDF with all heatmaps arranged in 2x2 grid
    pdf(paste0(output_path, "/aggregated_heatmaps_", context, ".pdf"), 
        width = 20, height = 16)
    
    # Arrange the four heatmaps in a 2x2 grid
    grid.arrange(p1[[4]], p2[[4]], p3[[4]], p4[[4]], 
                ncol = 2, nrow = 2,
                top = paste("EM-seq Methylation Analysis - ", context, " Context", sep = ""))
    
    dev.off()
    
    cat("✓ Generated aggregated heatmap for", context, "context\n")
}

# Generate combined correlation plots for all contexts in a single PDF
cat("Generating combined correlation plots...\n")
pdf(paste0(output_path, "/chromosomal_correlation_plots.pdf"), 
    width = 20, height = 16)

# Set up 2x2 layout for the four contexts
par(mfrow = c(2, 2))

for(context in names(heatmap_data)) {
    # Plot correlation for original data
    correlation_matrix <- cor(t(heatmap_data[[context]]$original))
    corrplot(correlation_matrix, 
             method = "circle", 
             type = "upper",
             addCoef.col = "white",
             number.cex = 0.6,
             number.digits = 2,
             title = paste("Pearson Correlation: ", context, " Context", sep = ""),
             mar = c(0, 0, 2, 0))
}

dev.off()
cat("✓ Generated combined correlation plots\n")
