rm(list = ls())

library(openxlsx)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

# If each element of enrichment_list is a data.frame with columns:
# - "group_id"
# - value column (e.g., "enrichment")
# Set value_col accordingly.
build_enrichment_matrix <- function(enrichment_list, value_col = "enrichment") {
  if (length(enrichment_list) == 0) stop("enrichment_list is empty")

  # Ensure samples are named
  if (is.null(names(enrichment_list)) || any(names(enrichment_list) == "")) {
    stop("enrichment_list must be a named list of samples")
  }

  # Detect structure: data.frame with family_id + value, or named numeric vector
  is_df <- vapply(enrichment_list, function(x) is.data.frame(x), logical(1))
  is_vec <- vapply(enrichment_list, function(x) is.numeric(x) && !is.null(names(x)), logical(1))

  if (!all(is_df | is_vec)) {
    stop("Each element of enrichment_list must be either a data.frame(family_id, value) or a named numeric vector")
  }

  # Collect all family IDs
  families <- unique(unlist(
    lapply(names(enrichment_list), function(samp) {
      x <- enrichment_list[[samp]]
      if (is.data.frame(x)) {
        if (!("group_id" %in% names(x))) stop("Missing 'group_id' column in some data.frames")
        if (!(value_col %in% names(x))) stop(paste0("Missing '", value_col, "' column in some data.frames"))
        x$group_id
      } else {
        names(x)
      }
    })
  ))

  # Initialize matrix with NA
  mat <- matrix(NA_real_, nrow = length(families), ncol = length(enrichment_list),
                dimnames = list(families, names(enrichment_list)))

  # Fill matrix
  for (samp in names(enrichment_list)) {
    x <- enrichment_list[[samp]]
    if (is.data.frame(x)) {
      vals <- setNames(x[[value_col]], x$group_id)
    } else {
      vals <- x
    }
    common <- intersect(names(vals), rownames(mat))
    mat[common, samp] <- vals[common]
  }

  # Return as data.frame with group_id as rownames
  as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
}


# Main analysis starts here ##########################################################
project_dir <- "C:/PROJECTS/Shane/Harding_250611/wo_chrY"
subfolder <- "cutoff0.05"
# Collect all heatmaps
all_heatmaps <- list()

for(dmr_method in c("dmrseq", "TE_targeted", "DSS")){
    base_dir <- file.path(project_dir, dmr_method)
    dmr_folder <- "DMRs"

    if(dmr_method == "TE_targeted") dmr_folder <- "Targeted"
    if(dmr_method == "dmrseq"){
        base_dir <- file.path(base_dir, subfolder)
    } 
    sample_dirs <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)

    if(dmr_method == "dmrseq"){
        sample_names <- gsub("dmrseq_", "", sample_dirs)
        sample_names <- gsub("_", "_vs_", sample_names)
        
    }else{
        sample_names <- sample_dirs
    }

    names(sample_dirs) <- sample_names

    hyper_TE_list <- list()
    hypo_TE_list <- list()

    # keep sample names in specific order
    sample_names <- c("IR2Gy6d_vs_NIR", "IR2Gy24h_vs_NIR", "IR10Gy6d_vs_NIR", "IR10Gy24h_vs_NIR", "IR2Gy6d_vs_IR2Gy24h", "IR10Gy6d_vs_IR10Gy24h", "IR10Gy6d_vs_IR2Gy6d", "IR10Gy24h_vs_IR2Gy24h")

    for(sample_name in sample_names){
        message("Processing ", dmr_method, " ", sample_name, " ###################################################################")

        # import DMRs and TE annotations
        working_dir <- file.path(base_dir, sample_dirs[sample_name], dmr_folder)

        hyper_enrichment <- read.xlsx(paste0(working_dir, "/TE_family_id_enrichment_hypermethylated_DMRs.xlsx")) %>%
            dplyr::mutate(enrichment = -log10(padj)) %>%
            dplyr::select(group_id, enrichment)
        hypo_enrichment <- read.xlsx(paste0(working_dir, "/TE_family_id_enrichment_hypomethylated_DMRs.xlsx")) %>%
            dplyr::mutate(enrichment = -log10(padj)) %>%
            dplyr::select(group_id, enrichment)

        hyper_TE_list[[sample_name]] <- hyper_enrichment
        hypo_TE_list[[sample_name]] <- hypo_enrichment
    }

    # create a data frame of hyper enrichments, with family_id as row names and sample names as column names
    hyper_df <- build_enrichment_matrix(hyper_TE_list, value_col = "enrichment")
    hypo_df <- build_enrichment_matrix(hypo_TE_list, value_col = "enrichment")

    head(hyper_df)
    head(hypo_df)
    # set NAs to 0
    hyper_df[is.na(hyper_df)] <- 0
    hypo_df[is.na(hypo_df)] <- 0
    # set Inf to 150
    hyper_df[sapply(hyper_df, is.infinite)] <- 150
    hypo_df[sapply(hypo_df, is.infinite)] <- 150

    # remove rows with all values < -log10(0.05)
    hyper_df <- hyper_df[rowSums(hyper_df > -log10(0.05)) > 0, ]
    hypo_df <- hypo_df[rowSums(hypo_df > -log10(0.05)) > 0, ]
    # plot a heatmap of hyper enrichments with custom color breaks
    # Define color breaks: values 0-3 get minimal color, values >3 get visible colors
    color_breaks <- c(0, 1.3, 3, seq(5, 150, length.out = 7))
    colors <- colorRampPalette(c("white", "lightblue", "blue", "darkblue", "red"))(length(color_breaks)-1)
    
    # Create clean matrices without object names
    hyper_mat <- matrix(as.numeric(as.matrix(hyper_df)), 
                       nrow = nrow(hyper_df), ncol = ncol(hyper_df),
                       dimnames = list(rownames(hyper_df), colnames(hyper_df)))
    hypo_mat <- matrix(as.numeric(as.matrix(hypo_df)), 
                      nrow = nrow(hypo_df), ncol = ncol(hypo_df),
                      dimnames = list(rownames(hypo_df), colnames(hypo_df)))
    
    # Use ComplexHeatmap for better legend control
    # For colorRamp2, we need breaks and colors to have same length
    # Take representative breaks for the color mapping
    color_breaks_for_ramp <- c(0, 1.3, 3, 10, 25, 50, 75, 100, 150)
    col_fun <- colorRamp2(color_breaks_for_ramp, colors)
    
    # Create heatmaps with proper legend titles and collect them
    ht_hyper <- Heatmap(hyper_mat, 
                       name = "-log10(padj)",
                       col = col_fun,
                       cluster_rows = TRUE, 
                       cluster_columns = FALSE,
                       show_row_names = TRUE, 
                       show_column_names = TRUE,
                       column_title = paste0(dmr_method, ": Hypermethylated DMRs"),
                       heatmap_legend_param = list(title = "-log10(padj)"))
    
    ht_hypo <- Heatmap(hypo_mat, 
                      name = "-log10(padj)",
                      col = col_fun,
                      cluster_rows = TRUE, 
                      cluster_columns = FALSE,
                      show_row_names = TRUE, 
                      show_column_names = TRUE,
                      column_title = paste0(dmr_method, ": Hypomethylated DMRs"),
                      heatmap_legend_param = list(title = "-log10(padj)"))
    
    # Add to collection
    all_heatmaps[[paste0(dmr_method, "_hyper")]] <- ht_hyper
    all_heatmaps[[paste0(dmr_method, "_hypo")]] <- ht_hypo
    
    # Also save individual heatmaps for this method
    pdf(file.path(base_dir, paste0("TE_family_heatmap_", dmr_method, ".pdf")))
    draw(ht_hyper)
    draw(ht_hypo)
    dev.off()
    cat("Individual heatmap saved to:", file.path(base_dir, paste0("TE_family_heatmap_", dmr_method, ".pdf")), "\n")
}

# Create combined heatmap plot and save in project directory
if(length(all_heatmaps) > 0) {
    pdf(file.path(project_dir, paste0("TE_family_heatmap_all_methods_", subfolder, ".pdf")), width = 20, height = 20)
    
    # Create separate plots for each method to avoid row mismatch issues
    plot_count <- 0
    
    for(method in c("dmrseq", "TE_targeted", "DSS")) {
        hyper_key <- paste0(method, "_hyper")
        hypo_key <- paste0(method, "_hypo")
        
        if(hyper_key %in% names(all_heatmaps) && hypo_key %in% names(all_heatmaps)) {
            plot_count <- plot_count + 1
            
            # Draw each method's heatmaps separately on the same page
            if(plot_count == 1) {
                # First method - create new page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(3, 2)))
            }
            
            # Draw hypermethylated heatmap
            pushViewport(viewport(layout.pos.row = plot_count, layout.pos.col = 1))
            draw(all_heatmaps[[hyper_key]], newpage = FALSE)
            popViewport()
            
            # Draw hypomethylated heatmap  
            pushViewport(viewport(layout.pos.row = plot_count, layout.pos.col = 2))
            draw(all_heatmaps[[hypo_key]], newpage = FALSE)
            popViewport()
        }
    }
    
    if(plot_count > 0) {
        popViewport()
        # Add overall title
        grid.text("TE Family Enrichment Across DMR Methods", 
                 x = 0.5, y = 0.995, 
                 gp = gpar(fontsize = 16, fontface = "bold"))
    }
    
    dev.off()
    cat("Combined heatmap saved to:", file.path(project_dir, paste0("TE_family_heatmap_all_methods_", subfolder, ".pdf")), "\n")
}

