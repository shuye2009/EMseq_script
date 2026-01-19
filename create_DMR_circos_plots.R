# Circos plots for DMR data using circlize package
# This script creates comprehensive circular genomic visualizations of DMRs

library(circlize)
library(GenomicRanges)
library(openxlsx)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(grid)

# Function to prepare chromosome data for circlize
prepare_chromosome_data <- function() {
    # Human chromosome lengths (hg38/GRCh38)
    chr_lengths <- data.frame(
        chr = paste0("chr", c(1:22, "X", "Y")),
        length = c(248956422, 242193529, 198295559, 190214555, 181538259,
                   170805979, 159345973, 145138636, 138394717, 133797422,
                   135086622, 133275309, 114364328, 107043718, 101991189,
                   90338345, 83257441, 80373285, 58617616, 64444167,
                   46709983, 50818468, 156040895, 57227415),
        stringsAsFactors = FALSE
    )
    return(chr_lengths)
}

# Function to create DMR density bins (similar to DML approach)
create_dmr_density_bins <- function(dmr_data, chr_data, bin_size = 1e6) {
    density_data <- data.frame()
    
    for(chr in chr_data$chr) {
        chr_dmrs <- dmr_data[dmr_data[["DMR_seqnames"]] == chr, ]
        
        if(nrow(chr_dmrs) > 0) {
            chr_length <- chr_data$length[chr_data$chr == chr]
            bins <- seq(0, chr_length, by = bin_size)
            bin_counts <- rep(0, length(bins) - 1)
            
            # Count DMRs in each bin based on their start position
            for(i in seq_len(length(bins) - 1)) {
                bin_counts[i] <- sum(chr_dmrs[["DMR_start"]] >= bins[i] & 
                                   chr_dmrs[["DMR_start"]] < bins[i + 1])
            }
            
            chr_density <- data.frame(
                chr = chr, 
                start = bins[-length(bins)], 
                end = bins[-1],
                count = bin_counts, 
                density = bin_counts / (bin_size / 1e6)  # DMRs per Mb
            )
            density_data <- rbind(density_data, chr_density)
        }
    }
    return(density_data)
}

# Function to create basic circos plot with DMRs
create_dmr_circos <- function(dmr_data, sample_name, output_dir, 
                              show_te = TRUE, show_expression = TRUE) {
    
    cat("Creating Circos plot for", sample_name, "\n")
    
    # Prepare chromosome data
    chr_data <- prepare_chromosome_data()
    
    # Filter DMR data to standard chromosomes
    dmr_filtered <- dmr_data[dmr_data[["DMR_seqnames"]] %in% chr_data$chr, ]
    
    if(nrow(dmr_filtered) == 0) {
        cat("No DMRs found on standard chromosomes for", sample_name, "\n")
        return(NULL)
    }
    
    # Create output filename
    pdf_file <- file.path(output_dir, paste0("circos_DMR_", sample_name, ".pdf"))
    
    # Start PDF device with white background for better color contrast
    pdf(pdf_file, width = 12, height = 12, bg = "white")
    
    # Initialize circos plot
    circos.clear()
    circos.par("start.degree" = 90, "gap.degree" = 2)
    
    # Initialize with chromosome data
    circos.initialize(factors = chr_data$chr, xlim = cbind(rep(0, nrow(chr_data)), chr_data$length))
    
    # Track 1: Chromosome ideogram
    circos.track(factors = chr_data$chr, ylim = c(0, 1), 
                 panel.fun = function(x, y) {
                     circos.text(CELL_META$xcenter, CELL_META$ycenter, 
                                CELL_META$sector.index, cex = 0.8, facing = "clockwise", 
                                niceFacing = TRUE)
                 }, track.height = 0.05, bg.border = NA, bg.col = "grey90")
    
    # Track 2: DMR density bins
    if(nrow(dmr_filtered) > 0) {
        # Create density bins for DMRs
        density_data <- create_dmr_density_bins(dmr_filtered, chr_data, bin_size = 1e6)
        max_density <- max(density_data$density, na.rm = TRUE)
        
        if(max_density > 0) {
            circos.track(factors = chr_data$chr, ylim = c(0, max_density), 
                         panel.fun = function(x, y) {
                             chr <- CELL_META$sector.index
                             chr_density <- density_data[density_data$chr == chr, ]
                             
                             if(nrow(chr_density) > 0) {
                                 for(i in seq_len(nrow(chr_density))) {
                                     if(chr_density$density[i] > 0) {
                                         circos.rect(chr_density$start[i], 0, 
                                                    chr_density$end[i], chr_density$density[i], 
                                                    col = "darkblue", border = NA)
                                     }
                                 }
                             }
                         }, track.height = 0.1, bg.border = "#b1b1b1")
            
            circos.yaxis(side = "left", at = pretty(c(0, max_density)), 
                         labels.cex = 0.6, tick.length = 0.1)
            
            lgd_density <- Legend(at = c("DMR Density"), 
                                legend_gp = gpar(fill = "darkblue"),
                                title = "Track 1: DMR Density",
                                title_gp = gpar(fontsize = 10, fontface = "bold"),
                                labels_gp = gpar(fontsize = 8))
        }
    }
    
    # Track 2: Individual DMRs colored by methylation change
    if("DMR_diff.Methy" %in% colnames(dmr_filtered)) {
        # Create color mapping for methylation changes
        methy_range <- range(dmr_filtered$DMR_diff.Methy, na.rm = TRUE)
        col_fun <- colorRamp2(c(methy_range[1], 0, methy_range[2]), 
                             c("blue", "white", "red"))
        
        circos.track(factors = chr_data$chr, ylim = c(-1, 1), 
                     panel.fun = function(x, y) {
                         chr <- CELL_META$sector.index
                         chr_dmrs <- dmr_filtered[dmr_filtered$DMR_seqnames == chr, ]
                         
                         if(nrow(chr_dmrs) > 0) {
                             for(i in seq_len(nrow(chr_dmrs))) {
                                 start_pos <- chr_dmrs$DMR_start[i]
                                 end_pos <- chr_dmrs$DMR_end[i]
                                 methy_val <- chr_dmrs$DMR_diff.Methy[i]
                                 
                                 if(!is.na(methy_val)) {
                                     y_pos <- ifelse(methy_val > 0, 0.5, -0.5)
                                     circos.rect(start_pos, 0, end_pos, y_pos, 
                                                col = col_fun(methy_val), border = NA)
                                 }
                             }
                         }
                     }, track.height = 0.1, bg.border = "#b1b1b1")
        
        # Add legend for methylation changes
        lgd_methy <- Legend(at = c(methy_range[1], 0, methy_range[2]),
                           col_fun = col_fun, 
                           title = "Track2: Methylation Change",
                           title_gp = gpar(fontsize = 10, fontface = "bold"),
                           labels_gp = gpar(fontsize = 8))
    }
    
    # Track 4: TE family enrichment (if TE data available)
    if(show_te && "TE_family_id" %in% colnames(dmr_filtered)) {
        # Get top TE families
        te_counts <- table(dmr_filtered$TE_family_id)
        top_tes <- names(sort(te_counts, decreasing = TRUE))[seq_len(min(10, length(te_counts)))]
        
        # Create color palette for TE families
        te_colors <- rainbow(length(top_tes))
        names(te_colors) <- top_tes
        
        circos.track(factors = chr_data$chr, ylim = c(0, 1), 
                     panel.fun = function(x, y) {
                         chr <- CELL_META$sector.index
                         chr_dmrs <- dmr_filtered[dmr_filtered$DMR_seqnames == chr, ]
                         
                         if(nrow(chr_dmrs) > 0) {
                             for(i in seq_len(nrow(chr_dmrs))) {
                                 te_family <- chr_dmrs$TE_family_id[i]
                                 if(!is.na(te_family) && te_family %in% top_tes) {
                                     start_pos <- chr_dmrs$DMR_start[i]
                                     end_pos <- chr_dmrs$DMR_end[i]
                                     circos.rect(start_pos, 0, end_pos, 0.8, 
                                                col = te_colors[te_family], border = NA)
                                 }
                             }
                         }
                     }, track.height = 0.1, bg.border = "#b1b1b1")
        
        # Add legend for TE families
        lgd_te <- Legend(at = top_tes, legend_gp = gpar(fill = te_colors),
                        title = "Track3: Top TE Families",
                        title_gp = gpar(fontsize = 10, fontface = "bold"),
                        labels_gp = gpar(fontsize = 8))
    }
    
    # Track 5: Gene expression changes (if available)
    if(show_expression && "log2FoldChange" %in% colnames(dmr_filtered)) {
        # Filter for DMRs with expression data
        expr_dmrs <- dmr_filtered[!is.na(dmr_filtered$log2FoldChange), ]
        
        if(nrow(expr_dmrs) > 0) {
            expr_range <- range(expr_dmrs$log2FoldChange, na.rm = TRUE)
            expr_col_fun <- colorRamp2(c(expr_range[1], 0, expr_range[2]), 
                                      c("green", "white", "orange"))
            
            circos.track(factors = chr_data$chr, ylim = c(-1, 1), 
                         panel.fun = function(x, y) {
                             chr <- CELL_META$sector.index
                             chr_dmrs <- expr_dmrs[expr_dmrs$DMR_seqnames == chr, ]
                             
                             if(nrow(chr_dmrs) > 0) {
                                 for(i in seq_len(nrow(chr_dmrs))) {
                                     start_pos <- chr_dmrs$DMR_start[i]
                                     end_pos <- chr_dmrs$DMR_end[i]
                                     expr_val <- chr_dmrs$log2FoldChange[i]
                                     
                                     if(!is.na(expr_val)) {
                                         y_pos <- ifelse(expr_val > 0, 0.5, -0.5)
                                         circos.rect(start_pos, 0, end_pos, y_pos, 
                                                    col = expr_col_fun(expr_val), border = NA)
                                     }
                                 }
                             }
                         }, track.height = 0.1, bg.border = "#b1b1b1")
            
            # Add legend for expression changes
            lgd_expr <- Legend(at = c(expr_range[1], 0, expr_range[2]),
                              col_fun = expr_col_fun, 
                              title = "Track4: Expression Change (log2FC)",
                              title_gp = gpar(fontsize = 10, fontface = "bold"),
                              labels_gp = gpar(fontsize = 8))
        }
    }
    
    # Add title
    title(paste("DMR circos Plot:", sample_name), cex.main = 1.5)
    
    # Add legends - ensure we're in the right graphics context
    legend_list <- list()
    if(exists("lgd_density")) legend_list <- append(legend_list, list(lgd_density))
    if(exists("lgd_methy")) legend_list <- append(legend_list, list(lgd_methy))
    if(exists("lgd_te")) legend_list <- append(legend_list, list(lgd_te))
    if(exists("lgd_expr")) legend_list <- append(legend_list, list(lgd_expr))
    
    if(length(legend_list) > 0) {
        # Create a new page for legends
        #grid.newpage()
        pushViewport(viewport(width = 1, height = 1))
        draw(packLegend(list = legend_list), x = unit(0.5, "npc"), y = unit(0.5, "npc"))
        popViewport()
    }
    
    # Close PDF device
    dev.off()
    
    cat("Circos plot saved to:", pdf_file, "\n")
    
    # Clear circos parameters
    circos.clear()
    
    return(pdf_file)
}

# Function to create comparative circos plot (multiple samples)
create_comparative_circos <- function(dmr_list, sample_names, output_dir) {
    
    cat("Creating comparative Circos plot for", length(sample_names), "samples\n")
    
    # Prepare chromosome data
    chr_data <- prepare_chromosome_data()
    
    # Create output filename
    pdf_file <- file.path(output_dir, "circos_DMR_comparative.pdf")
    
    # Start PDF device
    pdf(pdf_file, width = 12, height = 12)
    
    # Initialize circos plot
    circos.clear()
    circos.par("start.degree" = 90, "gap.degree" = 2)
    
    # Initialize with chromosome data
    circos.initialize(factors = chr_data$chr, xlim = cbind(rep(0, nrow(chr_data)), chr_data$length))
    
    # Track 1: Chromosome ideogram
    circos.track(factors = chr_data$chr, ylim = c(0, 1), 
                 panel.fun = function(x, y) {
                     circos.text(CELL_META$xcenter, CELL_META$ycenter, 
                                CELL_META$sector.index, cex = 0.8, facing = "clockwise", 
                                niceFacing = TRUE)
                 }, track.height = 0.05, bg.border = NA, bg.col = "grey90")
    
    # Create tracks for each sample
    colors <- brewer.pal(min(length(sample_names), 8), "Set2")
    if(length(sample_names) > 8) {
        colors <- rainbow(length(sample_names))
    }
    names(colors) <- sample_names
    
    for(i in seq_along(dmr_list)) {
        dmr_data <- dmr_list[[i]]
        current_sample <- sample_names[i]
        
        # Filter to standard chromosomes
        dmr_filtered <- dmr_data[dmr_data[["DMR_seqnames"]] %in% chr_data$chr, ]
        
        if(nrow(dmr_filtered) > 0) {
            # Create density bins for comparative plot
            density_data <- create_dmr_density_bins(dmr_filtered, chr_data, bin_size = 2e6)
            max_density <- max(density_data$density, na.rm = TRUE)
            
            if(max_density > 0) {
                circos.track(factors = chr_data$chr, ylim = c(0, max_density), 
                             panel.fun = function(x, y) {
                                 chr <- CELL_META$sector.index
                                 chr_density <- density_data[density_data$chr == chr, ]
                                 
                                 if(nrow(chr_density) > 0) {
                                     for(j in seq_len(nrow(chr_density))) {
                                         if(chr_density$density[j] > 0) {
                                             circos.rect(chr_density$start[j], 0, 
                                                        chr_density$end[j], chr_density$density[j], 
                                                        col = colors[current_sample], border = NA)
                                         }
                                     }
                                 }
                             }, track.height = 0.05, bg.border = "#b1b1b1")
            }
        }
    }
    
    # Add title
    title("Comparative Circos Plot: DMR Distribution", cex.main = 1.5)
    
    # Create legend for samples positioned at top-right corner
    if(length(sample_names) > 0) {
        lgd_samples <- Legend(at = sample_names, 
                             legend_gp = gpar(fill = colors[sample_names]),
                             title = "Comparisons",
                             title_gp = gpar(fontsize = 12, fontface = "bold"),
                             labels_gp = gpar(fontsize = 9),
                             grid_height = unit(0.8, "cm"),
                             grid_width = unit(0.8, "cm"))
        
        # Create a new page for legend
        #grid.newpage()
        pushViewport(viewport(width = 1, height = 1))
        draw(lgd_samples, x = unit(0.85, "npc"), y = unit(0.9, "npc"), just = c("left"))
        popViewport()
    }
    
    # Close PDF device
    dev.off()
    
    cat("Comparative Circos plot saved to:", pdf_file, "\n")
    
    # Clear circos parameters
    circos.clear()
    
    return(pdf_file)
}

# Main function to generate all circos plots
generate_dmr_circos_plots <- function(base_dir = "C:/PROJECTS/Shane/Harding_250611/dss", sample_names = c("IR2Gy6d_vs_NIR", "IR2Gy24h_vs_NIR"), dmr_folder = "DMRs") {
    
    dmr_data_list <- list()
    
    # Load data for each sample
    for(sample_name in sample_names) {
        working_dir <- file.path(base_dir, sample_name, dmr_folder)
        dmr_file <- file.path(working_dir, "DMRs_annotated_with_TE_dge.xlsx")
        output_dir <- file.path(base_dir, sample_name)
        if(file.exists(dmr_file)) {
            cat("Loading data for", sample_name, "\n")
            dmr_data <- read.xlsx(dmr_file)
            dmr_data_list[[sample_name]] <- dmr_data
            
            # Create individual circos plot
            create_dmr_circos(dmr_data, sample_name, output_dir)
            
        } else {
            cat("Warning: File not found:", dmr_file, "\n")
        }
    }
    
    # Create comparative plot if we have multiple samples
    if(length(dmr_data_list) > 1) {
        create_comparative_circos(dmr_data_list, names(dmr_data_list), base_dir)
    }
    
    cat("\nCircos plot generation completed!\n")
    cat("Individual plots saved in respective sample directories\n")
    if(length(dmr_data_list) > 1) {
        cat("Comparative plot saved in base directory:", base_dir, "\n")
    }
}

# Example usage:
# generate_dmr_circos_plots()

cat("Circos plotting functions loaded successfully!\n")
cat("Run generate_dmr_circos_plots() to create all plots\n")
cat("Or use create_dmr_circos() for individual sample plots\n")
