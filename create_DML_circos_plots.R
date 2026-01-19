# DML Circos plots using circlize
library(circlize)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(grid)

# Chromosome data
prepare_chromosome_data <- function() {
    data.frame(
        chr = paste0("chr", c(1:22, "X", "Y")),
        length = c(248956422, 242193529, 198295559, 190214555, 181538259,
                   170805979, 159345973, 145138636, 138394717, 133797422,
                   135086622, 133275309, 114364328, 107043718, 101991189,
                   90338345, 83257441, 80373285, 58617616, 64444167,
                   46709983, 50818468, 156040895, 57227415)
    )
}

# DML density bins
create_dml_density_bins <- function(dml_data, chr_data, bin_size = 1e6) {
    density_data <- data.frame()
    for(chr in chr_data$chr) {
        chr_dmls <- dml_data[dml_data$chr == chr, ]
        if(nrow(chr_dmls) > 0) {
            chr_length <- chr_data$length[chr_data$chr == chr]
            bins <- seq(0, chr_length, by = bin_size)
            bin_counts <- rep(0, length(bins) - 1)
            for(i in seq_len(length(bins) - 1)) {
                bin_counts[i] <- sum(chr_dmls$pos >= bins[i] & chr_dmls$pos < bins[i + 1])
            }
            chr_density <- data.frame(
                chr = chr, start = bins[-length(bins)], end = bins[-1],
                count = bin_counts, density = bin_counts / (bin_size / 1e6)
            )
            density_data <- rbind(density_data, chr_density)
        }
    }
    return(density_data)
}

# Main DML Circos function
create_dml_circos <- function(dml_data, sample_name, output_dir) {
    cat("Creating DML Circos plot for", sample_name, "\n")
    
    chr_data <- prepare_chromosome_data()
    
    # Add chr prefix if needed
    if(!any(grepl("^chr", dml_data$chr))) {
        dml_data$chr <- paste0("chr", dml_data$chr)
    }
    
    # Filter to standard chromosomes
    dml_filtered <- dml_data[dml_data$chr %in% chr_data$chr, ]
    if(nrow(dml_filtered) == 0) return(NULL)
    
    # Filter significant DMLs
    if("fdr" %in% colnames(dml_filtered)) {
        significant_dmls <- dml_filtered[dml_filtered$fdr < 0.05, ]
    } else {
        significant_dmls <- dml_filtered
    }
    
    # Create PDF
    pdf_file <- file.path(output_dir, paste0("circos_DML_", sample_name, ".pdf"))
    pdf(pdf_file, width = 12, height = 12)
    
    # Initialize circos
    circos.clear()
    circos.par("start.degree" = 90, "gap.degree" = 2)
    circos.initialize(factors = chr_data$chr, xlim = cbind(rep(0, nrow(chr_data)), chr_data$length))
    
    # Track 1: Chromosomes
    circos.track(factors = chr_data$chr, ylim = c(0, 1), 
                 panel.fun = function(x, y) {
                     circos.text(CELL_META$xcenter, CELL_META$ycenter, 
                                CELL_META$sector.index, cex = 0.8, facing = "clockwise", niceFacing = TRUE)
                 }, track.height = 0.05, bg.col = "grey90")
    
    # Track 2: DML density
    density_data <- create_dml_density_bins(dml_filtered, chr_data)
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
    }
    
    # Track 3: Methylation differences
    if("diff" %in% colnames(significant_dmls)) {
        # Sample for performance
        if(nrow(significant_dmls) > 3000) {
            set.seed(42)
            plot_dmls <- significant_dmls[sample(nrow(significant_dmls), 3000), ]
        } else {
            plot_dmls <- significant_dmls
        }
        
        diff_range <- range(plot_dmls$diff, na.rm = TRUE)
        col_fun <- colorRamp2(c(diff_range[1], 0, diff_range[2]), c("blue", "white", "red"))
        
        circos.track(factors = chr_data$chr, ylim = c(-1, 1), 
                     panel.fun = function(x, y) {
                         chr <- CELL_META$sector.index
                         chr_dmls <- plot_dmls[plot_dmls$chr == chr, ]
                         if(nrow(chr_dmls) > 0) {
                             for(i in seq_len(nrow(chr_dmls))) {
                                 pos <- chr_dmls$pos[i]
                                 diff_val <- chr_dmls$diff[i]
                                 if(!is.na(diff_val)) {
                                     y_pos <- ifelse(diff_val > 0, 0.7, -0.7)
                                     circos.points(pos, y_pos, col = col_fun(diff_val), pch = 16, cex = 0.2)
                                 }
                             }
                         }
                     }, track.height = 0.1, bg.border = "#b1b1b1")
    }
    
    title(paste("DML Circos Plot:", sample_name))
    
    # Add legends for tracks
    legend_list <- list()
    
    # Legend for DML density track
    if(max_density > 0) {
        lgd_density <- Legend(at = c("DML Density"), 
                             legend_gp = gpar(fill = "darkblue"),
                             title = "Track 1: DML Density",
                             title_gp = gpar(fontsize = 10, fontface = "bold"),
                             labels_gp = gpar(fontsize = 9))
        legend_list <- append(legend_list, list(lgd_density))
    }
    
    # Legend for methylation differences
    if("diff" %in% colnames(significant_dmls) && exists("col_fun")) {
        lgd_diff <- Legend(at = c(round(diff_range[1], 3), 0, round(diff_range[2], 3)),
                          col_fun = col_fun, 
                          title = "Track 2: Methylation Difference",
                          title_gp = gpar(fontsize = 10, fontface = "bold"),
                          labels_gp = gpar(fontsize = 9),
                          legend_height = unit(3, "cm"))
        legend_list <- append(legend_list, list(lgd_diff))
    }
    
    # Draw legends if any exist
    if(length(legend_list) > 0) {
        combined_legend <- packLegend(list = legend_list, direction = "vertical", gap = unit(1, "cm"))
        draw(combined_legend, x = unit(0.5, "npc"), y = unit(0.5, "npc"))
    }
    
    dev.off()
    circos.clear()
    
    cat("DML Circos plot saved to:", pdf_file, "\n")
    return(pdf_file)
}

# Comparative DML Circos plot
create_comparative_dml_circos <- function(base_dir, sample_names, output_dir = NULL) {
    cat("Creating comparative DML Circos plot\n")
    
    chr_data <- prepare_chromosome_data()
    if(is.null(output_dir)) output_dir <- base_dir
    
    # Load all DML data
    dml_list <- list()
    for(sample_name in sample_names) {
        dml_file <- file.path(base_dir, sample_name, paste0("DML_", sample_name, ".tsv"))
        if(file.exists(dml_file)) {
            dml_data <- read.table(dml_file, header = TRUE, sep = "\t")
            if(!any(grepl("^chr", dml_data$chr))) dml_data$chr <- paste0("chr", dml_data$chr)
            dml_filtered <- dml_data[dml_data$chr %in% chr_data$chr, ]
            if("fdr" %in% colnames(dml_filtered)) dml_filtered <- dml_filtered[dml_filtered$fdr < 0.05, ]
            dml_list[[sample_name]] <- dml_filtered
        }
    }
    
    if(length(dml_list) == 0) return(NULL)
    
    # Create plot
    pdf_file <- file.path(output_dir, "circos_DML_comparative.pdf")
    pdf(pdf_file, width = 12, height = 12)
    
    circos.clear()
    circos.par("start.degree" = 90, "gap.degree" = 2)
    circos.initialize(factors = chr_data$chr, xlim = cbind(rep(0, nrow(chr_data)), chr_data$length))
    
    # Chromosome track
    circos.track(factors = chr_data$chr, ylim = c(0, 1), 
                 panel.fun = function(x, y) {
                     circos.text(CELL_META$xcenter, CELL_META$ycenter, 
                                CELL_META$sector.index, cex = 0.8, facing = "clockwise", niceFacing = TRUE)
                 }, track.height = 0.05, bg.col = "grey90")
    
    # Sample tracks
    colors <- brewer.pal(min(length(sample_names), 8), "Set2")
    if(length(sample_names) > 8) {
        colors <- rainbow(length(sample_names))
    }
    names(colors) <- sample_names
    
    for(i in seq_along(dml_list)) {
        sample_name <- names(dml_list)[i]
        dml_data <- dml_list[[i]]
        
        if(nrow(dml_data) > 0) {
            density_data <- create_dml_density_bins(dml_data, chr_data, bin_size = 2e6)
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
                                                        col = colors[sample_name], border = NA)
                                         }
                                     }
                                 }
                             }, track.height = 0.05, bg.border = "#b1b1b1")
            }
        }
    }
    
    title("Comparative Circos Plot: DML Distribution", cex.main = 1.5)
    
    # Create legend for samples positioned at top-right corner
    if(length(sample_names) > 0) {
        lgd_samples <- Legend(at = sample_names, 
                             legend_gp = gpar(fill = colors[sample_names]),
                             title = "Comparisons",
                             title_gp = gpar(fontsize = 12, fontface = "bold"),
                             labels_gp = gpar(fontsize = 9),
                             grid_height = unit(0.8, "cm"),
                             grid_width = unit(0.8, "cm"))
        
        # Position legend at top-right corner
        draw(lgd_samples, x = unit(0.85, "npc"), y = unit(0.9, "npc"), just = c("left"))
    }
    
    dev.off()
    circos.clear()
    
    cat("Comparative DML plot saved to:", pdf_file, "\n")
    return(pdf_file)
}

# Main function to generate all DML circos plots
generate_dml_circos_plots <- function(base_dir = "C:/PROJECTS/Shane/Harding_250611/dss", 
                                      sample_names = c("IR2Gy6d_vs_NIR", "IR2Gy24h_vs_NIR", "IR10Gy6d_vs_NIR", "IR10Gy24h_vs_NIR")) {
    
    cat("Generating DML Circos plots for", length(sample_names), "samples\n")
    
    # Create individual DML circos plots for each sample
    for(sample_name in sample_names) {
        working_dir <- file.path(base_dir, sample_name)
        dml_file <- file.path(working_dir, paste0("DML_", sample_name, ".tsv"))
        
        if(file.exists(dml_file)) {
            cat("Loading DML data for", sample_name, "\n")
            dml_data <- read.table(dml_file, header = TRUE, sep = "\t") %>%
                dplyr::filter(fdr < 0.05, abs(diff) > 0.1)
            
            # Create individual circos plot
            create_dml_circos(dml_data, sample_name, working_dir)
            
        } else {
            cat("Warning: DML file not found:", dml_file, "\n")
        }
    }
    
    # Create comparative plot
    create_comparative_dml_circos(base_dir, sample_names)
    
    cat("\nDML Circos plot generation completed!\n")
    cat("Individual plots saved in respective sample directories\n")
    cat("Comparative plot saved in base directory:", base_dir, "\n")
}

# Example usage:
# generate_dml_circos_plots()

cat("DML Circos functions loaded!\n")
