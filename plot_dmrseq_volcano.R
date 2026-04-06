
library(data.table)
library(ggplot2)
library(dplyr)
library(openxlsx)

# ==============================================================================
# Settings and Paths
# ==============================================================================

dmrseq_base_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/dmrseq"
output_dir <- file.path(dmrseq_base_dir, "dmrseq_Volcano_Plots")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# Process Comparisons
# ==============================================================================

# Get list of comparisons from dmrseq directory
comparisons <- list.dirs(dmrseq_base_dir, full.names = FALSE, recursive = FALSE)
# Filter to keep only likely comparison folders (avoiding files or standard folders like 'plot')
comparisons <- comparisons[grepl("_vs_", comparisons)]

message("Found comparisons: ", paste(comparisons, collapse = ", "))

for (comp in comparisons) {
    message("\nProcessing: ", comp)
    
    dmr_file <- file.path(dmrseq_base_dir, comp, "DMRs", "DMRs_annotated.xlsx")
    
    if (file.exists(dmr_file)) {
        # Read Data
        dmrs <- read.xlsx(dmr_file)
        
        # Check columns
        if (!all(c("betaCoefficient", "p.value", "q.value") %in% colnames(dmrs))) {
            warning("  Missing required columns (betaCoefficient, p.value, q.value) in ", comp)
            next
        }
        
        # Prepare Plot Data
        # Handle p=0 by replacing with minimum non-zero p-value / 10 or similar, or just cap
        min_nonzero_p <- min(dmrs$p.value[dmrs$p.value > 0], na.rm = TRUE)
        plot_data <- dmrs %>%
            mutate(
                plot_pval = ifelse(p.value == 0, min_nonzero_p * 0.1, p.value),
                logP = -log10(plot_pval),
                Direction = case_when(
                    p.value < 0.05 & betaCoefficient > 0.5 ~ "Hypermethylated",
                    p.value < 0.05 & betaCoefficient < -0.5 ~ "Hypomethylated",
                    TRUE ~ "Not Significant"
                )
            )
        
        # Count significant
        n_hyper <- sum(plot_data$Direction == "Hypermethylated")
        n_hypo  <- sum(plot_data$Direction == "Hypomethylated")
        message("  Hyper: ", n_hyper, " | Hypo: ", n_hypo)
        
        # Plot
        p <- ggplot(plot_data, aes(x = betaCoefficient, y = logP, color = Direction)) +
            geom_point(alpha = 0.6, size = 1.5) +
            scale_color_manual(values = c(
                "Hypermethylated" = "#E74C3C", # Red
                "Hypomethylated"  = "#3498DB", # Blue
                "Not Significant" = "grey70"
            )) +
            labs(title = paste("DMR Volcano Plot:", comp),
                 subtitle = paste("Hyper:", n_hyper, "| Hypo:", n_hypo, "(p < 0.05, |beta| > 0.5)"),
                 x = "BetaCoefficient (Effect Size)",
                 y = "-log10(p-value)") +
            theme_bw() +
            theme(
                plot.title = element_text(hjust = 0.5, face = "bold"),
                plot.subtitle = element_text(hjust = 0.5),
                legend.position = "top"
            ) +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.3) + # Raw p=0.05 line visual guide
            geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", alpha = 0.3)
            
        # Save
        out_prefix <- file.path(output_dir, paste0("Volcano_beta_pval_", comp))
        ggsave(paste0(out_prefix, ".pdf"), p, width = 8, height = 6)
        ggsave(paste0(out_prefix, ".png"), p, width = 8, height = 6, dpi = 300)
        
        message("  Saved plots to ", output_dir)
        
    } else {
        warning("  DMR table not found: ", dmr_file)
    }
}

message("\nDone.")
