rm(list = ls())

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(openxlsx)
library(ggplot2)
library(cowplot)
library(GenomicRanges)
library(plyranges)
library(GenomicFeatures)


expression_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/RNAseq/deseq2/single_factor_analysis_gene_level"


promoter_folder <- "DSS" #"dmrseq", "TE_targeted" "DSS"


emseq_dir <- file.path("C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq", promoter_folder)

target_folder <- "GSEA"
sample_dirs <- list.dirs(emseq_dir, full.names = FALSE, recursive = FALSE)

names(sample_dirs) <- sample_dirs

# redefine order of samples
sample_names <- c("IR2Gy24h_vs_NIR", "IR10Gy24h_vs_NIR", "IR10Gy24h_vs_IR2Gy24h", "IR2Gy6d_vs_NIR", "IR10Gy6d_vs_NIR", "IR10Gy6d_vs_IR2Gy6d", "IR2Gy6d_vs_IR2Gy24h", "IR10Gy6d_vs_IR10Gy24h")

sample_dirs <- sample_dirs[sample_names] 

gene_expression_promoter_plot_list <- list() # all regulation in promoter
gene_expression_in_regulated_promoter_plot_list <- list() # only up-regulated and down-regulated in promoter


for(sample_name in sample_names){
    message("Processing ", sample_name, " ##############################################################################################")
    # remove leading "IR" from sample name
    SA <- gsub("^IR", "", sample_name)
    SA <- gsub("_IR", "_", SA)
    
    working_dir <- file.path(emseq_dir, sample_dirs[sample_name], target_folder)
    output_dir <- file.path(working_dir, "Integrate_gene_expression_with_Promoter_GSEA")
    if(!dir.exists(output_dir)) dir.create(output_dir)
    # import differential transcript expression results
    dge <- read.csv(paste0(expression_dir, "/",SA, "/", SA, "_DESeq2_results_regular.csv")) %>%
        dplyr::filter(!is.na(padj)) %>%
        mutate(regulation = case_when(
            log2FoldChange > 1  & padj < 0.05 ~ "Up-regulated",
            log2FoldChange < -1 & padj < 0.05 ~ "Down-regulated",
            TRUE ~ "No-change"
        )) %>%
        dplyr::select(-direction)
    head(dge)
    

    # import DMRs and TE annotations
    promoter <- read_tsv(file.path(working_dir, paste0("GSEA_result_", sample_name, ".tsv"))) %>%
        dplyr::filter(!is.na(padj)) %>%
        mutate(direction = case_when(
            NES > 0 & padj < 0.01 ~ "Hypermethylated", 
            NES < 0 & padj < 0.01 ~ "Hypomethylated",
            TRUE ~ "Not_significant")) %>%
        mutate(gene_name = sapply(strsplit(pathway, "_"), function(x) x[1])) %>%
        mutate(gene_id = sapply(strsplit(pathway, "_"), function(x) x[2]))

    
    
    gene_promoter_overlap <- inner_join(dge, promoter, by = "gene_name")

    head(gene_promoter_overlap)

    # plot correlation between gene expression and promoter methylation
    p_cor <- ggplot(gene_promoter_overlap, aes(x = NES, y = log2FoldChange)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        labs(title = "Gene Expression vs Promoter Methylation",
             y = "Log2 Fold Change",
             x = "NES") +
        theme_minimal()
    
    ggsave(paste0(output_dir, "/gene_expression_vs_promoter_methylation.png"), p_cor, width = 10, height = 8)
    ggsave(paste0(output_dir, "/gene_expression_vs_promoter_methylation.pdf"), p_cor, width = 10, height = 8)
    
    #####################################################################################################
    ######### Count gene names by annotation for up and down regulated gene that overlap with DMRs #####
    #####################################################################################################

    
    
    write.xlsx(gene_promoter_overlap, paste0(output_dir, "/gene_expression_in_methylated_Promoter.xlsx"))
    message("Rows after filtering: ", nrow(gene_promoter_overlap))
    
    if(nrow(gene_promoter_overlap) == 0){
        message("No significant gene overlap with DMRs for ", sample_name)
        next
    }
    # Count distinct gene by regulation, and methylation status
    promoter_gene_counts <- gene_promoter_overlap %>%
        group_by(regulation, direction) %>%
        summarise(count = n_distinct(gene_name), .groups = "drop") %>%
        mutate(regulation = as.character(regulation),
               direction = as.character(direction))
    
    # Create complete combinations to ensure all facets show all families
    all_combinations <- expand.grid(
        regulation = c("Up-regulated", "Down-regulated", "No-change"),
        direction = c("Hypermethylated", "Hypomethylated", "Not_significant"),
        stringsAsFactors = FALSE
    )
    
    promoter_gene_counts_complete <- all_combinations %>%
        left_join(promoter_gene_counts, by = c("regulation", "direction")) %>%
        replace_na(list(count = 0)) %>%
        arrange(regulation, direction)
    
    
    # Create faceted bar chart
    p_promoter <- ggplot(promoter_gene_counts_complete, 
                    aes(x = direction, y = count, fill = regulation)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = c("Up-regulated" = "#E74C3C", "Down-regulated" = "#3498DB", "No-change" = "#95A5A6")) +
        labs(title = sample_name,
             x = "",
             y = "Number of genes regulated",
             fill = "Regulation") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "top",
              strip.background = element_rect(fill = "gray90"),
              strip.text = element_text(face = "bold")) 
    
    # Save plot
    ggsave(paste0(output_dir, "/gene_expression_in_methylated_Promoter_barplot.pdf"), 
           plot = p_promoter, width = 14, height = 6)
    ggsave(paste0(output_dir, "/gene_expression_in_methylated_Promoter_barplot.png"), 
           plot = p_promoter, width = 14, height = 6, dpi = 300)
    gene_expression_promoter_plot_list[[sample_name]] <- p_promoter
    message("Created DMR overlap bar chart and saved counts for ", sample_name)
    
    # Create contingency table: regulation (rows) vs methylation status (columns)
    cor_count <- data.frame(promoter_gene_counts_complete) %>%
        pivot_wider(names_from = direction, values_from = count)

    count_matrix <- as.matrix(cor_count %>%
        dplyr::select(-regulation))
    rownames(count_matrix) <- cor_count$regulation
    count_matrix <- count_matrix[c("Up-regulated", "Down-regulated"), c("Hypomethylated", "Hypermethylated")]

    message("\nContingency table for ", sample_name, ":")
    print(count_matrix)
    
    # Perform Fisher's exact test on the full contingency table
    # Tests if there's an association between regulation type and methylation status
    if(all(rowSums(count_matrix) > 0) && all(colSums(count_matrix) > 0)) {
        fisher_test_result <- fisher.test(count_matrix, simulate.p.value = TRUE)
        message("\nFisher's Exact Test Results:")
        message("P-value: ", format(fisher_test_result$p.value, scientific = TRUE, digits = 3))
        print(fisher_test_result)
        
        # Save the test results
        fisher_results_df <- data.frame(
            comparison = sample_name,
            p_value = fisher_test_result$p.value,
            method = fisher_test_result$method
        )
        write.xlsx(list(
            contingency_table = as.data.frame(count_matrix),
            fisher_test = fisher_results_df
        ), paste0(output_dir, "/promoter_methylation_expression_fisher_test.xlsx"), rowNames = TRUE)
    }

    # only plot regulated TE expression in DMR
    regulated_promoter <- promoter_gene_counts_complete %>%
        filter(regulation %in% c("Up-regulated", "Down-regulated"), direction %in% c("Hypomethylated", "Hypermethylated"))
    
    # Create faceted bar chart
    p_regulated_promoter <- ggplot(regulated_promoter, 
                    aes(x = direction, y = count, fill = regulation)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = c("Up-regulated" = "#E74C3C", "Down-regulated" = "#3498DB")) +
        labs(title = sample_name,
             x = "",
             y = "Number of regulated genes Overlapping with DMRs",
             fill = "Regulation") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "top",
              strip.background = element_rect(fill = "gray90"),
              strip.text = element_text(face = "bold")) 
    
    # Save plot
    ggsave(paste0(output_dir, "/gene_regulated_in_methylated_Promoter_barplot.pdf"), 
           plot = p_regulated_promoter, width = 14, height = 6)
    ggsave(paste0(output_dir, "/gene_regulated_in_methylated_Promoter_barplot.png"), 
           plot = p_regulated_promoter, width = 14, height = 6, dpi = 300)
    gene_expression_in_regulated_promoter_plot_list[[sample_name]] <- p_regulated_promoter
    message("Created regulated DMR overlap bar chart and saved counts for ", sample_name)
}

# Save all plots use cowplot as 2 x 4 multiplot
if(length(gene_expression_promoter_plot_list) > 0){
    multiplot_promoter <- plot_grid(plotlist = gene_expression_promoter_plot_list, ncol = 3, nrow = 3, labels = "AUTO")
    ggsave(paste0(emseq_dir, "/gene_expression_in_methylated_Promoter_barplot_multiplot.pdf"), 
           plot = multiplot_promoter, width = 20, height = 15)
    ggsave(paste0(emseq_dir, "/gene_expression_in_methylated_Promoter_barplot_multiplot.png"), 
           plot = multiplot_promoter, width = 20, height = 15, dpi = 300)
    message("Saved combined gene expression in methylated promoter plot with ", length(gene_expression_promoter_plot_list), " panels")
}

if(length(gene_expression_in_regulated_promoter_plot_list) > 0){
    multiplot_promoter <- plot_grid(plotlist = gene_expression_in_regulated_promoter_plot_list, ncol = 3, nrow = 3, labels = "AUTO")
    ggsave(paste0(emseq_dir, "/gene_expression_in_methylated_Promoter_barplot_regulated_multiplot.pdf"), 
           plot = multiplot_promoter, width = 20, height = 15)
    ggsave(paste0(emseq_dir, "/gene_expression_in_methylated_Promoter_barplot_regulated_multiplot.png"), 
           plot = multiplot_promoter, width = 20, height = 15, dpi = 300)
    message("Saved combined gene expression in methylated promoter plot with ", length(gene_expression_in_regulated_promoter_plot_list), " panels")
}


