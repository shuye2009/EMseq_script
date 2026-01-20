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
library(ChIPseeker)
library(org.Hs.eg.db)

# Create TxDb from GTF
gtf_file <- "C:/PROJECTS/resource/T2T_CHM13/Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf"
txdb <- txdbmaker::makeTxDbFromGFF(gtf_file,
                        format = "gtf",
                        organism = "Homo sapiens")

# Read GTF to create gene_id to gene_name mapping
gtf_data <- rtracklayer::import(gtf_file)
gene_mapping <- as.data.frame(gtf_data) %>%
    filter(type == "gene") %>%
    dplyr::select(gene_id, gene_name) %>%
    distinct()
message("Created gene mapping with ", nrow(gene_mapping), " genes")

expression_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/RNAseq/deseq2/single_factor_analysis_gene_level"


dmr_method <- "dmrseq" #"dmrseq", "TE_targeted" "DSS"


emseq_dir <- file.path("C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq", dmr_method)

dmr_folder <- "DMR"
sample_dirs <- list.dirs(emseq_dir, full.names = FALSE, recursive = FALSE)

names(sample_dirs) <- sample_dirs

# redefine order of samples
sample_names <- c("IR2Gy24h_vs_NIR", "IR10Gy24h_vs_NIR", "IR10Gy24h_vs_IR2Gy24h", "IR2Gy6d_vs_NIR", "IR10Gy6d_vs_NIR", "IR10Gy6d_vs_IR2Gy6d", "IR2Gy6d_vs_IR2Gy24h", "IR10Gy6d_vs_IR10Gy24h")

sample_dirs <- sample_dirs[sample_names] 

gene_expression_dmr_plot_list <- list() # all regulation in DMR
gene_expression_in_dmr_plot_list <- list() # only up-regulated and down-regulated in DMR

regulated_count <- data.frame()

for(sample_name in sample_names){
    message("Processing ", sample_name, " ##############################################################################################")
    # remove leading "IR" from sample name
    SA <- gsub("^IR", "", sample_name)
    SA <- gsub("_IR", "_", SA)
    
    working_dir <- file.path(emseq_dir, sample_dirs[sample_name], dmr_folder)
    output_dir <- file.path(working_dir, "Integrate_gene_expression_with_DMR")
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
    
    # count up-regulated and down-regulated genes
    up <- sum(dge$regulation == "Up-regulated")
    down <- sum(dge$regulation == "Down-regulated")

    arow <- data.frame(Up_regulated = up, Down_regulated = down, comparison = sample_name)
    regulated_count <- rbind(regulated_count, arow)

    # import DMRs and TE annotations
    dmr <- read_tsv(file.path(working_dir, "DMR_table.tsv")) %>%
        mutate(direction = ifelse(stat > 0, "Hypermethylated", "Hypomethylated"))

    # annotate DMRs
    dmr_gr <- makeGRangesFromDataFrame(dmr, keep.extra.columns = TRUE)
    dmrAnno <- annotatePeak(dmr_gr, 
                         tssRegion=c(-1000, 0),
                         TxDb=txdb, 
                         annoDb="org.Hs.eg.db")

    # View annotation summary
    print(dmrAnno)
    annotated_dmr <- as.data.frame(dmrAnno) %>%
        mutate(annotation = gsub(" \\(.+\\)", "", annotation))
    
    # Map T2T gene IDs to gene symbols using the GTF mapping
    annotated_dmr <- annotated_dmr %>%
        left_join(gene_mapping, by = c("geneId" = "gene_id"))
    
    message("Mapped ", sum(!is.na(annotated_dmr$gene_name)), " out of ", nrow(annotated_dmr), " DMRs to gene symbols")
    
    gene_dmr_overlap <- full_join(dge, annotated_dmr, by = "gene_name")

    head(gene_dmr_overlap)

    #####################################################################################################
    ######### Count gene names by annotation for up and down regulated gene that overlap with DMRs #####
    #####################################################################################################
    # Debug: Check data before filtering
    message("Total rows in gene_dmr_overlap: ", nrow(gene_dmr_overlap))
    message("Rows with DMR overlap: ", sum(!is.na(gene_dmr_overlap$seqnames)))
    message("Rows with padj < 0.05: ", sum(gene_dmr_overlap$padj < 0.05, na.rm = TRUE))
    message("Rows with |log2FC| > 1: ", sum(abs(gene_dmr_overlap$log2FoldChange) > 1, na.rm = TRUE))
    
    # Filter for genes that overlap with DMRs (non-NA DMR columns indicate overlap)
    gene_dmr_overlap_filtered <- gene_dmr_overlap %>%
        filter(!is.na(seqnames)) 
    write.xlsx(gene_dmr_overlap_filtered, paste0(output_dir, "/gene_expression_in_DMR.xlsx"))
    message("Rows after filtering: ", nrow(gene_dmr_overlap_filtered))
    
    if(nrow(gene_dmr_overlap_filtered) == 0){
        message("No significant gene overlap with DMRs for ", sample_name)
        next
    }
    # Count distinct gene by regulation, and methylation status
    dmr_gene_counts <- gene_dmr_overlap_filtered %>%
        group_by(annotation, regulation, direction) %>%
        summarise(count = n_distinct(gene_name), .groups = "drop") %>%
        mutate(regulation = as.character(regulation),
               direction = as.character(direction))
    
    # Create complete combinations to ensure all facets show all families
    all_combinations <- expand.grid(
        annotation = unique(dmr_gene_counts$annotation),
        regulation = c("Up-regulated", "Down-regulated", "No-change"),
        direction = c("Hypermethylated", "Hypomethylated"),
        stringsAsFactors = FALSE
    )
    
    dmr_gene_counts_complete <- all_combinations %>%
        left_join(dmr_gene_counts, by = c("annotation", "regulation", "direction")) %>%
        replace_na(list(count = 0)) %>%
        arrange(desc(count))
    
    
    # Create faceted bar chart
    p_dmr <- ggplot(dmr_gene_counts_complete, 
                    aes(x = annotation, y = count, fill = regulation)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~direction, ncol = 2) +
        scale_fill_manual(values = c("Up-regulated" = "#E74C3C", "Down-regulated" = "#3498DB", "No-change" = "#95A5A6")) +
        labs(title = sample_name,
             x = "",
             y = "Number of genes Overlapping with DMRs",
             fill = "Regulation") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "top",
              strip.background = element_rect(fill = "gray90"),
              strip.text = element_text(face = "bold")) +
        ylim(0, 5000)
    
    # Save plot
    ggsave(paste0(output_dir, "/gene_expression_in_DMR_barplot_faceted.pdf"), 
           plot = p_dmr, width = 14, height = 6)
    ggsave(paste0(output_dir, "/gene_expression_in_DMR_barplot_faceted.png"), 
           plot = p_dmr, width = 14, height = 6, dpi = 300)
    
    # Save overlap table
    gene_dmr_overlap_filtered <- gene_dmr_overlap_filtered %>%
        dplyr::filter(regulation != "No-change") %>%
        dplyr::select(c("gene_name", "log2FoldChange", "padj", "regulation", "direction")) %>%
        dplyr::arrange(regulation, direction) %>%
        dplyr::distinct()
    write.xlsx(gene_dmr_overlap_filtered, paste0(output_dir, "/regulated_gene_expression_in_DMR.xlsx"))
    gene_expression_dmr_plot_list[[sample_name]] <- p_dmr
    
    message("Created DMR overlap bar chart and saved counts for ", sample_name)

    # only plot regulated TE expression in DMR
    regulated_dmr <- dmr_gene_counts_complete %>%
        filter(regulation != "No-change")
    
    # Create faceted bar chart
    p_regulated_dmr <- ggplot(regulated_dmr, 
                    aes(x = annotation, y = count, fill = regulation)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~direction, ncol = 2) +
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
              strip.text = element_text(face = "bold")) +
        ylim(0, 200)
    
    # Save plot
    ggsave(paste0(output_dir, "/gene_regulated_in_DMR_barplot_faceted.pdf"), 
           plot = p_regulated_dmr, width = 14, height = 6)
    ggsave(paste0(output_dir, "/gene_regulated_in_DMR_barplot_faceted.png"), 
           plot = p_regulated_dmr, width = 14, height = 6, dpi = 300)
    gene_expression_in_dmr_plot_list[[sample_name]] <- p_regulated_dmr
    message("Created regulated DMR overlap bar chart and saved counts for ", sample_name)
}

# Save all plots use cowplot as 2 x 4 multiplot
if(length(gene_expression_dmr_plot_list) > 0){
    multiplot_dmr <- plot_grid(plotlist = gene_expression_dmr_plot_list, ncol = 3, nrow = 3, labels = "AUTO")
    ggsave(paste0(emseq_dir, "/gene_expression_in_DMR_barplot_faceted_multiplot.pdf"), 
           plot = multiplot_dmr, width = 20, height = 15)
    ggsave(paste0(emseq_dir, "/gene_expression_in_DMR_barplot_faceted_multiplot.png"), 
           plot = multiplot_dmr, width = 20, height = 15, dpi = 300)
    message("Saved combined gene expression in DMR plot with ", length(gene_expression_dmr_plot_list), " panels")
}

if(length(gene_expression_in_dmr_plot_list) > 0){
    multiplot_dmr <- plot_grid(plotlist = gene_expression_in_dmr_plot_list, ncol = 3, nrow = 3, labels = "AUTO")
    ggsave(paste0(emseq_dir, "/gene_regulated_in_DMR_barplot_faceted_multiplot.pdf"), 
           plot = multiplot_dmr, width = 20, height = 15)
    ggsave(paste0(emseq_dir, "/gene_regulated_in_DMR_barplot_faceted_multiplot.png"), 
           plot = multiplot_dmr, width = 20, height = 15, dpi = 300)
    message("Saved combined regulated gene expression in DMR plot with ", length(gene_expression_in_dmr_plot_list), " panels")
}

# plot barchart of regulated_count
regulated_count_long <- regulated_count %>%
    pivot_longer(cols = -comparison, names_to = "regulation", values_to = "count") %>%
    mutate(regulation = as.character(regulation)) %>%
    mutate(comparison = factor(comparison, levels=sample_names))

p_regulated_count <- ggplot(regulated_count_long, 
                    aes(x = comparison, y = count, fill = regulation)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Up_regulated" = "#E74C3C", "Down_regulated" = "#3498DB")) +
    labs(title = "Differential Gene expression",
         x = "",
         y = "Number of genes",
         fill = "Regulation") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
          legend.position = "top",
          legend.text = element_text(size = 16),
          axis.title.y = element_text(size = 16)) 

# save plots

ggsave(paste0(emseq_dir, "/genes_regulated.pdf"), 
           plot = p_regulated_count, width = 20, height = 15)
ggsave(paste0(emseq_dir, "/genes_regulated.png"), 
           plot = p_regulated_count, width = 20, height = 15, dpi = 300)
    
# save regulated_count to csv
write.csv(regulated_count, paste0(emseq_dir, "/genes_regulated.csv"))
