rm(list=ls())

library(clusterProfiler)
library(bit64)
library(dplyr)
library(ggplot2)
# to plot and transform data
library(readxl)
library(purrr)
library(tibble)
library(tidyr)
library(readr)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(enrichplot)
library(GOSemSim)
library(GenomicRanges)
library(plyranges)
# Function to create volcano plot with labeled top genes
library(ggrepel)  # For better gene label positioning
library(EnhancedVolcano)


source("C:/PROJECTS/Repositories/Functional_enrichment/Compare_cluster_lib.R")
source("C:/PROJECTS/Repositories/Functional_enrichment/Function_analysis_for_RNAseq_lib.R")
source("C:/PROJECTS/Repositories/Functional_enrichment/DESeq2_RNAseq_RSEM_lib.R")
source("C:/PROJECTS/Repositories/Functional_enrichment/GO_KEGG_enrichment_lib.R")

genome <- "hg38"
TE <- TRUE
GSEA <- FALSE
TE_gtf <- "C:/PROJECTS/Shane/RNAseq/GRCh38_GENCODE_rmsk_TE.gtf/GRCh38_GENCODE_rmsk_TE.gtf"

TE_gr <- plyranges::read_gff(TE_gtf)
TE_df <- data.frame(gene_id = TE_gr$gene_id, 
                    gene_name = TE_gr$gene_name, 
                    transcript_id = TE_gr$transcript_id, 
                    family_id = TE_gr$family_id,
                    class_id = TE_gr$class_id)

########################################################################################################
############ Multi-factor DESeq2 analysis for time and dose effects #####################################
########################################################################################################

if(TRUE){
    # Multi-factor DESeq2 analysis for time and dose effects
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(dplyr)

    # Load count data
    wd <- "C:/PROJECTS/Shane/RNAseq/20250604_LH00244_0317_A22WFMGLT3_Harding_Shreya_RNAseq/results_TE/deseq2"
    setwd(wd)
    rsem_count <- read.delim("rsem_isoform_level_counts.tsv", 
                            header = TRUE, row.names = 1)
    if(TE){
        outdir <- file.path(wd, "multi_factor_analysis_TE_transcript") 
    }
    else{
        outdir <- file.path(wd, "multi_factor_analysis_transcript") 
    }
    if(!dir.exists(outdir)) dir.create(outdir)
    setwd(outdir)

    # Prepare count matrix (remove gene name column, keep only count data)
    count_matrix <- rsem_count %>%
     dplyr::select(starts_with("X"))
    count_matrix <- round(as.matrix(count_matrix))  # Round to integers for DESeq2

    # Create sample design
    sample_design <- data.frame(time=rep(rep(c("24h", "6d"), each=3), times=2),
                                dose=rep(c("2Gy", "10Gy"), each=6))
    rownames(sample_design) <- colnames(count_matrix)

    sample_design$time <- factor(sample_design$time, levels = c("24h", "6d"))
    sample_design$dose <- factor(sample_design$dose, levels = c("2Gy", "10Gy"))

    # Test for time dose and interaction effect (time:dose)
    cat("\n=== Testing Interaction Effect ===\n")
    dds_interaction <- DESeqDataSetFromMatrix(countData = count_matrix,
                                            colData = sample_design,
                                            design = ~ time + dose + time:dose)
    
    keep <- rowSums(counts(dds_interaction) >= 10) >= 3
    dds_interaction <- dds_interaction[keep, ]
    dds_interaction <- DESeq(dds_interaction)

    # Extract results for main effects
    # Time effect (comparing different time points)
    res_time <- results(dds_interaction, name = "time_6d_vs_24h")  # 6 days vs 24 hours

    # Dose effect (comparing different doses)
    res_dose <- results(dds_interaction, name = "dose_10Gy_vs_2Gy")   # 10 Gy vs 2 Gy

    # Summary of results
    cat("=== Time Effect (6d vs 24h) ===\n")
    summary(res_time, alpha = 0.05)

    cat("\n=== Dose Effect (10Gy vs 2Gy) ===\n")
    summary(res_dose, alpha = 0.05)

    # Create results data frames with gene names
    res_time_df <- data.frame(res_time)
    res_time_df$transcript_id <- rownames(res_time_df)
    res_time_df <- res_time_df[order(res_time_df$padj), ]
    if(TE){
        res_time_df_TE <- res_time_df %>%
            dplyr::filter(transcript_id %in% TE_df$transcript_id) %>%
            dplyr::left_join(TE_df[, c("transcript_id", "gene_id", "gene_name", "family_id", "class_id")], by = "transcript_id") %>%
            dplyr::distinct(transcript_id, .keep_all = TRUE)
        res_time_df <- res_time_df %>%
            dplyr::filter(!transcript_id %in% TE_df$transcript_id)
    }

    res_dose_df <- data.frame(res_dose)
    res_dose_df$transcript_id <- rownames(res_dose_df)
    res_dose_df <- res_dose_df[order(res_dose_df$padj), ]
    if(TE){
        res_dose_df_TE <- res_dose_df %>%
            dplyr::filter(transcript_id %in% TE_df$transcript_id) %>%
            dplyr::left_join(TE_df[, c("transcript_id", "gene_id", "gene_name", "family_id", "class_id")], by = "transcript_id") %>%
            dplyr::distinct(transcript_id, .keep_all = TRUE)
        res_dose_df <- res_dose_df %>%
            dplyr::filter(!transcript_id %in% TE_df$transcript_id)
    }

    # Save results
    write.csv(res_time_df, "time_6d_vs_24h_DESeq2_results.csv", row.names = FALSE)
    write.csv(res_dose_df, "dose_10Gy_vs_2Gy_DESeq2_results.csv", row.names = FALSE)
    if(TE){
        write.csv(res_time_df_TE, "time_6d_vs_24h_DESeq2_results_TE.csv", row.names = FALSE)
        write.csv(res_dose_df_TE, "dose_10Gy_vs_2Gy_DESeq2_results_TE.csv", row.names = FALSE)
    }

    # PCA plot
    vst <- vst(dds_interaction, blind = FALSE)
    pca_data <- plotPCA(vst, intgroup = c("time", "dose"), returnData = TRUE)
    percentVar <- round(100 * attr(pca_data, "percentVar"))

    ggplot(pca_data, aes(PC1, PC2, color = time, shape = dose)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = name), size = 3, max.overlaps = Inf) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA - Time and Dose Effects") +
    theme_minimal()

    ggsave("PCA_time_dose_effects.pdf", width = 8, height = 6, dpi = 300)

    # Heatmap of top variable genes
    vst_df <- assay(vst)
    if(TE){vst_df <- vst_df[row.names(vst_df) %in% TE_df$transcript_id, ]}
    top_var_genes <- head(order(rowVars(vst_df), decreasing = TRUE), 5000)
    heatmap_matrix <- vst_df[top_var_genes, ]
    #rownames(heatmap_matrix) <- rownames(vst_df)[top_var_genes]

    # Annotation for heatmap
    annotation_col <- data.frame(
    Time = sample_design$time,
    Dose = sample_design$dose,
    row.names = colnames(heatmap_matrix)
    )

    pheatmap(heatmap_matrix,
            annotation_col = annotation_col,
            scale = "row",
            clustering_distance_rows = "correlation",
            clustering_distance_cols = "correlation",
            filename = "heatmap_top5000_variable_genes.pdf",
            width = 10, height = 8)

     # Likelihood ratio test for interaction
    dds_lrt <- DESeq(dds_interaction, test = "LRT", reduced = ~ time + dose)
    resultsNames(dds_lrt)
    # Extract results for interaction effect
    res_interaction <- results(dds_lrt)

    cat("Genes with significant time:dose interaction (p < 0.05):\n")
    sig_interaction <- sum(res_interaction$padj < 0.05, na.rm = TRUE)
    cat(paste("Number of genes:", sig_interaction, "\n"))

    # Save interaction results
    res_interaction_df <- data.frame(res_interaction)
    res_interaction_df$transcript_id <- rownames(res_interaction_df)
    res_interaction_df <- res_interaction_df[order(res_interaction_df$padj), ]
    if(TE){
        res_interaction_df_TE <- res_interaction_df %>%
            dplyr::filter(transcript_id %in% TE_df$transcript_id) %>%
            dplyr::left_join(TE_df[, c("transcript_id", "gene_id", "gene_name", "family_id", "class_id")], by = "transcript_id") %>%
            dplyr::distinct(transcript_id, .keep_all = TRUE)
        res_interaction_df <- res_interaction_df %>%
            dplyr::filter(!transcript_id %in% TE_df$transcript_id)
    }
    write.csv(res_interaction_df, "time_dose_interaction_DESeq2_results.csv", row.names = FALSE)
    if(TE){
        write.csv(res_interaction_df_TE, "time_dose_interaction_DESeq2_results_TE.csv", row.names = FALSE)
    }
    
    # Generate volcano plots for the main effects and interaction
    cat("\n=== Generating Volcano Plots ===\n")
    volcano_time <- create_volcano_plot(res_time_df, "time_6d_vs_24h", 50, genome="hg38", label_col = "transcript_id", plot_gene_type = FALSE)
    volcano_dose <- create_volcano_plot(res_dose_df, "dose_10Gy_vs_2Gy", 50, genome="hg38", label_col = "transcript_id", plot_gene_type = FALSE)
    volcano_interaction <- create_volcano_plot(res_interaction_df, "time_dose_interaction", 50, genome="hg38", label_col = "transcript_id", plot_gene_type = FALSE)
    if(TE){
        volcano_time_TE <- create_volcano_plot(res_time_df_TE, "time_6d_vs_24h_TE", 50, genome="hg38", label_col = "transcript_id", plot_gene_type = FALSE)
        volcano_dose_TE <- create_volcano_plot(res_dose_df_TE, "dose_10Gy_vs_2Gy_TE", 50, genome="hg38", label_col = "transcript_id", plot_gene_type = FALSE)
        volcano_interaction_TE <- create_volcano_plot(res_interaction_df_TE, "time_dose_interaction_TE", 50, genome="hg38", label_col = "transcript_id", plot_gene_type = FALSE)
    }

    cat("Volcano plots saved!\n")

    # GSEA plot
    if(GSEA){
        create_gsea_plot(res_time_df, "Time_Effect_6d_vs_24h")
        create_gsea_plot(res_dose_df, "Dose_Effect_10Gy_vs_2Gy")
        create_gsea_plot(res_interaction_df, "Time_Dose_Interaction")
    }

    cat("\n=== Analysis Complete ===\n")
    cat("Results saved:\n")
    cat("- time_6d_vs_24h_DESeq2_results.csv\n")
    cat("- dose_10Gy_vs_2Gy_DESeq2_results.csv\n")
    cat("- time_dose_interaction_DESeq2_results.csv\n")
    cat("- volcano_time_6d_vs_24h.pdf\n")
    cat("- volcano_dose_10Gy_vs_2Gy.pdf\n")
    cat("- volcano_interaction_time_dose.pdf\n")
    cat("- PCA_time_dose_effects.pdf\n")
    cat("- heatmap_top5000_variable_genes.pdf\n")
}


########################################################################################################
############ Single-factor DESeq2 analysis for condition effects #####################################
########################################################################################################

if(TRUE){
    # Single-factor DESeq2 analysis for condition effects
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(dplyr)

    # Load count data
    wd <- "C:/PROJECTS/Shane/RNAseq/20250604_LH00244_0317_A22WFMGLT3_Harding_Shreya_RNAseq/results_TE/deseq2"
    setwd(wd)
    rsem_count <- read.delim("rsem_isoform_level_counts.tsv", 
                            header = TRUE, row.names = 1)
    if(TE){
        outdir <- file.path(wd, "single_factor_analysis_TE_transcript")
    }
    else{
        outdir <- file.path(wd, "single_factor_analysis_transcript")
    }
    if(!dir.exists(outdir)) dir.create(outdir)

    # Prepare count matrix (remove gene name column, keep only count data)
    count_matrix <- rsem_count[, -1]
    count_matrix <- round(as.matrix(count_matrix))  # Round to integers for DESeq2

    # Remove leading X from column names
    colnames(count_matrix) <- gsub("^X", "", colnames(count_matrix))

    # Create sample design
    sample_design <- read.delim("C:/PROJECTS/Shane/RNAseq/20250604_LH00244_0317_A22WFMGLT3_Harding_Shreya_RNAseq/config/design.tsv", header = TRUE)
    rownames(sample_design) <- sample_design$sample

    sample_design$condition <- factor(sample_design$condition, levels = c("NIR", "2Gy24h", "2Gy6d", "10Gy24h", "10Gy6d"))

    # Define contrasts for DESeq2 (will be used in results() function)
    contrast_list <- list(
      "2Gy24h_vs_NIR" = c("condition", "2Gy24h", "NIR"),
      "2Gy6d_vs_NIR" = c("condition", "2Gy6d", "NIR"),
      "10Gy24h_vs_NIR" = c("condition", "10Gy24h", "NIR"),
      "10Gy6d_vs_NIR" = c("condition", "10Gy6d", "NIR"),
      "2Gy6d_vs_2Gy24h" = c("condition", "2Gy6d", "2Gy24h"),
      "10Gy24h_vs_2Gy24h" = c("condition", "10Gy24h", "2Gy24h"),
      "10Gy6d_vs_2Gy6d" = c("condition", "10Gy6d", "2Gy6d"),
      "10Gy6d_vs_10Gy24h" = c("condition", "10Gy6d", "10Gy24h")
    )

    # Test for condition effect
    cat("\n=== Testing Condition Effect ===\n")
    dds_condition <- DESeqDataSetFromMatrix(countData = count_matrix,
                                            colData = sample_design,
                                            design = ~ condition)
    
    keep <- rowSums(counts(dds_condition) >= 10) >= 3
    dds_condition <- dds_condition[keep, ]
    dds_condition <- DESeq(dds_condition)

    resultsNames(dds_condition)

    

    # PCA plot
    vst <- vst(dds_condition, blind = FALSE)
    pca_data <- plotPCA(vst, intgroup = c("condition"), returnData = TRUE, ntop = 5000)
    percentVar <- round(100 * attr(pca_data, "percentVar"))

    library(ggrepel)  # For better text positioning
    
    ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = name), size = 3, max.overlaps = Inf) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA - Condition Effects") +
    theme_minimal()

    ggsave(file.path(outdir, "PCA_condition_effects.pdf"), width = 8, height = 6, dpi = 300)

    # Heatmap of top variable genes
    vst_df <- assay(vst)
    if(TE){
        vst_df <- vst_df[row.names(vst_df) %in% TE_df$transcript_id, ]
    }
    top_var_genes <- head(order(rowVars(vst_df), decreasing = TRUE), 5000)
    heatmap_matrix <- vst_df[top_var_genes, ]
    #rownames(heatmap_matrix) <- rownames(vst_df)[top_var_genes]

    # Annotation for heatmap
    annotation_col <- data.frame(
    condition = sample_design$condition,
    row.names = colnames(heatmap_matrix)
    )

    pheatmap(heatmap_matrix,
            annotation_col = annotation_col,
            scale = "row",
            clustering_distance_rows = "correlation",
            clustering_distance_cols = "correlation",
            filename = file.path(outdir, "heatmap_top5000_variable_genes.pdf"),
            width = 10, height = 8)

   # Extract results for all contrasts
    results_list <- list()
    for(contrast_name in names(contrast_list)) {
        cat(paste("\n=== Processing contrast:", contrast_name, "==="))
        contrast_dir <- file.path(outdir, contrast_name)
        if(!dir.exists(contrast_dir)) dir.create(contrast_dir)
        setwd(contrast_dir)
        # Extract results using contrast specification
        res <- results(dds_condition, contrast = contrast_list[[contrast_name]])
        
        # Summary of results
        cat(paste("\n=== Summary for", contrast_name, "==="))
        summary(res, alpha = 0.05)
        
        # Create results data frame with gene names
        res_df <- data.frame(res)
        res_df$transcript_id <- rownames(res_df)
        res_df <- res_df[order(res_df$padj), ]
        if(TE){
            res_df_TE <- res_df %>%
                dplyr::filter(transcript_id %in% TE_df$transcript_id) %>%
                dplyr::left_join(TE_df[, c("transcript_id", "gene_id", "gene_name", "family_id", "class_id")], by = "transcript_id") %>%
                dplyr::distinct(transcript_id, .keep_all = TRUE)
            res_df <- res_df %>%
                dplyr::filter(!transcript_id %in% TE_df$transcript_id)
        }
        
        # Store in results list
        results_list[[contrast_name]] <- res_df
        
        # Save results
        output_file <- paste0(contrast_name, "_DESeq2_results.csv")
        write.csv(res_df, output_file, row.names = FALSE)
        if(TE){
            write.csv(res_df_TE, paste0(contrast_name, "_DESeq2_results_TE.csv"), row.names = FALSE)
        }
        cat(paste("\nSaved:", basename(output_file)))

         # Generate volcano plots for the main effects and interaction
        cat("\n=== Generating Volcano Plots ===\n")
        
        # contrast effect volcano plot
        volcano_contrast <- create_volcano_plot(res_df, contrast_name, 50, genome="hg38", label_col = "transcript_id", plot_gene_type = FALSE)
        if(TE){
            volcano_contrast_TE <- create_volcano_plot(res_df_TE, paste0(contrast_name, "_TE"), 50, genome="hg38", label_col = "transcript_id", plot_gene_type = FALSE)
        }
        
        # GSEA plot
        if(GSEA){
            create_gsea_plot(res_df, contrast_name)
        }
    } 
   
}

