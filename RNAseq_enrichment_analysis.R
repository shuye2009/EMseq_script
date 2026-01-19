rm(list=ls())

library(clusterProfiler)
library(bit64)
library(dplyr)
library(ggplot2)
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
library(ggrepel)  # For better gene label positioning
library(EnhancedVolcano)
library(DESeq2)
library(pheatmap)



source("C:/PROJECTS/Repositories/Functional_enrichment/Compare_cluster_lib.R")
source("C:/PROJECTS/Repositories/Functional_enrichment/Function_analysis_for_RNAseq_lib.R")
source("C:/PROJECTS/Repositories/Functional_enrichment/DESeq2_RNAseq_RSEM_lib.R")
source("C:/PROJECTS/Repositories/Functional_enrichment/GO_KEGG_enrichment_lib.R")

genome <- "hg38"
TE <- TRUE # set to TRUE if the data have TE gene expression, FALSE for basic annotations
GSEA <- TRUE # set to TRUE if GSEA is to be performed

TE_gtf <- "C:/PROJECTS/resource/GRCh38_GENCODE_rmsk_TE.gtf"

TE_gr <- plyranges::read_gff(TE_gtf)
TE_df <- data.frame(gene_id = TE_gr$gene_id, 
                    gene_name = TE_gr$gene_name, 
                    transcript_id = TE_gr$transcript_id, 
                    family_id = TE_gr$family_id,
                    class_id = TE_gr$class_id)

# Convert TE_gr to BED format
TE_bed <- data.frame(
    chrom = as.character(seqnames(TE_gr)),
    chromStart = start(TE_gr) - 1,  # BED format is 0-based
    chromEnd = end(TE_gr),          # BED format end is 1-based (exclusive)
    name = TE_gr$transcript_id,
    score = width(TE_gr),           # Length as score
    strand = as.character(strand(TE_gr))
)

# Write BED file
write.table(TE_bed, "C:/PROJECTS/resource/GRCh38_GENCODE_rmsk_TE.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("TE annotations saved to GRCh38_GENCODE_rmsk_TE.bed\n")

# Function to perform family enrichment analysis
perform_family_enrichment <- function(res_df, analysis_name, padj_threshold = 0.05, log2fc_threshold = 1) {
    cat(paste("\n=== Performing Family Enrichment Analysis for", analysis_name, "===\n"))
    
    # Check if we have the required columns and data
    if(!("family_id" %in% colnames(res_df)) || nrow(res_df) == 0) {
        cat("Skipping family enrichment: no family_id column or empty dataframe\n")
        return(NULL)
    }
    
    # Classify genes as up-regulated, down-regulated, or not significant
    res_df$regulation <- "Not_significant"
    res_df$regulation[res_df$padj < padj_threshold & res_df$log2FoldChange > log2fc_threshold] <- "Up_regulated"
    res_df$regulation[res_df$padj < padj_threshold & res_df$log2FoldChange < -log2fc_threshold] <- "Down_regulated"
    
    # Get all TE families in the dataset
    all_families <- unique(res_df$family_id[!is.na(res_df$family_id)])
    
    if(length(all_families) == 0) {
        cat("No families found in the dataset\n")
        return(NULL)
    }
    
    # Initialize results dataframe
    family_enrichment <- data.frame(
        family_id = character(),
        regulation = character(),
        family_count = integer(),
        total_family_genes = integer(),
        background_count = integer(),
        total_background_genes = integer(),
        expected = numeric(),
        fold_enrichment = numeric(),
        p_value = numeric(),
        padj = numeric(),
        stringsAsFactors = FALSE
    )
    
    # Perform enrichment analysis for each family and regulation direction
    for(family in all_families){
        for(reg_type in c("Up_regulated", "Down_regulated")){
            # Count genes in this family with this regulation
            family_reg_count <- sum(res_df$family_id == family & res_df$regulation == reg_type, na.rm = TRUE)
            
            # Total genes in this family
            total_family_genes <- sum(res_df$family_id == family, na.rm = TRUE)
            
            # Count all genes with this regulation (background)
            background_reg_count <- sum(res_df$regulation == reg_type, na.rm = TRUE)
            
            # Total genes in background
            total_background_genes <- nrow(res_df)
            
            # Expected count under null hypothesis
            expected <- (total_family_genes * background_reg_count) / total_background_genes
            
            # Fold enrichment
            fold_enrichment <- ifelse(expected > 0, family_reg_count / expected, NA)
            
            # Fisher's exact test
            contingency_matrix <- matrix(c(
                family_reg_count,
                background_reg_count - family_reg_count,
                total_family_genes - family_reg_count,
                total_background_genes - background_reg_count - total_family_genes + family_reg_count
            ), nrow = 2)
            
            # Only perform test if we have sufficient counts
            if(all(contingency_matrix >= 0) && sum(contingency_matrix) > 0){
                fisher_test <- fisher.test(contingency_matrix, alternative = "greater")
                p_value <- fisher_test$p.value
            } else {
                p_value <- 1
            }
            
            # Add to results
            family_enrichment <- rbind(family_enrichment, data.frame(
                family_id = family,
                regulation = reg_type,
                family_count = family_reg_count,
                total_family_genes = total_family_genes,
                background_count = background_reg_count,
                total_background_genes = total_background_genes,
                expected = expected,
                fold_enrichment = fold_enrichment,
                p_value = p_value,
                padj = NA,
                stringsAsFactors = FALSE
            ))
        }
    }
    
    # Adjust p-values for multiple testing
    if(nrow(family_enrichment) > 0){
        family_enrichment$padj <- p.adjust(family_enrichment$p_value, method = "BH")
        family_enrichment <- family_enrichment[order(family_enrichment$padj, family_enrichment$p_value), ]
        
        # Save family enrichment results
        family_output_file <- paste0(analysis_name, "_family_enrichment_analysis.csv")
        write.csv(family_enrichment, family_output_file, row.names = FALSE)
        cat(paste("Saved family enrichment analysis:", basename(family_output_file), "\n"))
        
        # Print summary of significant enrichments
        sig_enrichments <- family_enrichment[family_enrichment$padj < 0.05 & family_enrichment$fold_enrichment > 1, ]
        if(nrow(sig_enrichments) > 0){
            cat(paste("Found", nrow(sig_enrichments), "significantly enriched family-regulation combinations (padj < 0.05):\n"))
            for(i in 1:min(10, nrow(sig_enrichments))){
                cat(paste(" -", sig_enrichments$family_id[i], "(", sig_enrichments$regulation[i], 
                        "): Fold enrichment =", round(sig_enrichments$fold_enrichment[i], 2), 
                        ", padj =", format(sig_enrichments$padj[i], scientific = TRUE, digits = 3), "\n"))
            }
        } else {
            cat("No significantly enriched family-regulation combinations found (padj < 0.05)\n")
        }
        
        return(family_enrichment)
    } else {
        cat("No enrichment results generated\n")
        return(NULL)
    }
}

# read in data


########################################################################################################
############ Multi-factor DESeq2 analysis for time and dose effects #####################################
########################################################################################################

if(FALSE){
    wd <- "C:/PROJECTS/Shane/RNAseq/20250604_LH00244_0317_A22WFMGLT3_Harding_Shreya_RNAseq/results_TE/deseq2"
    setwd(wd)
    rsem_count <- read.delim("rsem_gene_level_counts.tsv", 
                            header = TRUE, row.names = 1)
                            # Average counts for duplicated GENENAME using group_by and summarize
    rsem_count <- rsem_count %>%
    group_by(GENENAME) %>%
    summarize(across(starts_with("X"), mean), .groups = "drop") %>%
    as.data.frame()

    # Prepare count matrix (remove gene name column, keep only count data)
    rownames(rsem_count) <- rsem_count$GENENAME
    count_matrix <- rsem_count[, -1]  # Remove GENENAME column
    count_matrix <- round(as.matrix(count_matrix))  # Round to integers for DESeq2

    # Remove leading X from column names
    colnames(count_matrix) <- gsub("^X", "", colnames(count_matrix))
    
    if(TE){
        outdir <- file.path(wd, "multi_factor_analysis_TE_gene_level")
    }else{
        outdir <- file.path(wd, "multi_factor_analysis_gene_level")
    }
    if(!dir.exists(outdir)) dir.create(outdir)
    setwd(outdir)

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
    res_time_df$gene_name <- rownames(res_time_df)
    res_time_df <- res_time_df[order(res_time_df$padj), ]
    if(TE){
        res_time_df_TE <- res_time_df %>%
          dplyr::filter(gene_name %in% TE_df$gene_name) %>%
          dplyr::left_join(TE_df[, c("gene_name", "gene_id", "family_id", "class_id")], by = "gene_name") %>%
          dplyr::distinct(gene_name, .keep_all = TRUE)

        res_time_df <- res_time_df %>%
          dplyr::filter(!gene_name %in% TE_df$gene_name)
    }

    res_dose_df <- data.frame(res_dose)
    res_dose_df$gene_name <- rownames(res_dose_df)
    res_dose_df <- res_dose_df[order(res_dose_df$padj), ]
    if(TE){
        res_dose_df_TE <- res_dose_df %>%
        dplyr::filter(gene_name %in% TE_df$gene_name) %>%
        dplyr::left_join(TE_df[, c("gene_name", "gene_id", "family_id", "class_id")], by = "gene_name") %>%
        dplyr::distinct(gene_name, .keep_all = TRUE)

        res_dose_df <- res_dose_df %>%
          dplyr::filter(!gene_name %in% TE_df$gene_name)
    }

    # Save results
    write.csv(res_time_df, "time_6d_vs_24h_DESeq2_results.csv", row.names = FALSE)
    write.csv(res_dose_df, "dose_10Gy_vs_2Gy_DESeq2_results.csv", row.names = FALSE)
    if(TE) write.csv(res_time_df_TE, "time_6d_vs_24h_DESeq2_results_TE.csv", row.names = FALSE)
    if(TE) write.csv(res_dose_df_TE, "dose_10Gy_vs_2Gy_DESeq2_results_TE.csv", row.names = FALSE)

    if(TE && nrow(res_time_df_TE) > 0){
        perform_family_enrichment(res_time_df_TE, "time_6d_vs_24h_TE")
    }
    if(TE && nrow(res_dose_df_TE) > 0){
        perform_family_enrichment(res_dose_df_TE, "dose_10Gy_vs_2Gy_TE")
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
    if(TE){
        vst_df <- vst_df[row.names(vst_df) %in% TE_df$gene_name, ]
    }
    top_var_genes <- head(order(rowVars(vst_df), decreasing = TRUE), min(5000, nrow(vst_df)))
    heatmap_matrix <- vst_df[top_var_genes, ]

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
    res_interaction_df$gene_name <- rownames(res_interaction_df)
    res_interaction_df <- res_interaction_df[order(res_interaction_df$padj), ]
    if(TE){
        res_interaction_df_TE <- res_interaction_df %>%
            dplyr::filter(gene_name %in% TE_df$gene_name) %>%
            dplyr::left_join(TE_df[, c("gene_name", "gene_id", "family_id", "class_id")], by = "gene_name") %>%
            dplyr::distinct(gene_name, .keep_all = TRUE)

        res_interaction_df <- res_interaction_df %>%
          dplyr::filter(!gene_name %in% TE_df$gene_name)
    }
    write.csv(res_interaction_df, "time_dose_interaction_DESeq2_results.csv", row.names = FALSE)
    if(TE) write.csv(res_interaction_df_TE, "time_dose_interaction_DESeq2_results_TE.csv", row.names = FALSE)

    if(TE && nrow(res_interaction_df_TE) > 0){
        perform_family_enrichment(res_interaction_df_TE, "time_dose_interaction_TE")
    }
    
    # Generate volcano plots for the main effects and interaction
    cat("\n=== Generating Volcano Plots ===\n")
    volcano_time <- create_volcano_plot(res_time_df, "time_6d_vs_24h", 50, genome="hg38", plot_gene_type = TRUE)
    volcano_dose <- create_volcano_plot(res_dose_df, "dose_10Gy_vs_2Gy", 50, genome="hg38", plot_gene_type = TRUE)
    volcano_interaction <- create_volcano_plot(res_interaction_df, "time_dose_interaction", 50, genome="hg38", plot_gene_type = TRUE)

    if(TE){
        volcano_time_TE <- create_volcano_plot(res_time_df_TE, "time_6d_vs_24h_TE", 50, genome="hg38", plot_gene_type = FALSE)
        volcano_dose_TE <- create_volcano_plot(res_dose_df_TE, "dose_10Gy_vs_2Gy_TE", 50, genome="hg38", plot_gene_type = FALSE)
        volcano_interaction_TE <- create_volcano_plot(res_interaction_df_TE, "time_dose_interaction_TE", 50, genome="hg38", plot_gene_type = FALSE) 
    }

    cat("Volcano plots saved!")

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

if(FALSE){
    # Load count data
    wd <- "C:/PROJECTS/Shane/RNAseq/20250604_LH00244_0317_A22WFMGLT3_Harding_Shreya_RNAseq/results_TE/deseq2"
    setwd(wd)
    rsem_count <- read.delim("rsem_gene_level_counts.tsv", 
                            header = TRUE, row.names = 1)
                            # Average counts for duplicated GENENAME using group_by and summarize
    rsem_count <- rsem_count %>%
    group_by(GENENAME) %>%
    summarize(across(everything(), mean), .groups = "drop") %>%
    as.data.frame()

    # Prepare count matrix (remove gene name column, keep only count data)
    rownames(rsem_count) <- rsem_count$GENENAME
    count_matrix <- rsem_count[, -1]  # Remove GENENAME column
    count_matrix <- round(as.matrix(count_matrix))  # Round to integers for DESeq2

    # Remove leading X from column names
    colnames(count_matrix) <- gsub("^X", "", colnames(count_matrix))
    if(TE){
        outdir <- file.path(wd, "single_factor_analysis_TE_gene_level")
    }
    else{
        outdir <- file.path(wd, "single_factor_analysis_gene_level")
    }
    if(!dir.exists(outdir)) dir.create(outdir)

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
        vst_df <- vst_df[row.names(vst_df) %in% TE_df$gene_name, ]
    }
    top_var_genes <- head(order(rowVars(vst_df), decreasing = TRUE), min(5000, nrow(vst_df)))
    heatmap_matrix <- vst_df[top_var_genes, ]

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
        res_df$gene_name <- rownames(res_df)
        res_df <- res_df[order(res_df$padj), ]
        if(TE){
            res_df_TE <- res_df %>%
            dplyr::filter(gene_name %in% TE_df$gene_name) %>%
            dplyr::left_join(TE_df[, c("gene_name", "gene_id", "family_id", "class_id")], by = "gene_name") %>%
            dplyr::distinct(gene_name, .keep_all = TRUE)

            res_df <- res_df %>%
            dplyr::filter(!gene_name %in% TE_df$gene_name)
        }
        
        # Store in results list
        results_list[[contrast_name]] <- res_df
        
        # Save results
        output_file <- paste0(contrast_name, "_DESeq2_results.csv")
        write.csv(res_df, output_file, row.names = FALSE)
        if(TE) write.csv(res_df_TE, paste0(contrast_name, "_DESeq2_results_TE.csv"), row.names = FALSE)
        cat(paste("\nSaved:", basename(output_file)))

        if(TE && nrow(res_df_TE) > 0){
            perform_family_enrichment(res_df_TE, contrast_name)
        }

         # Generate volcano plots for the main effects and interaction
        cat("\n=== Generating Volcano Plots ===\n")
        
        # contrast effect volcano plot
        volcano_contrast <- create_volcano_plot(res_df, contrast_name, 50, genome="hg38", plot_gene_type = TRUE)
        if(TE) volcano_contrast_TE <- create_volcano_plot(res_df_TE, paste0(contrast_name, "_TE"), 50, genome="hg38", plot_gene_type = FALSE)
        
        # GSEA plot
        if(GSEA){
            create_gsea_plot(res_df, contrast_name)
        }
    } 
   
}

if(FALSE){
  id2name_file <- "C:/PROJECTS/Ben/RNAsplicing/results/geneid2name.tsv"
  
  hd <- file.path("C:/PROJECTS/Shane/RNAseq/20250604_LH00244_0317_A22WFMGLT3_Harding_Shreya_RNAseq/results_rfstranded/deseq2")
  sample_pair_df <- read.delim("C:/PROJECTS/Shane/RNAseq/20250604_LH00244_0317_A22WFMGLT3_Harding_Shreya_RNAseq/config/comparisons.tsv", header = T)
  sample_pair <- apply(sample_pair_df, 1, function(x)paste(x, collapse="-"))
  directions <- c("up", "down")
  comparisons <- c("2Gy24h", "2Gy6d", "10Gy24h", "10Gy6d", "2Gy", "10Gy", "24h", "6d")

  gene_list <- collect_genes(hd, comparisons, sample_pair, directions)
  compareCluster_enrichment(gene_list=gene_list, comparisons=comparisons, c="RSEM", height = 10, width = 12)

}

if(TRUE){
    methylation_list <- "SETD9/SLC25A26/ANTKMT/KMT2B/ATPSCKMT/EEF1AKMT1/SMYD4/METTL6/METTL9/TFB2M/NSUN5/TGS1/COQ3/CARM1/TRMT112/LARP7/TYW3/CARNMT1/SETD4/SETD1A/PRDM4/MTO1/HNMT/SETMAR/NSUN4/MTR/METTL23/TRMT10B/PRDM12/METTL25/PRDM15/NOP2/LCMT2/TRMT1L/METTL15/TRMT61B/SETDB2/SETD2/METTL5/METTL4/DALRD3/ARMT1/EMG1/FDXACB1/NDUFAF5/TRDMT1/TRMT2B/NSUN3/MEPCE/SETD6/PRMT7/MRM1/SNRPD1/EEF2KMT/PRMT3/DIMT1/SPOUT1/PRDM5/WDR4/SMYD3/ZCCHC4/FAM86B1/DNMT1/METTL14/THUMPD2/SUV39H2/FAM86C2P/SNRPD2/METTL21A/BUD23/PRMT6/DNMT3B/TRMT10A/METTL16/RBM15B/SNRPG/EHMT1/MGMT/TYMS/GTPBP3/DPH5/FAM98B/SMYD2/PCIF1/EHMT2/MRM3/COQ5/NSD1/WDR6/TFB1M/THUMPD3/KMT5A/SNRPB/RAMACL/THADA/TRMT9B/N6AMT1/NTMT1/EZH2/SNRPE/COMTD1/EEF1AKMT2/SETDB1/AKT1/PRMT1/SUV39H1/TARBP1/RBM15/NNMT/TRMT11/TRMT5/PEMT/DOT1L/NSD2/FBL/CSKMT/TMT1A" 

    cenpa_list <- "CENPA/H4C3/H4C2/H2AC8/H4C1/H2BC11/H4C4/H4C14/H4C15/H4C8/H2AC4/H4C12/H4C16"

    apoptosis_list <- "EDA2R/CDIP1/MDM2/RRM2B/CDKN1A/TP73/PMAIP1/FBXW7/RPS27L/FYN/WFS1/TREM2/HIPK2/CD24/PHLDA3/APAF1/STYXL1/LRRK2/BAX/CYLD/TNFRSF10B/PLEKHF1/CHEK2/TMEM117/NHERF1/CD44/P4HB/PML/AEN/MOAP1/PDK2/TOPORS/MLH1/FGF2/TRIAP1/EIF2AK3/IER3/HYOU1/IFI27L2/SGMS1/YBX3/POLB/HDAC1/PINK1/RTKN2/CREB3/PTPN1/ITPR1/NKX3-1/BECN1/BBC3/PRKCD/GPX1/MELK/SELENOK/DDIAS/DNAJC10/BCL2L1/SOD1/BRCA1/USP15/JAK2/BRCA2/ATAD5/NLRP1/TXNDC12/NACC2/IFI6/GSKIP/MAPK7/ZNF622/CHAC1"

    dnmt_list <- "DNMT1/DNMT3A/DNMT3B/TET1/TET2/TET3/TDG"

    data_dir <- "C:/PROJECTS/Shane/RNAseq/20250604_LH00244_0317_A22WFMGLT3_Harding_Shreya_RNAseq/results_TE/deseq2"
   
    rsem_abundance <- read.delim(file.path(data_dir, "rsem_gene_level_abundance.tsv"), 
                        header = TRUE, row.names = 1)
                        # Average counts for duplicated GENENAME using group_by and summarize
    rsem_abundance <- rsem_abundance %>%
    group_by(GENENAME) %>%
    summarize(across(everything(), mean), .groups = "drop") %>%
    as.data.frame()

    # Prepare abundance matrix (remove gene name column, keep only count data)
    rownames(rsem_abundance) <- rsem_abundance$GENENAME
    abundance_matrix <- rsem_abundance[, -1]  # Remove GENENAME column

    # Remove leading X from column names
    colnames(abundance_matrix) <- gsub("^X", "", colnames(abundance_matrix))

    NIR_10Gy6d_abundance_matrix <- abundance_matrix[, grepl("NIR|10Gy6d", colnames(abundance_matrix))]
    NIR_2Gy6d_abundance_matrix <- abundance_matrix[, grepl("NIR|2Gy6d", colnames(abundance_matrix))]
    NIR_2Gy24h_abundance_matrix <- abundance_matrix[, grepl("NIR|2Gy24h", colnames(abundance_matrix))]
    NIR_10Gy24h_abundance_matrix <- abundance_matrix[, grepl("NIR|10Gy24h", colnames(abundance_matrix))]
    NIR_10Gy_abundance_matrix <- abundance_matrix[, grepl("NIR|10Gy24h|10Gy6d", colnames(abundance_matrix))]
    NIR_2Gy_abundance_matrix <- abundance_matrix[, grepl("NIR|2Gy24h|2Gy6d", colnames(abundance_matrix))]

    outdir <- "C:/PROJECTS/Shane/RNAseq/20250604_LH00244_0317_A22WFMGLT3_Harding_Shreya_RNAseq/results_TE/deseq2/single_factor_analysis_TE_gene_level"
    if(!dir.exists(outdir)) dir.create(outdir)
    plot_gene_list <- function(gene_list, label_name, abundance_matrix, outdir, folder, width = 6, height = 18){
        gene_list <- unlist(strsplit(gene_list, "/"))
        sub_abundance_matrix <- abundance_matrix[gene_list, ]
        # replace NAs with 0
        sub_abundance_matrix[is.na(sub_abundance_matrix)] <- 0
        # remove rows with all 0s
        sub_abundance_matrix <- sub_abundance_matrix[rowSums(sub_abundance_matrix) > 0, ]
        
        p <- pheatmap(sub_abundance_matrix, cluster_cols = FALSE)
        ggsave(file.path(outdir, folder, paste0(label_name, "_TPM_heatmap.pdf")), plot = p, width = width, height = height)

        sub_abundance_zscore <- t(apply(sub_abundance_matrix, 1, function(x) (x - mean(x)) / sd(x)))
        p <- pheatmap(sub_abundance_zscore, cluster_cols = FALSE)
        
        ggsave(file.path(outdir, folder, paste0(label_name, "_zscore_heatmap.pdf")), plot = p, width = width, height = height)
    }

    plot_gene_list(apoptosis_list, "apoptosis_genes", NIR_10Gy6d_abundance_matrix, outdir, "10Gy6d_vs_NIR")
    plot_gene_list(methylation_list, "methylation_genes", NIR_10Gy6d_abundance_matrix, outdir, "10Gy6d_vs_NIR")
    plot_gene_list(cenpa_list, "cenpa_genes", NIR_10Gy6d_abundance_matrix, outdir, "10Gy6d_vs_NIR")
 
    plot_gene_list(apoptosis_list, "apoptosis_genes", NIR_10Gy24h_abundance_matrix, outdir, "10Gy24h_vs_NIR")
    plot_gene_list(methylation_list, "methylation_genes", NIR_10Gy24h_abundance_matrix, outdir, "10Gy24h_vs_NIR")

    plot_gene_list(dnmt_list, "dnmt_genes", abundance_matrix, data_dir, "single_factor_analysis_TE_gene_level", width = 12, height = 8)
    plot_gene_list(dnmt_list, "dnmt_genes_10Gy", NIR_10Gy_abundance_matrix, data_dir, "single_factor_analysis_TE_gene_level", width = 12, height = 8)
    plot_gene_list(dnmt_list, "dnmt_genes_2Gy", NIR_2Gy_abundance_matrix, data_dir, "single_factor_analysis_TE_gene_level", width = 12, height = 8)

    # Plot lineplot for DNMTs with t-test and significance
    gene_list <- unlist(strsplit(dnmt_list, "/"))
    sub_abundance_df <- as.data.frame(NIR_10Gy_abundance_matrix[gene_list, ])
    # replace NAs with 0
    sub_abundance_df[is.na(sub_abundance_df)] <- 0
    # remove rows with all 0s
    sub_abundance_df <- sub_abundance_df[rowSums(sub_abundance_df) > 0, ]
    
    sub_abundance_zscore <- as.data.frame(t(apply(sub_abundance_df, 1, function(x) (x - mean(x)) / sd(x)))) %>%
        mutate(gene=rownames(sub_abundance_df))

    sub_abundance_zscore_long <- pivot_longer(sub_abundance_zscore, cols = -gene, names_to = "sample", values_to = "zscore") %>%
        mutate(condition = factor(sapply(sample, function(x)strsplit(x, "_")[[1]][1]), levels = c("NIR", "10Gy24h", "10Gy6d")))

    sub_abundance_df <- sub_abundance_df %>%
        mutate(gene=rownames(sub_abundance_df))
    sub_abundance_long <- pivot_longer(sub_abundance_df, cols = -gene, names_to = "sample", values_to = "TPM") %>%
        mutate(condition = factor(sapply(sample, function(x)strsplit(x, "_")[[1]][1]), levels = c("NIR", "10Gy24h", "10Gy6d")))

    # Perform Fisher's LSD test for each gene
    library(agricolae)
    lsd_results <- data.frame()
    anova_pvalues <- data.frame()
    for(g in unique(sub_abundance_long$gene)){
        gene_data <- sub_abundance_long %>% filter(gene == g)
        
        # Perform ANOVA first
        anova_model <- aov(TPM ~ condition, data = gene_data)
        anova_result <- anova(anova_model)
        anova_pvalue <- anova_result$`Pr(>F)`[1]
        
        # Store ANOVA p-value for each gene
        anova_pvalues <- rbind(anova_pvalues, data.frame(
            gene = g,
            anova_pvalue = anova_pvalue,
            anova_label = paste0("ANOVA p = ", format(anova_pvalue, digits = 3, scientific = TRUE))
        ))

        # Perform Fisher's LSD post-hoc test
        lsd_test <- LSD.test(anova_model, "condition", p.adj = "none", group = FALSE)
        lsd_comparisons <- lsd_test$comparison
        
        # Extract NIR vs 10Gy24h comparison
        comp_24h_nir <- lsd_comparisons[rownames(lsd_comparisons) == "10Gy24h - NIR", ]
        if(length(comp_24h_nir) > 0){
            lsd_results <- rbind(lsd_results, data.frame(
                gene = g,
                condition = "10Gy24h",
                pvalue = comp_24h_nir$pvalue,
                significance = ifelse(comp_24h_nir$pvalue < 0.001, "***", 
                                    ifelse(comp_24h_nir$pvalue < 0.01, "**",
                                         ifelse(comp_24h_nir$pvalue < 0.05, "*", "ns")))
            ))
        }
        
        # Extract NIR vs 10Gy6d comparison
        comp_6d_nir <- lsd_comparisons[rownames(lsd_comparisons) == "10Gy6d - NIR", ]
        if(length(comp_6d_nir) > 0){
            lsd_results <- rbind(lsd_results, data.frame(
                gene = g,
                condition = "10Gy6d",
                pvalue = comp_6d_nir$pvalue,
                significance = ifelse(comp_6d_nir$pvalue < 0.001, "***", 
                                    ifelse(comp_6d_nir$pvalue < 0.01, "**",
                                         ifelse(comp_6d_nir$pvalue < 0.05, "*", "ns")))
            ))
        }
    }

    # Add y-position for significance labels
    summary_stats <- sub_abundance_long %>%
        group_by(gene, condition) %>%
        summarise(mean_TPM = mean(TPM), 
                  se = sd(TPM)/sqrt(dplyr::n()),
                  .groups = "drop")
    
    max_y <- summary_stats %>%
        group_by(gene) %>%
        summarise(max_y = max(mean_TPM + se), .groups = "drop")
    
    min_y <- summary_stats %>%
        group_by(gene) %>%
        summarise(min_y = min(mean_TPM - se), .groups = "drop")
    
    lsd_results <- lsd_results %>%
        left_join(max_y, by = "gene") %>%
        mutate(y_pos = max_y + 0.2)
    
    # Add position for ANOVA labels at bottom of each facet
    anova_pvalues <- anova_pvalues %>%
        left_join(min_y, by = "gene") %>%
        mutate(x_pos = 2,  # Middle of three conditions
               y_pos = min_y - 0.5)

    p <- ggplot(sub_abundance_long, aes(x = condition, y = TPM)) +
        geom_line(stat = "summary", fun = "mean", aes(group = 1)) +
        geom_point(stat = "summary", fun = "mean", aes(group = 1), size = 3) +
        geom_errorbar(stat = "summary", fun.data = "mean_se", aes(group = 1), width = 0.1) +
        geom_text(data = lsd_results, aes(x = condition, y = y_pos, label = significance), 
                  size = 5, vjust = 0) +
        geom_text(data = anova_pvalues, aes(x = x_pos, y = y_pos, label = anova_label),
                  size = 3, hjust = 0.5, vjust = 1, fontface = "italic") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        labs(x = "Condition", y = "TPM", title = "DNMT Gene Expression (Fisher's LSD, n=3)") +
        facet_wrap(~ gene, scales = "free_y")
    
    ggsave(file.path(data_dir, "single_factor_analysis_TE_gene_level", paste0("dnmt_genes_10Gy_lineplot.pdf")), 
           plot = p, width = 12, height = 12)

}
