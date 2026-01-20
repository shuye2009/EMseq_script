rm(list = ls())

library(GenomicRanges)
library(ChIPseeker)
library(org.Hs.eg.db)
library(readr)
library(openxlsx)
library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)

gtf_file <- "C:/PROJECTS/resource/T2T_CHM13/Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf"
txdb <- txdbmaker::makeTxDbFromGFF(gtf_file,
                        format = "gtf",
                        organism = "Homo sapiens")

gtf_data <- rtracklayer::import(gtf_file)
gene_mapping <- as.data.frame(gtf_data) %>%
    dplyr::filter(type == "gene") %>%
    dplyr::select(gene_id, gene_name) %>%
    dplyr::distinct()
message("Created gene mapping with ", nrow(gene_mapping), " genes")

expression_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/RNAseq/deseq2/single_factor_analysis_gene_level"

dmr_method <- "dmrseq"

emseq_dir <- file.path("C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq", dmr_method)

dmr_folder <- "DMR"

sample_dirs <- list.dirs(emseq_dir, full.names = FALSE, recursive = FALSE)

names(sample_dirs) <- sample_dirs

sample_names <- c("IR2Gy24h_vs_NIR", "IR10Gy24h_vs_NIR", "IR10Gy24h_vs_IR2Gy24h", 
                  "IR2Gy6d_vs_NIR", "IR10Gy6d_vs_NIR", "IR10Gy6d_vs_IR2Gy6d", 
                  "IR2Gy6d_vs_IR2Gy24h", "IR10Gy6d_vs_IR10Gy24h")

sample_dirs <- sample_dirs[sample_names] 

enrichment_results <- data.frame()

for(sample_name in sample_names){
    message("\n============================================================")
    message("Processing ", sample_name)
    message("============================================================")
    
    SA <- gsub("^IR", "", sample_name)
    SA <- gsub("_IR", "_", SA)
    
    working_dir <- file.path(emseq_dir, sample_dirs[sample_name], dmr_folder)
    output_dir <- file.path(working_dir, "Integrate_gene_expression_with_DMR")
    if(!dir.exists(output_dir)) dir.create(output_dir)
    
    dge <- read.csv(paste0(expression_dir, "/",SA, "/", SA, "_DESeq2_results_regular.csv")) %>%
        dplyr::filter(!is.na(padj)) %>%
        mutate(regulation = case_when(
            log2FoldChange > 1  & padj < 0.05 ~ "Up-regulated",
            log2FoldChange < -1 & padj < 0.05 ~ "Down-regulated",
            TRUE ~ "No-change"
        )) %>%
        dplyr::select(-direction) 
    
    dmr <- read_tsv(file.path(working_dir, "DMR_table.tsv")) %>%
        mutate(direction = case_when(
            stat > 0 & pval < 0.05 ~ "Hypermethylated", 
            stat < 0 & pval < 0.05 ~ "Hypomethylated",
            TRUE ~ "Not_significant"))

    dmr_gr <- makeGRangesFromDataFrame(dmr, keep.extra.columns = TRUE)
    dmrAnno <- annotatePeak(dmr_gr, 
                         tssRegion=c(-1000, 0),
                         TxDb=txdb, 
                         annoDb="org.Hs.eg.db")

    annotated_dmr <- as.data.frame(dmrAnno) %>%
        mutate(annotation = gsub(" \\(.+\\)", "", annotation)) %>%
        filter(annotation != "Distal Intergenic")
    
    annotated_dmr <- annotated_dmr %>%
        left_join(gene_mapping, by = c("geneId" = "gene_id"))
    
    message("Mapped ", sum(!is.na(annotated_dmr$gene_name)), " out of ", nrow(annotated_dmr), " DMRs to gene symbols")
    
    gene_dmr_overlap <- full_join(dge, annotated_dmr, by = "gene_name")
    
    gene_dmr_summary <- gene_dmr_overlap %>%
        filter(!is.na(seqnames), regulation %in% c("Up-regulated", "Down-regulated"), direction %in% c("Hypomethylated", "Hypermethylated")) %>%
        group_by(gene_name, regulation) %>%
        summarise(
            n_hypo = sum(direction == "Hypomethylated"),
            n_hyper = sum(direction == "Hypermethylated"),
            .groups = "drop"
        ) %>%
        mutate(methylation_status = case_when(
            n_hypo > n_hyper ~ "Hypomethylated",
            n_hyper > n_hypo ~ "Hypermethylated",
            n_hypo == n_hyper ~ "Mixed"
        )) %>%
        mutate(comparison = sample_name)
    
    hypomethylated_genes <- gene_dmr_summary %>%
        filter(methylation_status == "Hypomethylated") %>%
        pull(gene_name) %>%
        unique()
    
    hypermethylated_genes <- gene_dmr_summary %>%
        filter(methylation_status == "Hypermethylated") %>%
        pull(gene_name) %>%
        unique()
    
    mixed_genes <- gene_dmr_summary %>%
        filter(methylation_status == "Mixed") %>%
        pull(gene_name) %>%
        unique()
    
    upregulated_genes <- gene_dmr_summary %>% 
        filter(regulation == "Up-regulated",
               methylation_status %in% c("Hypomethylated", "Hypermethylated")) %>%
        pull(gene_name) %>%
        unique()
    
    downregulated_genes <- gene_dmr_summary %>% 
        filter(regulation == "Down-regulated",
               methylation_status %in% c("Hypomethylated", "Hypermethylated")) %>%
        pull(gene_name) %>%
        unique()
    
    total_genes <- length(upregulated_genes) + length(downregulated_genes)
    message("\nTotal genes: ", total_genes)
    message("Up-regulated genes (in hypo/hyper DMRs): ", length(upregulated_genes))
    message("Down-regulated genes (in hypo/hyper DMRs): ", length(downregulated_genes))
    message("Genes predominantly in hypomethylated DMRs: ", length(hypomethylated_genes))
    message("Genes predominantly in hypermethylated DMRs: ", length(hypermethylated_genes))
    message("Genes with equal hypo/hyper DMRs (mixed): ", length(mixed_genes))
    
    up_in_hypo <- length(intersect(upregulated_genes, hypomethylated_genes))
    up_not_hypo <- length(upregulated_genes) - up_in_hypo
    not_up_in_hypo <- length(hypomethylated_genes) - up_in_hypo
    not_up_not_hypo <- total_genes - length(upregulated_genes) - length(hypomethylated_genes) + up_in_hypo
    
    message("\n--- Upregulated vs Hypomethylated ---")
    message("Up-regulated & Hypomethylated: ", up_in_hypo)
    message("Up-regulated & NOT Hypomethylated: ", up_not_hypo)
    message("NOT Up-regulated & Hypomethylated: ", not_up_in_hypo)
    message("NOT Up-regulated & NOT Hypomethylated: ", not_up_not_hypo)
    
    if(any(c(up_in_hypo, up_not_hypo, not_up_in_hypo, not_up_not_hypo) < 0)) {
        message("ERROR: Negative values in contingency table!")
        message("Skipping Fisher test for this comparison")
        next
    }
    
    contingency_matrix_up_hypo <- matrix(c(up_in_hypo, up_not_hypo, 
                                            not_up_in_hypo, not_up_not_hypo),
                                          nrow = 2, byrow = TRUE,
                                          dimnames = list(
                                              c("Upregulated", "Not_Upregulated"),
                                              c("Hypomethylated", "Not_Hypomethylated")
                                          ))
    
    fisher_test_up_hypo <- fisher.test(contingency_matrix_up_hypo)
    
    message("\nFisher's Exact Test (Upregulated vs Hypomethylated):")
    message("Odds Ratio: ", round(fisher_test_up_hypo$estimate, 3))
    message("P-value: ", format(fisher_test_up_hypo$p.value, scientific = TRUE, digits = 3))
    message("95% CI: [", round(fisher_test_up_hypo$conf.int[1], 3), ", ", 
            round(fisher_test_up_hypo$conf.int[2], 3), "]")
    
    down_in_hypo <- length(intersect(downregulated_genes, hypomethylated_genes))
    down_not_hypo <- length(downregulated_genes) - down_in_hypo
    not_down_in_hypo <- length(hypomethylated_genes) - down_in_hypo
    not_down_not_hypo <- total_genes - length(downregulated_genes) - length(hypomethylated_genes) + down_in_hypo
    
    message("\n--- Downregulated vs Hypomethylated ---")
    message("Down-regulated & Hypomethylated: ", down_in_hypo)
    
    if(any(c(down_in_hypo, down_not_hypo, not_down_in_hypo, not_down_not_hypo) < 0)) {
        message("ERROR: Negative values in contingency table!")
        message("Skipping Fisher test for this comparison")
        next
    }
    
    contingency_matrix_down_hypo <- matrix(c(down_in_hypo, down_not_hypo, 
                                              not_down_in_hypo, not_down_not_hypo),
                                            nrow = 2, byrow = TRUE)
    
    fisher_test_down_hypo <- fisher.test(contingency_matrix_down_hypo)
    
    up_in_hyper <- length(intersect(upregulated_genes, hypermethylated_genes))
    up_not_hyper <- length(upregulated_genes) - up_in_hyper
    not_up_in_hyper <- length(hypermethylated_genes) - up_in_hyper
    not_up_not_hyper <- total_genes - length(upregulated_genes) - length(hypermethylated_genes) + up_in_hyper
    
    message("\n--- Upregulated vs Hypermethylated ---")
    message("Up-regulated & Hypermethylated: ", up_in_hyper)
    
    if(any(c(up_in_hyper, up_not_hyper, not_up_in_hyper, not_up_not_hyper) < 0)) {
        message("ERROR: Negative values in contingency table!")
        message("Skipping Fisher test for this comparison")
        next
    }
    
    contingency_matrix_up_hyper <- matrix(c(up_in_hyper, up_not_hyper, 
                                             not_up_in_hyper, not_up_not_hyper),
                                           nrow = 2, byrow = TRUE)
    
    fisher_test_up_hyper <- fisher.test(contingency_matrix_up_hyper)
    
    down_in_hyper <- length(intersect(downregulated_genes, hypermethylated_genes))
    down_not_hyper <- length(downregulated_genes) - down_in_hyper
    not_down_in_hyper <- length(hypermethylated_genes) - down_in_hyper
    not_down_not_hyper <- total_genes - length(downregulated_genes) - length(hypermethylated_genes) + down_in_hyper
    
    message("\n--- Downregulated vs Hypermethylated ---")
    message("Down-regulated & Hypermethylated: ", down_in_hyper)
    
    if(any(c(down_in_hyper, down_not_hyper, not_down_in_hyper, not_down_not_hyper) < 0)) {
        message("ERROR: Negative values in contingency table!")
        message("Skipping Fisher test for this comparison")
        next
    }
    
    contingency_matrix_down_hyper <- matrix(c(down_in_hyper, down_not_hyper, 
                                               not_down_in_hyper, not_down_not_hyper),
                                             nrow = 2, byrow = TRUE)
    
    fisher_test_down_hyper <- fisher.test(contingency_matrix_down_hyper)
    
    result_row <- data.frame(
        comparison = sample_name,
        total_genes = total_genes,
        upregulated_genes = length(upregulated_genes),
        downregulated_genes = length(downregulated_genes),
        hypomethylated_genes = length(hypomethylated_genes),
        hypermethylated_genes = length(hypermethylated_genes),
        up_in_hypo = up_in_hypo,
        up_in_hypo_odds_ratio = fisher_test_up_hypo$estimate,
        up_in_hypo_pvalue = fisher_test_up_hypo$p.value,
        up_in_hypo_CI_lower = fisher_test_up_hypo$conf.int[1],
        up_in_hypo_CI_upper = fisher_test_up_hypo$conf.int[2],
        down_in_hypo = down_in_hypo,
        down_in_hypo_odds_ratio = fisher_test_down_hypo$estimate,
        down_in_hypo_pvalue = fisher_test_down_hypo$p.value,
        up_in_hyper = up_in_hyper,
        up_in_hyper_odds_ratio = fisher_test_up_hyper$estimate,
        up_in_hyper_pvalue = fisher_test_up_hyper$p.value,
        down_in_hyper = down_in_hyper,
        down_in_hyper_odds_ratio = fisher_test_down_hyper$estimate,
        down_in_hyper_pvalue = fisher_test_down_hyper$p.value
    )
    
    enrichment_results <- rbind(enrichment_results, result_row)
    
    test_results <- data.frame(
        Test = c("Up in Hypomethylated", "Down in Hypomethylated", 
                 "Up in Hypermethylated", "Down in Hypermethylated"),
        Overlap = c(up_in_hypo, down_in_hypo, up_in_hyper, down_in_hyper),
        Odds_Ratio = c(fisher_test_up_hypo$estimate, fisher_test_down_hypo$estimate,
                       fisher_test_up_hyper$estimate, fisher_test_down_hyper$estimate),
        P_value = c(fisher_test_up_hypo$p.value, fisher_test_down_hypo$p.value,
                    fisher_test_up_hyper$p.value, fisher_test_down_hyper$p.value),
        CI_lower = c(fisher_test_up_hypo$conf.int[1], fisher_test_down_hypo$conf.int[1],
                     fisher_test_up_hyper$conf.int[1], fisher_test_down_hyper$conf.int[1]),
        CI_upper = c(fisher_test_up_hypo$conf.int[2], fisher_test_down_hypo$conf.int[2],
                     fisher_test_up_hyper$conf.int[2], fisher_test_down_hyper$conf.int[2])
    )
    
    write.xlsx(test_results, paste0(output_dir, "/enrichment_test_results.xlsx"))
    
    write.xlsx(list(
        "Upregulated_Hypomethylated" = as.data.frame(contingency_matrix_up_hypo),
        "Downregulated_Hypomethylated" = as.data.frame(contingency_matrix_down_hypo),
        "Upregulated_Hypermethylated" = as.data.frame(contingency_matrix_up_hyper),
        "Downregulated_Hypermethylated" = as.data.frame(contingency_matrix_down_hyper)
    ), paste0(output_dir, "/contingency_tables.xlsx"))
    
    upregulated_in_hypomethylated <- data.frame(
        gene_name = intersect(upregulated_genes, hypomethylated_genes)
    ) %>%
        left_join(dge %>% dplyr::select(gene_name, log2FoldChange, padj), by = "gene_name") %>%
        arrange(desc(log2FoldChange))
    
    write.xlsx(upregulated_in_hypomethylated, 
               paste0(output_dir, "/upregulated_genes_in_hypomethylated_DMRs.xlsx"))
}

write.xlsx(enrichment_results, paste0(emseq_dir, "/enrichment_test_summary_all_comparisons.xlsx"))
write.csv(enrichment_results, paste0(emseq_dir, "/enrichment_test_summary_all_comparisons.csv"), row.names = FALSE)

enrichment_results_long <- enrichment_results %>%
    dplyr::select(comparison, up_in_hypo_odds_ratio, down_in_hypo_odds_ratio, 
           up_in_hyper_odds_ratio, down_in_hyper_odds_ratio,
           up_in_hypo_pvalue, down_in_hypo_pvalue,
           up_in_hyper_pvalue, down_in_hyper_pvalue) %>%
    pivot_longer(cols = -comparison, 
                 names_to = c("test", ".value"),
                 names_pattern = "(.+)_(odds_ratio|pvalue)") %>%
    mutate(
        significance = case_when(
            pvalue < 0.001 ~ "***",
            pvalue < 0.01 ~ "**",
            pvalue < 0.05 ~ "*",
            TRUE ~ "ns"
        ),
        test_label = case_when(
            test == "up_in_hypo" ~ "Up in Hypomethylated",
            test == "down_in_hypo" ~ "Down in Hypomethylated",
            test == "up_in_hyper" ~ "Up in Hypermethylated",
            test == "down_in_hyper" ~ "Down in Hypermethylated"
        ),
        comparison = factor(comparison, levels = sample_names)
    )

p_odds_ratio <- ggplot(enrichment_results_long, 
                       aes(x = comparison, y = odds_ratio, fill = test_label)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = significance), 
              position = position_dodge(width = 0.8), 
              vjust = -0.5, size = 5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("Up in Hypomethylated" = "#E74C3C", 
                                  "Down in Hypomethylated" = "#3498DB",
                                  "Up in Hypermethylated" = "#F39C12",
                                  "Down in Hypermethylated" = "#9B59B6")) +
    labs(title = "Enrichment Analysis: Gene Expression vs Methylation",
         x = "Comparison",
         y = "Odds Ratio",
         fill = "Test") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          legend.position = "top",
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 14))

ggsave(paste0(emseq_dir, "/enrichment_odds_ratio_barplot.pdf"), 
       plot = p_odds_ratio, width = 16, height = 8)
ggsave(paste0(emseq_dir, "/enrichment_odds_ratio_barplot.png"), 
       plot = p_odds_ratio, width = 16, height = 8, dpi = 300)

p_log_odds <- ggplot(enrichment_results_long, 
                     aes(x = comparison, y = log2(odds_ratio), fill = test_label)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = significance), 
              position = position_dodge(width = 0.8), 
              vjust = ifelse(enrichment_results_long$odds_ratio > 1, -0.5, 1.5), 
              size = 5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("Up in Hypomethylated" = "#E74C3C", 
                                  "Down in Hypomethylated" = "#3498DB",
                                  "Up in Hypermethylated" = "#F39C12",
                                  "Down in Hypermethylated" = "#9B59B6")) +
    labs(title = "Enrichment Analysis: Gene Expression vs Methylation (Log2 Scale)",
         x = "Comparison",
         y = "Log2(Odds Ratio)",
         fill = "Test") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          legend.position = "top",
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 14))

ggsave(paste0(emseq_dir, "/enrichment_log2_odds_ratio_barplot.pdf"), 
       plot = p_log_odds, width = 16, height = 8)
ggsave(paste0(emseq_dir, "/enrichment_log2_odds_ratio_barplot.png"), 
       plot = p_log_odds, width = 16, height = 8, dpi = 300)

message("\n========================================")
message("ENRICHMENT ANALYSIS COMPLETE")
message("========================================")
message("Results saved to: ", emseq_dir)
message("\nKey files:")
message("- enrichment_test_summary_all_comparisons.xlsx")
message("- enrichment_odds_ratio_barplot.pdf")
message("- enrichment_log2_odds_ratio_barplot.pdf")
