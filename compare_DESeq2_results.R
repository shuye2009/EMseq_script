library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)

rnaseq_base_dir   <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/RNAseq/deseq2/single_factor_analysis_gene_level"

# Output Path
output_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/RNAseq/compare_DEG_GSEA"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Get list of comparisons from rnaseq directory
comparisons <- list.dirs(rnaseq_base_dir, full.names = FALSE, recursive = FALSE)
# Filter for actual comparison folders (contain "_vs_")
comparisons <- comparisons[grepl("_vs_", comparisons)]
# Filter out files or non-comparison folders if any (simple check: must be in targeted dir too usually)
comparisons <- comparisons[dir.exists(file.path(rnaseq_base_dir, comparisons))]

message("Found ", length(comparisons), " comparisons: ", paste(comparisons, collapse = ", "))

# Read DESeq2 results for each comparison
comparison_data <- list()
for (comp in comparisons) {
  comp_dir <- file.path(rnaseq_base_dir, comp)
  comp_file <- list.files(comp_dir, pattern = "DESeq2_results_regular.csv", full.names = TRUE)
  if (length(comp_file) > 0) {
    comparison_data[[comp]] <- as.data.frame(read.csv(comp_file))
  }
}

# Process 10Gy24h data
df_24h <- comparison_data[["10Gy24h_vs_NIR"]] %>%
  filter(!is.na(padj)) %>%
  mutate(rank = stat) %>%
  arrange(desc(rank))

ranked_genes_24h <- df_24h$rank
names(ranked_genes_24h) <- df_24h$gene_name

# Process 10Gy6d data
df_6d <- comparison_data[["10Gy6d_vs_NIR"]] %>%
  filter(!is.na(padj)) %>%
  mutate(rank = stat) %>%
  arrange(desc(rank))

ranked_genes_6d <- df_6d$rank
names(ranked_genes_6d) <- df_6d$gene_name

message("\nRunning GSEA for 10Gy24h...")
# Run GSEA for each timepoint
gsea_24h <- gseGO(geneList = ranked_genes_24h,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pvalueCutoff = 0.05,
                  keyType = "SYMBOL")

message("\nRunning GSEA for 10Gy6d...")
gsea_6d <- gseGO(geneList = ranked_genes_6d,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 pvalueCutoff = 0.05,
                 keyType = "SYMBOL")

message("\nMerging results for comparison...")
# Merge for comparison
comparison <- merge_result(list("10Gy_24h" = gsea_24h, 
                                "10Gy_6d" = gsea_6d))

# Visualize comparison
message("\nGenerating comparison plots...")
p1 <- dotplot(comparison, showCategory = 30, split = ".sign", label_format = 100) + 
  facet_grid(.~.sign)
ggsave(file.path(output_dir, "GSEA_10Gy6d_vs_10Gy24h_dotplot.pdf"), p1, width = 20, height = 20)

# Save results
write.csv(as.data.frame(gsea_24h), file.path(output_dir, "GSEA_10Gy24h_vs_NIR_results.csv"), row.names = FALSE)
write.csv(as.data.frame(gsea_6d), file.path(output_dir, "GSEA_10Gy6d_vs_NIR_results.csv"), row.names = FALSE)

message("\nDone! Results saved to: ", output_dir)

# compare difference in ranks of common genes
common_genes <- intersect(names(ranked_genes_24h), names(ranked_genes_6d))
rank_diff <- ranked_genes_6d[common_genes] - ranked_genes_24h[common_genes]
# sort rank_diff in descending order
rank_diff <- sort(rank_diff, decreasing = TRUE)
# GSEA on rank_diff
gsea_rank_diff <- gseGO(geneList = rank_diff,
                        OrgDb = org.Hs.eg.db,
                        ont = "BP",
                        pvalueCutoff = 0.05,
                        keyType = "SYMBOL")

# Filter for top 50 activated and top 50 suppressed by q-value
gsea_rank_diff_df <- as.data.frame(gsea_rank_diff)

# Get top 50 activated (positive NES) by q-value
activated <- gsea_rank_diff_df %>%
  filter(NES > 0) %>%
  arrange(p.adjust) %>%
  head(50)

# Get top 50 suppressed (negative NES) by q-value
suppressed <- gsea_rank_diff_df %>%
  filter(NES < 0) %>%
  arrange(p.adjust) %>%
  head(50)

# Combine top pathways
top_pathways <- rbind(activated, suppressed)
selected_ids <- top_pathways$ID

# Create filtered GSEA object for plotting
gsea_rank_diff_filtered <- gsea_rank_diff
gsea_rank_diff_filtered@result <- gsea_rank_diff@result[gsea_rank_diff@result$ID %in% selected_ids, ]

# Plot with filtered results
message("\nPlotting top 50 activated and top 50 suppressed pathways...")
p2 <- ridgeplot(gsea_rank_diff_filtered, showCategory = 100, label_format = 100) +
  ggtitle("Top 50 Activated and Top 50 Suppressed Pathways (by q-value)\n10Gy 6d vs 24h")
ggsave(file.path(output_dir, "GSEA_10Gy6d_vs_10Gy24h_rank_diff_top100_ridgeplot.pdf"), p2, width = 20, height = 20)


# Save all gsea_rank_diff results
write.csv(as.data.frame(gsea_rank_diff), file.path(output_dir, "GSEA_10Gy6d_vs_10Gy24h_rank_diff_ALL_results.csv"), row.names = FALSE)
# Save filtered top 100 results
write.csv(top_pathways, file.path(output_dir, "GSEA_10Gy6d_vs_10Gy24h_rank_diff_TOP100_results.csv"), row.names = FALSE)


