rm(list = ls())

# ==============================================================================
# SASP Analysis Script
# Analyzes expression of senescence-associated secretory phenotype (SASP) genes
# in RNA-seq samples: NIR vs IR groups (2Gy24h, 2Gy6d, 10Gy24h, 10Gy6d)
# ==============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(singscore)
library(SummarizedExperiment)
library(DESeq2)
library(openxlsx)
library(org.Hs.eg.db)

# ==============================================================================
# Paths
# ==============================================================================
rnaseq_base_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/RNAseq/deseq2"
design_file <- "C:/PROJECTS/Shane/RNAseq/20250604_LH00244_0317_A22WFMGLT3_Harding_Shreya_RNAseq/config/design.tsv"
deseq2_dir <- file.path(rnaseq_base_dir, "single_factor_analysis_gene_level")
outdir <- file.path(rnaseq_base_dir, "SASP_analysis")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ==============================================================================
# Define SASP gene list (curated from literature: SenMayo, SASP Atlas, Reactome)
# ==============================================================================
sasp_genes <- c(
  # Pro-inflammatory cytokines
  "IL1A", "IL1B", "IL6", "TNF", "IFNG",
  # Chemokines (CXC family)
  "CXCL8", "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL10", "CXCL12",
  # Chemokines (CC family)
  "CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8", "CCL13", "CCL20", "CCL26",
  # Growth factors & regulators
  "VEGFA", "FGF2", "IGF1", "TGFB1", "AREG", "CSF1", "CSF2", "CSF3",
  "HGF", "EGF", "PDGFA", "PDGFB", "IGFBP3", "IGFBP7",
  # Matrix metalloproteinases
  "MMP1", "MMP2", "MMP3", "MMP9", "MMP10", "MMP12", "MMP13", "MMP14",
  # TIMPs
  "TIMP1", "TIMP2",
  # Serine proteases & regulators
  "PLAU", "PLAT", "SERPINE1", "SERPINB2",
  # Cathepsins
  "CTSB", "CTSD", "CTSK",
  # Adhesion molecules
  "ICAM1", "VCAM1",
  # Other secreted factors
  "GDF15", "STC1", "PTGS2", "SPARC", "CHI3L1", "LGALS3", "CLU", "FSTL1",
  # Senescence markers (non-secreted, for reference)
  "CDKN1A", "CDKN2A", "LMNB1", "TP53"
)

cat("Total SASP genes defined:", length(sasp_genes), "\n")

# ==============================================================================
# Load sample design
# ==============================================================================
sample_design <- read.delim(design_file, header = TRUE)
rownames(sample_design) <- sample_design$sample
sample_design$condition <- factor(sample_design$condition,
                                   levels = c("NIR", "2Gy24h", "2Gy6d", "10Gy24h", "10Gy6d"))
sample_design$ir_group <- ifelse(sample_design$condition == "NIR", "NIR", "IR")
cat("\nSample design:\n")
print(sample_design)

# ==============================================================================
# Load TPM expression data
# ==============================================================================
tpm_raw <- read.csv(file.path(rnaseq_base_dir, "rsem_gene_level_abundance.tsv"),
                     row.names = 1, sep = "\t", check.names = FALSE)

gene_symbols <- tpm_raw[, 1]
tpm_matrix <- as.matrix(tpm_raw[, -1])
rownames(tpm_matrix) <- gene_symbols

cat("\nTPM matrix dimensions:", dim(tpm_matrix), "\n")

# Log2-transform TPM
log2tpm <- log2(tpm_matrix + 1)

# ==============================================================================
# Filter for SASP genes present in the data
# ==============================================================================
sasp_found <- intersect(sasp_genes, rownames(log2tpm))
sasp_missing <- setdiff(sasp_genes, rownames(log2tpm))

cat("\nSASP genes found in expression data:", length(sasp_found), "/", length(sasp_genes), "\n")
if (length(sasp_missing) > 0) {
  cat("Missing SASP genes:", paste(sasp_missing, collapse = ", "), "\n")
}

sasp_log2tpm <- log2tpm[sasp_found, , drop = FALSE]
sasp_tpm <- tpm_matrix[sasp_found, , drop = FALSE]

# ==============================================================================
# Load DESeq2 results for all IR vs NIR contrasts
# ==============================================================================
contrasts <- c("2Gy24h_vs_NIR", "2Gy6d_vs_NIR", "10Gy24h_vs_NIR", "10Gy6d_vs_NIR")
deseq2_list <- list()

for (contr in contrasts) {
  res_file <- file.path(deseq2_dir, contr, paste0(contr, "_DESeq2_results_regular.csv"))
  if (file.exists(res_file)) {
    res <- read.csv(res_file)
    res$contrast <- contr
    deseq2_list[[contr]] <- res
    cat("Loaded DESeq2 results for:", contr, "-", nrow(res), "genes\n")
  } else {
    cat("WARNING: DESeq2 results not found for:", contr, "\n")
  }
}

# Combine all DESeq2 results
deseq2_all <- bind_rows(deseq2_list)
deseq2_all$contrast <- factor(deseq2_all$contrast, levels = contrasts)

# Filter for SASP genes
deseq2_sasp <- deseq2_all %>% filter(gene_name %in% sasp_found)
cat("\nSASP genes with DESeq2 results:", length(unique(deseq2_sasp$gene_name)), "\n")

# ==============================================================================
# 1. DESeq2 Stats Plots
# ==============================================================================

# 1a. Volcano plots for each contrast, highlighting SASP genes
for (contr in contrasts) {
  res_contr <- deseq2_list[[contr]]
  if (is.null(res_contr)) next
  
  res_contr <- res_contr %>%
    mutate(is_sasp = gene_name %in% sasp_found,
           sig_sasp = is_sasp & padj < 0.05 & abs(log2FoldChange) > 1,
           neglog10_padj = -log10(padj))
  
  # Cap for plotting
  res_contr$neglog10_padj <- pmin(res_contr$neglog10_padj, 50)
  
  top_labels <- res_contr %>%
    filter(sig_sasp) %>%
    arrange(padj) %>%
    head(20)
  
  p <- ggplot(res_contr, aes(x = log2FoldChange, y = neglog10_padj)) +
    geom_point(aes(color = sig_sasp), size = 0.8, alpha = 0.6) +
    geom_point(data = res_contr %>% filter(sig_sasp),
               aes(x = log2FoldChange, y = neglog10_padj),
               color = "#E41A1C", size = 1.5) +
    geom_text_repel(data = top_labels,
                    aes(label = gene_name),
                    size = 3, max.overlaps = 20, color = "#E41A1C") +
    scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey70")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
    labs(title = paste(contr, "- SASP Genes Highlighted"),
         subtitle = paste0("Red: SASP genes (padj<0.05 & |log2FC|>1), n=",
                          sum(res_contr$sig_sasp, na.rm = TRUE)),
         x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  
  ggsave(file.path(outdir, paste0(contr, "_volcano_SASP.pdf")),
         p, width = 10, height = 8, dpi = 300)
  cat("Saved volcano plot:", contr, "\n")
}

# 1b. Combined log2FC bar plot of SASP genes across all contrasts
sasp_fc_summary <- deseq2_sasp %>%
  filter(padj < 0.05) %>%
  dplyr::select(gene_name, log2FoldChange, contrast) %>%
  pivot_wider(names_from = contrast, values_from = log2FoldChange, values_fill = 0) %>%
  pivot_longer(-gene_name, names_to = "contrast", values_to = "log2FoldChange") %>%
  mutate(contrast = factor(contrast, levels = contrasts))

# Heatmap-style dot plot of SASP gene log2FC
p_fc <- ggplot(sasp_fc_summary, aes(x = contrast, y = gene_name, fill = log2FoldChange)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "log2FC") +
  labs(title = "SASP Genes: log2 Fold Change (IR vs NIR)",
       subtitle = "Only genes with padj < 0.05 in at least one contrast",
       x = "", y = "") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9))

ggsave(file.path(outdir, "SASP_log2FC_heatmap.pdf"),
       p_fc, width = 8, height = max(8, length(unique(sasp_fc_summary$gene_name)) * 0.3),
       dpi = 300)

# 1c. Bar plot of number of significant SASP genes per contrast
sig_counts <- deseq2_sasp %>%
  filter(padj < 0.05) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
  dplyr::count(contrast, direction)

p_bar <- ggplot(sig_counts, aes(x = contrast, y = n, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(0.7), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Up" = "#B2182B", "Down" = "#2166AC")) +
  labs(title = "Significant SASP Genes per Contrast",
       subtitle = "padj < 0.05",
       x = "", y = "Number of genes") +
  theme_minimal(base_size = 12)

ggsave(file.path(outdir, "SASP_sig_gene_counts.pdf"),
       p_bar, width = 8, height = 6, dpi = 300)

# ==============================================================================
# 2. Singscore Analysis
# ==============================================================================

# Prepare gene set for singscore (up-regulated SASP signature)
sasp_up_genes <- sasp_found  # Use all SASP genes as the signature

# Get protein-coding genes as background
pc_genes <- keys(org.Hs.eg.db, keytype = "SYMBOL")
gene_info <- AnnotationDbi::select(org.Hs.eg.db, keys = pc_genes,
                                    columns = "GENETYPE", keytype = "SYMBOL")
pc_genes <- gene_info$SYMBOL[gene_info$GENETYPE == "protein-coding"]
cat("Protein-coding genes:", length(pc_genes), "\n")

# Build rank data from log2TPM using protein-coding genes as background
# singscore expects genes as rows, samples as columns
log2tpm_pc <- log2tpm[rownames(log2tpm) %in% pc_genes, , drop = FALSE]
log2tpm_num <- apply(log2tpm_pc, 2, as.numeric)
rownames(log2tpm_num) <- rownames(log2tpm_pc)
rank_data <- singscore::rankGenes(log2tpm_num)

cat("Rank data dimensions (protein-coding):", dim(rank_data), "\n")
cat("Rank data range:", range(rank_data, na.rm = TRUE), "\n")

# Calculate singscore using SASP genes as the upSet against protein-coding background
scoredf <- singscore::simpleScore(rankData = rank_data,
                                   upSet = sasp_up_genes,
                                   centerScore = TRUE)

# Merge with sample design
scoredf$sample <- rownames(scoredf)
scoredf <- scoredf %>%
  left_join(sample_design, by = "sample")

cat("\nSingscore summary:\n")
print(summary(scoredf$TotalScore))

# ==============================================================================
# 3. Singscore Boxplot: NIR vs IR groups
# ==============================================================================

# 3a. NIR vs all IR combined
p_box_all <- ggplot(scoredf, aes(x = ir_group, y = TotalScore, fill = ir_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  scale_fill_manual(values = c("NIR" = "#4DAF4A", "IR" = "#E41A1C")) +
  labs(title = "SASP Singscore: NIR vs All IR Groups",
       x = "", y = "Singscore (TotalScore)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave(file.path(outdir, "SASP_singscore_NIR_vs_IR.pdf"),
       p_box_all, width = 6, height = 6, dpi = 300)

# 3b. NIR vs each IR condition separately
p_box_cond <- ggplot(scoredf, aes(x = condition, y = TotalScore, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  scale_fill_manual(values = c("NIR" = "#4DAF4A", "2Gy24h" = "#FFD92F",
                                "2Gy6d" = "#FC8D62", "10Gy24h" = "#E41A1C",
                                "10Gy6d" = "#984EA3")) +
  labs(title = "SASP Singscore by Condition",
       x = "", y = "Singscore (TotalScore)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(outdir, "SASP_singscore_by_condition.pdf"),
       p_box_cond, width = 8, height = 6, dpi = 300)

# ==============================================================================
# 4. Statistical Tests: Singscore NIR vs IR groups
# ==============================================================================

cat("\n========== Statistical Tests ==========\n")

# 4a. NIR vs all IR combined (Wilcoxon rank-sum test)
nir_scores <- scoredf$TotalScore[scoredf$ir_group == "NIR"]
ir_scores <- scoredf$TotalScore[scoredf$ir_group == "IR"]
ttest_all <- t.test(ir_scores, nir_scores, alternative = "two.sided")
cat(sprintf("\nNIR vs All IR (t-test): p = %.4e\n", ttest_all$p.value))

# 4b. NIR vs each IR condition (t-test)
test_results <- data.frame(
  contrast = character(),
  p_value = numeric(),
  mean_NIR = numeric(),
  mean_IR = numeric(),
  stringsAsFactors = FALSE
)

ir_conditions <- c("2Gy24h", "2Gy6d", "10Gy24h", "10Gy6d")
for (cond in ir_conditions) {
  ir_cond_scores <- scoredf$TotalScore[scoredf$condition == cond]
  tt <- t.test(ir_cond_scores, nir_scores, alternative = "two.sided")
  test_results <- rbind(test_results, data.frame(
    contrast = paste0(cond, "_vs_NIR"),
    p_value = tt$p.value,
    mean_NIR = mean(nir_scores),
    mean_IR = mean(ir_cond_scores),
    stringsAsFactors = FALSE
  ))
  cat(sprintf("%s vs NIR (t-test): p = %.4e, mean NIR=%.3f, mean %s=%.3f\n",
              cond, tt$p.value, mean(nir_scores), cond, mean(ir_cond_scores)))
}

# 4c. Kruskal-Wallis test across all 5 conditions
kw <- kruskal.test(TotalScore ~ condition, data = scoredf)
cat(sprintf("\nKruskal-Wallis across all conditions: p = %.4e\n", kw$p.value))

# Save test results
write.csv(test_results, file.path(outdir, "SASP_singscore_tests.csv"), row.names = FALSE)

# 4d. Add p-value annotations above each IR group (t-test vs NIR)
ir_conditions <- c("2Gy24h", "2Gy6d", "10Gy24h", "10Gy6d")
y_max <- max(scoredf$TotalScore, na.rm = TRUE)
y_range <- diff(range(scoredf$TotalScore, na.rm = TRUE))

pval_labels <- sapply(ir_conditions, function(cond) {
  tt <- t.test(scoredf$TotalScore[scoredf$condition == cond],
                scoredf$TotalScore[scoredf$condition == "NIR"],
                alternative = "two.sided")
  ifelse(tt$p.value < 0.001, sprintf("p=%.2e", tt$p.value), sprintf("p=%.3f", tt$p.value))
})

p_box_cond_stats <- ggplot(scoredf, aes(x = condition, y = TotalScore, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  scale_fill_manual(values = c("NIR" = "#4DAF4A", "2Gy24h" = "#FFD92F",
                                "2Gy6d" = "#FC8D62", "10Gy24h" = "#E41A1C",
                                "10Gy6d" = "#984EA3")) +
  annotate("text", x = 2:5, y = y_max + y_range * 0.06,
           label = pval_labels, size = 3.2, color = "grey30") +
  labs(title = "SASP Singscore by Condition",
       subtitle = sprintf("Kruskal-Wallis p = %.2e  |  p-values: t-test vs NIR", kw$p.value),
       x = "", y = "Singscore (TotalScore)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(outdir, "SASP_singscore_by_condition_stats.pdf"),
       p_box_cond_stats, width = 8, height = 6, dpi = 300)

# ==============================================================================
# 5. Heatmap of log2TPM for SASP genes
# ==============================================================================

# Sort samples by condition
sample_order <- rownames(sample_design)[order(sample_design$condition)]
sasp_log2tpm_ordered <- sasp_log2tpm[, sample_order, drop = FALSE]

# Remove genes with near-zero variance (causes NaN in correlation distance)
gene_vars <- apply(sasp_log2tpm_ordered, 1, var, na.rm = TRUE)
sasp_log2tpm_ordered <- sasp_log2tpm_ordered[gene_vars > 1e-6, , drop = FALSE]
cat("SASP genes retained for heatmap after variance filter:", nrow(sasp_log2tpm_ordered), "\n")

# Annotation for heatmap
annotation_col <- data.frame(
  Condition = sample_design[sample_order, "condition"],
  IR_Group = sample_design[sample_order, "ir_group"],
  row.names = sample_order
)

ann_colors <- list(
  Condition = c("NIR" = "#4DAF4A", "2Gy24h" = "#FFD92F",
                "2Gy6d" = "#FC8D62", "10Gy24h" = "#E41A1C",
                "10Gy6d" = "#984EA3"),
  IR_Group = c("NIR" = "#4DAF4A", "IR" = "#E41A1C")
)

# Heatmap with row scaling (z-score)
pheatmap(sasp_log2tpm_ordered,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         border_color = NA,
         fontsize_row = 8,
         fontsize_col = 9,
         main = "SASP Genes: log2(TPM+1) Z-score",
         filename = file.path(outdir, "SASP_log2TPM_heatmap.pdf"),
         width = 10, height = max(8, length(sasp_found) * 0.28))

# Heatmap without row scaling (raw log2TPM)
pheatmap(sasp_log2tpm_ordered,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         scale = "none",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         border_color = NA,
         fontsize_row = 8,
         fontsize_col = 9,
         main = "SASP Genes: log2(TPM+1)",
         filename = file.path(outdir, "SASP_log2TPM_heatmap_raw.pdf"),
         width = 10, height = max(8, length(sasp_found) * 0.28))

cat("\nHeatmaps saved.\n")

# ==============================================================================
# 6. Export Results
# ==============================================================================

# Export SASP DESeq2 results
write.csv(deseq2_sasp, file.path(outdir, "SASP_DESeq2_results_all_contrasts.csv"), row.names = FALSE)

# Export singscore data
write.csv(scoredf, file.path(outdir, "SASP_singscore_data.csv"), row.names = FALSE)

# Export SASP TPM matrix
write.csv(as.data.frame(sasp_tpm), file.path(outdir, "SASP_TPM_matrix.csv"))

# Create Excel workbook with all results
wb <- createWorkbook()
addWorksheet(wb, "DESeq2_SASP")
writeData(wb, "DESeq2_SASP", deseq2_sasp)
addWorksheet(wb, "Singscore")
writeData(wb, "Singscore", scoredf)
addWorksheet(wb, "Statistical_Tests")
writeData(wb, "Statistical_Tests", test_results)
addWorksheet(wb, "SASP_TPM")
writeData(wb, "SASP_TPM", as.data.frame(sasp_tpm) %>% rownames_to_column("gene"))
saveWorkbook(wb, file.path(outdir, "SASP_analysis_results.xlsx"), overwrite = TRUE)

cat("\n========== SASP Analysis Complete ==========\n")
cat("Results saved to:", outdir, "\n")
cat("Files generated:\n")
for (f in list.files(outdir)) {
  cat("  ", f, "\n")
}
