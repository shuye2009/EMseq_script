# Script to check for ERV in total RNAseq data

library(dplyr)
library(ggplot2)
library(purrr)
library(readr)
library(tibble)
library(tidyr)
library(ComplexHeatmap)
library(pheatmap)
library(VennDiagram)
library(RColorBrewer)
library(edgeR)

enrichment_repeats <- function(de_table, repeat_class){
  # Extract repeat family from row names
  extract_class <- function(x) {
    parts <- strsplit(x, "\\|")
    sapply(parts, function(p) ifelse(length(p) >= 2, p[2], NA))
  }
  
  # Add repeat family column if not already present
  if (!("repeat_class" %in% colnames(de_table))) {
    de_table$repeat_class <- extract_class(rownames(de_table))
  }
  
  # Filter out rows with NA repeat family
  de_table_filtered <- de_table[!is.na(de_table$repeat_class), ]
  
  # Total number of repeats with family information
  total_repeats <- nrow(de_table_filtered)
  
  # Count of the specific repeat family in the entire dataset
  class_count <- sum(de_table_filtered$repeat_class == repeat_class, na.rm = TRUE)
  
  # Count of up-regulated repeats
  up_regulated <- sum(de_table_filtered$is_significant == "Up", na.rm = TRUE)
  
  # Count of down-regulated repeats
  down_regulated <- sum(de_table_filtered$is_significant == "Down", na.rm = TRUE)
  
  # Count of the specific repeat family in up-regulated repeats
  class_up <- sum(de_table_filtered$repeat_class == repeat_class & 
                    de_table_filtered$is_significant == "Up", na.rm = TRUE)
  
  # Count of the specific repeat family in down-regulated repeats
  class_down <- sum(de_table_filtered$repeat_class == repeat_class & 
                      de_table_filtered$is_significant == "Down", na.rm = TRUE)
  
  # Create contingency tables for Fisher's exact test
  # For up-regulation
  up_table <- matrix(c(
    class_up,                                # Family & Up-regulated
    class_count - class_up,                 # Family & Not up-regulated
    up_regulated - class_up,                 # Not family & Up-regulated
    total_repeats - class_count - up_regulated + class_up  # Not family & Not up-regulated
  ), nrow = 2, byrow = TRUE)
  
  # For down-regulation
  down_table <- matrix(c(
    class_down,                                  # Family & Down-regulated
    class_count - class_down,                   # Family & Not down-regulated
    down_regulated - class_down,                 # Not family & Down-regulated
    total_repeats - class_count - down_regulated + class_down  # Not family & Not down-regulated
  ), nrow = 2, byrow = TRUE)
  
  # Perform Fisher's exact test
  fisher_up <- fisher.test(up_table, alternative = "greater")
  fisher_down <- fisher.test(down_table, alternative = "greater")
  
  # Calculate fold enrichment
  fold_enrichment_up <- ifelse(up_regulated > 0, 
                              (class_up / up_regulated) / (class_count / total_repeats), 
                              NA)
  
  fold_enrichment_down <- ifelse(down_regulated > 0, 
                                (class_down / down_regulated) / (class_count / total_repeats), 
                                NA)
  
  # Return results as a data frame
  result <- data.frame(
    repeat_class = repeat_class,
    total_repeats = total_repeats,
    class_count = class_count,
    up_regulated = up_regulated,
    down_regulated = down_regulated,
    class_up = class_up,
    class_down = class_down,
    p_value_up = fisher_up$p.value,
    p_value_down = fisher_down$p.value,
    odds_ratio_up = fisher_up$estimate,
    odds_ratio_down = fisher_down$estimate,
    fold_enrichment_up = fold_enrichment_up,
    fold_enrichment_down = fold_enrichment_down
  )
  
  return(result)
}

setwd("C:/PROJECTS/Shane/totalRNAseq")
outDir <- "C:/PROJECTS/Shane/totalRNAseq/plots"
if(!dir.exists(outDir)) dir.create(outDir)

# Read in the data
r_counts <- as.data.frame(readr::read_tsv("rCounts.txt", col_names = TRUE, skip = 1)) 

# Clean up the gene ids and remove columns Chr, Start, End, Strand, remove finalbam/ and .sort.bam from column names
r_counts <- r_counts %>%
  dplyr::mutate(Geneid = gsub("\\#|\\/", "|", Geneid)) %>%
  dplyr::mutate(Geneid = ifelse(grepl("Simple_repeat", Geneid, fixed = FALSE), paste0(Geneid, "|Simple_repeat"), Geneid)) %>%
  dplyr::rename_with(~gsub("finalbam/|\\.sort\\.bam", "", .), -c(Geneid, Length)) %>%
  dplyr::select(-c(Chr, Start, End, Strand))
head(r_counts)
summary(r_counts)

# Select columns starting with A and convert to matrix and set row names to Geneid
r_counts_matrix <- r_counts %>%
  dplyr::select(starts_with("A")) %>%
  as.matrix()
# replace NA with 0
r_counts_matrix[is.na(r_counts_matrix)] <- 0
rownames(r_counts_matrix) <- r_counts$Geneid
summary(r_counts_matrix)

# Scale the data for better visualization
r_counts_scaled <- log(r_counts_matrix + 1)
# Read in metadata and create a data frame as column annotations
metadata <- as.data.frame(readr::read_tsv("SH07_metadata.txt", col_names = TRUE, skip = 0))
metadata <- metadata %>%
  dplyr::mutate(sampleID = gsub("finalbam\\.|\\.sort\\.bam", "", sampleID)) %>%
  dplyr::select(-replicate)
head(metadata)  

# Create column annotation data frame - ensure row names match column names in matrix
column_anno <- metadata %>%
  dplyr::filter(sampleID %in% colnames(r_counts_matrix)) %>%
  dplyr::arrange(match(sampleID, colnames(r_counts_matrix))) %>%
  tibble::column_to_rownames(var = "sampleID")


# Use a colorblind-friendly palette with distinct colors from multiple palettes
# Combine multiple palettes to get at least 21 distinct colors
set1_colors <- brewer.pal(9, "Set1")
set2_colors <- brewer.pal(8, "Set2")
set3_colors <- brewer.pal(12, "Set3")
distinct_colors <- c(set1_colors, set2_colors, set3_colors)
# Ensure we have at least 21 unique colors
distinct_colors <- distinct_colors[1:21]
# Create a named vector for condition colors
condition_colors <- setNames(distinct_colors[1:length(unique(column_anno$condition))], 
                            unique(column_anno$condition))


# Create the heatmap
pdf(file.path(outDir, "total_rnaseq_readCounts_in_repeats_heatmap.pdf"), width = 8, height = 12)
p <- pheatmap(
  r_counts_scaled,
  annotation_col = column_anno,
  annotation_colors = list(condition = condition_colors),
  show_rownames = FALSE,  # Too many rows to display names
  fontsize_col = 8,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Total RNAseq Read Counts Heatmap",
  raster = TRUE
)
print(p)
if (!is.null(dev.list())) dev.off()


# Perform edgeR analysis to contrast treatment(such aspIR.3d) vs control(=pNIR) ####

# Create DGEList object
dge <- DGEList(counts = r_counts_matrix)

# Filter out low count genes
keep <- filterByExpr(dge, group = column_anno$condition,
                    min.count = 3, min.prop = 0.2, min.total.count = 6)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalize for library size
dge <- calcNormFactors(dge)

# Create design matrix focusing on the contrast of interest
# Extract samples belonging to treatment and control groups
control <- "pNIR"
treatments <- unique(column_anno$condition[column_anno$condition != control])

# store the results of the enrichment test in a list
enrichment_results <- list()
for(treatment in treatments){

  treatmentDir <- file.path(outDir, treatment)
  if (!dir.exists(treatmentDir)) dir.create(treatmentDir)
  
  pir_pnir_samples <- column_anno %>%
    filter(condition %in% c(treatment, control)) %>%
    rownames()

  # Subset DGE object to include only these samples
  dge_contrast <- dge[, pir_pnir_samples]

  # Create a factor for the groups
  group <- factor(column_anno[pir_pnir_samples, "condition"])

  # Create design matrix
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)

  # Estimate dispersion
  dge_contrast <- estimateDisp(dge_contrast, design)

  # Fit the model
  fit <- glmQLFit(dge_contrast, design)

  # Define the contrast (treatment vs control)
  contrast_formula <- paste0("`", treatment, "`-`", control, "`")
  my_contrast <- makeContrasts(contrasts=contrast_formula, levels=design)

  # Perform the test
  qlf <- glmQLFTest(fit, contrast = my_contrast)

  # Get the results
  de_results <- topTags(qlf, n = Inf)
  de_table <- as.data.frame(de_results)

  # Add a column to indicate if the gene is significantly DE
  de_table$is_significant <- ifelse(de_table$FDR < 0.05,
                                    ifelse(de_table$logFC > 1, "Up", 
                                    ifelse(de_table$logFC < -1, "Down", "Not Significant")), 
                                    "Not Significant") 

  # Write results to file
  write.csv(de_table, file = file.path(treatmentDir, paste0("edgeR_", treatment, "_vs_", control, ".csv")))

  # Test enrichment of Simple_repeat
  Simple_repeat_enrichment <- enrichment_repeats(de_table, "Simple_repeat")
  enrichment_results[[treatment]] <- Simple_repeat_enrichment

  # Create MA plot
  pdf(file.path(treatmentDir, paste0("MA_plot_", treatment, "_vs_", control, ".pdf")))
  with(de_table, plot(logCPM, logFC, 
                    pch = 20, 
                    col = ifelse(is_significant == "Up", "red", 
                                ifelse(is_significant == "Down", "blue", "gray")),
                    main = paste0("MA Plot: ", treatment, " vs ", control),
                    xlab = "Average logCPM",
                    ylab = "logFC"))
  abline(h = 0, col = "black", lty = 2)
  legend("topright", 
        legend = c("Up", "Down", "Not Significant"), 
        col = c("red", "blue", "gray"), 
        pch = 20)
  dev.off()

  de_genes <- rownames(de_table[de_table$is_significant != "Not Significant", ])
  erv_de_genes <- de_genes[grepl("ERV", de_genes)]
  erv_up_genes <- erv_de_genes[de_table[erv_de_genes, "is_significant"] == "Up"]
  erv_down_genes <- erv_de_genes[de_table[erv_de_genes, "is_significant"] == "Down"]

  # Partition all genes into up, down, not significant, then make a row annotation
  gene_anno <- data.frame(
    Repeat = rownames(de_table),
    DE_status = de_table$is_significant,
    stringsAsFactors = TRUE
  )
  rownames(gene_anno) <- gene_anno$Repeat

  # Create volcano plot
  pdf(file.path(treatmentDir, paste0("Volcano_plot_", treatment, "_vs_", control, ".pdf")), width = 12, height = 10)
  # Basic plot
  with(de_table, plot(logFC, -log10(PValue), 
                    pch = 20, 
                    col = ifelse(is_significant == "Up", "red", 
                                ifelse(is_significant == "Down", "blue", "gray")),
                    main = paste0("Volcano Plot: ", treatment, " vs ", control),
                    xlab = "logFC",
                    ylab = "-log10(P-value)"))
  abline(v = -1, col = "black", lty = 2)
  abline(v = 1, col = "black", lty = 2)
  abline(h = -log10(0.05), col = "black", lty = 2)
  legend("topright", 
        legend = c("Up", "Down", "Not Significant"), 
        col = c("red", "blue", "gray"), 
        pch = 20)

  # Add labels for top DE ERV genes
  if(length(erv_de_genes) > 0) {
    # Get coordinates for the top DE ERV genes
    erv_de_coords <- de_table[erv_de_genes, c("logFC", "PValue")]
    
    # Add text labels with gene names
    # Use a small subset if there are too many genes to label clearly
    max_labels <- min(50, length(erv_de_genes))
    # Sort by significance to label the most significant ones
    top_labels <- erv_de_genes[order(de_table[erv_de_genes, "PValue"])[1:max_labels]]
    
    # Add labels with slight offset to avoid overlapping with points
    text(de_table[top_labels, "logFC"], 
        -log10(de_table[top_labels, "PValue"]), 
        labels = top_labels, 
        cex = 0.7,  # Smaller text size
        pos = 3,    # Position above the point
        offset = 0.5)
  }

  dev.off()

  # Create a heatmap of top DE genes

  if (length(de_genes) > 0) {
    # Extract counts for top DE genes
    top_counts <- r_counts_scaled[de_genes, pir_pnir_samples]
    de_anno <- gene_anno %>%
      dplyr::select(-Repeat)
    # Create heatmap
    pdf(file.path(treatmentDir, paste0("DE_genes_heatmap_", treatment, "_vs_", control, ".pdf")), width = 10, height = 10)
    pheatmap(
      top_counts,
      annotation_row = de_anno,
      annotation_col = column_anno,
      annotation_colors = list(
        DE_status = c("Up" = "red", "Down" = "blue", "Not Significant" = "gray"),
        condition = condition_colors
      ),
      scale = "row",  # Scale by row to show relative expression
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_rownames = ifelse(length(de_genes) <= 100, TRUE, FALSE),
      fontsize_row = 8,
      fontsize_col = 8,
      main = paste0("DE Genes: ", treatment, " vs ", control)
    )
    dev.off()
  }

  # Create a heatmap of top DE ERV genes

  if (length(erv_de_genes) > 0) {
    # Extract counts for top DE ERV genes
    top_counts <- r_counts_scaled[erv_de_genes, pir_pnir_samples]
    evr_anno <- gene_anno %>%
      dplyr::select(-Repeat)
    # Create heatmap
    pdf(file.path(treatmentDir, paste0("ERV_DE_genes_heatmap_", treatment, "_vs_", control, ".pdf")), width = 10, height = 10)
    pheatmap(
      top_counts,
      annotation_row = evr_anno,
      annotation_col = column_anno,
      annotation_colors = list(
        DE_status = c("Up" = "red", "Down" = "blue", "Not Significant" = "gray"),
        condition = condition_colors
      ),
      scale = "row",  # Scale by row to show relative expression
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_rownames = ifelse(length(erv_de_genes) <= 100, TRUE, FALSE),
      fontsize_row = 8,
      fontsize_col = 8,
      main = paste0("ERV DE Genes: ", treatment, " vs ", control)
    )
    dev.off()
  }


  # Get hyper and hypo methylated ERVs in three contexts
  for(context in c("CG")){
    hyper_anno_file <- paste0("C:\\PROJECTS\\Shane\\Harding_combined\\h4h\\DSS\\twoGroupComp\\", context, "\\WT-IR_vs_WT-NIR\\HOMER\\hyper\\DMRs_hyper_annot.tab")
    hypo_anno_file <- paste0("C:\\PROJECTS\\Shane\\Harding_combined\\h4h\\DSS\\twoGroupComp\\", context, "\\WT-IR_vs_WT-NIR\\HOMER\\hypo\\DMRs_hypo_annot.tab")

    hyper_file <- paste0("C:\\PROJECTS\\Shane\\Harding_combined\\h4h\\DSS\\twoGroupComp\\", context, "\\WT-IR_vs_WT-NIR\\DMRs\\DMR_hyper.bed")
    hypo_file <- paste0("C:\\PROJECTS\\Shane\\Harding_combined\\h4h\\DSS\\twoGroupComp\\", context, "\\WT-IR_vs_WT-NIR\\DMRs\\DMR_hypo.bed")
    
    contextDir <- file.path(treatmentDir, context)
    if(!dir.exists(contextDir)) dir.create(contextDir)
    # Read the file with parameters to handle inconsistent columns
    hyper_anno <- read.table(hyper_anno_file, 
                      header = TRUE, 
                      sep = "\t", 
                      fill = TRUE,       # Fill columns with NA if missing
                      quote = "",        # No quote characters
                      check.names = FALSE, # Don't modify column names
                      comment.char = "")  # Don't treat any character as comments
    # Convert column names to lower case
    colnames(hyper_anno) <- tolower(colnames(hyper_anno))

    dim(hyper_anno)

    hypo_anno <- read.table(hypo_anno_file, 
                      header = TRUE, 
                      sep = "\t", 
                      fill = TRUE,       # Fill columns with NA if missing
                      quote = "",        # No quote characters
                      check.names = FALSE, # Don't modify column names
                      comment.char = "")  # Don't treat any character as comments
    # Convert column names to lower case
    colnames(hypo_anno) <- tolower(colnames(hypo_anno))

    dim(hypo_anno)

    hyper <- read.table(hyper_file, 
                      header = TRUE, 
                      sep = "\t", 
                      fill = TRUE,       # Fill columns with NA if missing
                      quote = "",        # No quote characters
                      check.names = FALSE, # Don't modify column names
                      comment.char = "")  %>%
                      dplyr::select(-strand)

    dim(hyper)

    hypo <- read.table(hypo_file, 
                      header = TRUE, 
                      sep = "\t", 
                      fill = TRUE,       # Fill columns with NA if missing
                      quote = "",        # No quote characters
                      check.names = FALSE, # Don't modify column names
                      comment.char = "")  %>%
                      dplyr::select(-strand)  

    dim(hypo)

    # Extract the repeat names, class, and family by splitting on the pipe character
    hyper_repeat_info <- hyper_anno %>%
      dplyr::filter(grepl("\\|", `detailed annotation`)) %>%
      dplyr::mutate(
        start = start - 1, # adjust start position for 0-based indexing
        repeat_full = `detailed annotation`,
        repeat_parts = strsplit(`detailed annotation`, "\\|"),
        repeat_name = sapply(repeat_parts, function(x) ifelse(length(x) >= 1, x[1], NA)),
        repeat_class = sapply(repeat_parts, function(x) ifelse(length(x) >= 2, x[2], NA)),
        repeat_family = sapply(repeat_parts, function(x) ifelse(length(x) >= 3, x[3], NA)),
        DE_status = de_table[repeat_full, "is_significant"]
      ) %>%
      dplyr::select(-repeat_parts)

    hypo_repeat_info <- hypo_anno %>%
      dplyr::filter(grepl("\\|", `detailed annotation`)) %>%
      dplyr::mutate(
        start = start - 1, # adjust start position for 0-based indexing
        repeat_full = `detailed annotation`,
        repeat_parts = strsplit(`detailed annotation`, "\\|"),
        repeat_name = sapply(repeat_parts, function(x) ifelse(length(x) >= 1, x[1], NA)),
        repeat_class = sapply(repeat_parts, function(x) ifelse(length(x) >= 2, x[2], NA)),
        repeat_family = sapply(repeat_parts, function(x) ifelse(length(x) >= 3, x[3], NA)),
        DE_status = de_table[repeat_full, "is_significant"]
      ) %>%
      dplyr::select(-repeat_parts)

    # Join hyper_repleat_info with hyper
    hyper_repeat_info <- hyper_repeat_info %>%
      dplyr::left_join(hyper, by = c("chr", "start", "end"))

    # Filter for significant repeats
    hyper_repeat_DE <- hyper_repeat_info %>%
      #dplyr::filter(DE_status == "Up" | DE_status == "Down") %>%
      dplyr::select(c(chr, start, end, strand, repeat_name, repeat_class, repeat_family, stat, length, nCG, ratio, DE_status)) %>%
      dplyr::arrange(desc(abs(stat)))
    
    # Join hypo_repleat_info with hypo
    hypo_repeat_info <- hypo_repeat_info %>%
      dplyr::left_join(hypo, by = c("chr", "start", "end"))

    # Filter for significant repeats
    hypo_repeat_DE <- hypo_repeat_info %>%
      #dplyr::filter(DE_status == "Up" | DE_status == "Down") %>%
      dplyr::select(c(chr, start, end, strand, repeat_name, repeat_class, repeat_family, stat, length, nCG, ratio, DE_status)) %>%
      dplyr::arrange(desc(abs(stat)))
    
    # Output the significant repeats
    write.table(hyper_repeat_DE, file.path(contextDir, "hyper_significantDE_repeats.tab"), sep = "\t", row.names = FALSE)
    write.table(hypo_repeat_DE, file.path(contextDir, "hypo_significantDE_repeats.tab"), sep = "\t", row.names = FALSE)

    
    # Count occurrences of each repeat class
    hyper_class_counts <- table(hyper_repeat_info$repeat_class)
    hypo_class_counts <- table(hypo_repeat_info$repeat_class)

    # Count occurrences of each repeat family
    hyper_class_counts <- table(hyper_repeat_info$repeat_class)
    hypo_class_counts <- table(hypo_repeat_info$repeat_class)

    # Display the top repeat classes
    print("Top repeat classes in hyper DMRs:")
    head(sort(hyper_class_counts, decreasing = TRUE), 10)

    print("Top repeat classes in hypo DMRs:")
    head(sort(hypo_class_counts, decreasing = TRUE), 10)

    # Display the top repeat families
    print("Top repeat families in hyper DMRs:")
    head(sort(hyper_class_counts, decreasing = TRUE), 10)

    print("Top repeat families in hypo DMRs:")
    head(sort(hypo_class_counts, decreasing = TRUE), 10)
    hyper_repeats <- hyper_repeat_info %>% 
      dplyr::select(repeat_full)

    hypo_repeats <- hypo_repeat_info %>% 
      dplyr::select(repeat_full)

    # Count the number of repeats found
    print(paste("Number of repeats in hyper DMRs:", nrow(hyper_repeats)))
    print(paste("Number of repeats in hypo DMRs:", nrow(hypo_repeats)))

    # Get unique repeats
    hyper_repeats_unique <- unique(hyper_repeats$repeat_full)
    hypo_repeats_unique <- unique(hypo_repeats$repeat_full) 
    repeats <- unique(c(hyper_repeats_unique, hypo_repeats_unique)) 

    # Create a Venn diagram of the unique repeats
    # Set up the directory for the Venn diagram output


    # Create the Venn diagram
    venn.plot <- VennDiagram::venn.diagram(
      x = list(
        "Hyper DM repeats" = hyper_repeats_unique,
        "Hypo DM repeats" = hypo_repeats_unique
      ),
      filename = file.path(contextDir, "repeat_overlap.png"),
      output = TRUE,
      
      # Customize the appearance
      col = c("red", "blue"),
      fill = c(alpha("red", 0.3), alpha("blue", 0.3)),
      alpha = 0.5,
      
      # Add labels
      cat.col = c("red", "blue"),
      cat.cex = 1.5,
      cat.pos = c(0, 0),
      cat.dist = c(0.05, 0.05),
      
      # Set the title
      main = "Overlap of Repeats",
      main.cex = 2
    )

    # Print the counts for the Venn diagram
    print(paste("Number of unique repeats in hyper DMRs:", length(hyper_repeats_unique)))
    print(paste("Number of unique repeats in hypo DMRs:", length(hypo_repeats_unique)))
    print(paste("Number of repeats in both hyper and hypo DMRs:", 
                length(intersect(hyper_repeats_unique, hypo_repeats_unique))))
    print(paste("Total number of unique repeats:", length(repeats)))

    # Get repeats unique to hyper and hypo
    hyper_unique <- setdiff(hyper_repeats_unique, hypo_repeats_unique)
    hypo_unique <- setdiff(hypo_repeats_unique, hyper_repeats_unique)
    # Get repeats shared by hyper and hypo
    shared_repeats <- intersect(hyper_repeats_unique, hypo_repeats_unique)
    # Create a row annotation based on the unique and shared repeats
    row_anno <- data.frame(
      Repeat = c(hyper_unique, hypo_unique, shared_repeats),
      Type = c(rep("Hyper Unique", length(hyper_unique)), 
      rep("Hypo Unique", length(hypo_unique)), 
      rep("Shared", length(shared_repeats))),
      stringsAsFactors = FALSE
    ) 

    row_anno <- left_join(row_anno, gene_anno, by = "Repeat") %>%
    tibble::column_to_rownames(var = "Repeat") 

    # Create a heatmap of the repeat counts
    tmp <- r_counts_scaled[rownames(r_counts_scaled) %in% rownames(row_anno), pir_pnir_samples]
    summary(tmp)
    # Remove rows with all zeros
    tmp <- tmp[rowSums(tmp) != 0, ]
    # Remove rows with identical values across all samples
    # Calculate row variance - rows with identical values will have variance of 0
    row_var <- apply(tmp, 1, var)
    # Keep only rows with non-zero variance
    tmp <- tmp[row_var > 0, ]
    # Update row_anno to match the filtered tmp matrix
    row_anno <- row_anno[rownames(tmp), , drop = FALSE]
    
    pdf(file.path(contextDir, "Read_counts_in_repeat_heatmap.pdf"), width = 8, height = 16)
    p <- pheatmap(
      tmp,
      annotation_row = row_anno,
      annotation_col = column_anno,
      annotation_colors = list(
        Type = c("Hyper Unique" = "#E41A1C", "Hypo Unique" = "#377EB8", "Shared" = "#4DAF4A"),
        DE_status = c("Up" = "magenta", "Down" = "cyan", "Not Significant" = "gray"),
        condition = condition_colors
      ),
      scale = "row",  # Scale by row to show relative expression
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      fontsize_row = 6,
      fontsize_col = 8,
      main = "Read Counts in DM Repeats",
      border = TRUE
    )
    print(p)
    dev.off()

    # Look for ERV-specific repeats (filter for classes containing "ERV" or "LTR")
    hyper_erv <- hyper_repeat_info %>%
      dplyr::filter(grepl("ERV", repeat_family, ignore.case = TRUE))

    hypo_erv <- hypo_repeat_info %>%
      dplyr::filter(grepl("ERV", repeat_family, ignore.case = TRUE))

    # Count ERV/LTR classes
    hyper_erv_family_counts <- table(hyper_erv$repeat_family)
    hypo_erv_family_counts <- table(hypo_erv$repeat_family)

    print("ERV/LTR classes in hyper DMRs:")
    print(hyper_erv_family_counts)

    print("ERV/LTR classes in hypo DMRs:")
    print(hypo_erv_family_counts)

    # Get unique repeats
    hyper_erv_unique <- unique(hyper_erv$repeat_full)
    hypo_erv_unique <- unique(hypo_erv$repeat_full)

    # Create a Venn diagram of the unique ERV-specific repeats
    venn.plot.erv <- VennDiagram::venn.diagram(
      x = list(
        "Hyper DM ERV" = hyper_erv_unique,
        "Hypo DM ERV" = hypo_erv_unique
      ),
      filename = file.path(contextDir, "erv_repeat_overlap.png"),
      output = TRUE,
      
      # Customize the appearance
      col = c("red", "blue"),
      fill = c(alpha("red", 0.3), alpha("blue", 0.3)),
      alpha = 0.5,
      
      # Add labels
      cat.col = c("red", "blue"),
      cat.cex = 1.5,
      cat.pos = c(0, 0),
      cat.dist = c(0.05, 0.05),
      
      # Set the title
      main = "Overlap of DM ERV",
      main.cex = 2
    )

    # Print the counts for the Venn diagram
    print(paste("Number of unique DM ERV in hyper DMRs:", length(hyper_erv_unique)))
    print(paste("Number of unique DM ERV in hypo DMRs:", length(hypo_erv_unique)))
    print(paste("Number of DM ERV in both hyper and hypo DMRs:", 
                length(intersect(hyper_erv_unique, hypo_erv_unique))))
    print(paste("Total number of unique DM ERV:", length(unique(c(hyper_erv_unique, hypo_erv_unique)))))

    # Get ERV row annotation
    row_anno_erv <- row_anno[unique(c(hyper_erv_unique, hypo_erv_unique)), , drop = FALSE]

    # Create a heatmap of the repeat counts
    pdf(file.path(contextDir, "Read_counts_in_ERV_heatmap.pdf"))
    p <- pheatmap(
      r_counts_scaled[rownames(r_counts_scaled) %in% rownames(row_anno_erv), pir_pnir_samples],
      annotation_row = row_anno_erv,
      annotation_col = column_anno,
      annotation_colors = list(
        Type = c("Hyper Unique" = "#E41A1C", "Hypo Unique" = "#377EB8", "Shared" = "#4DAF4A"),
        DE_status = c("Up" = "magenta", "Down" = "cyan", "Not Significant" = "gray"),
        condition = condition_colors
      ),
      scale = "row",
      fontsize_row = 6,
      fontsize_col = 8,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      main = "Read Counts in DM ERV",
      border = TRUE
    )
    print(p)
    dev.off()

    # remove all ".log" files
    log_files <- list.files(path = contextDir, pattern = "\\.log$", full.names = TRUE)
    if(length(log_files) > 0) {
      file.remove(log_files)
      print(paste("Removed", length(log_files), "log files"))
    }

  }
}


# Convert the enrichment results to a data frame
enrichment_results_df <- do.call(rbind, enrichment_results)
write.csv(enrichment_results_df, file = file.path(outDir, "enrichment_results_Simple_repeat.csv"))

# Create a matrix of repeats by treatment showing DE status
# First, collect all unique repeats across all treatments
all_repeats <- list()
de_tables <- list()

# Rerun the loop to collect all de_tables
for(treatment in treatments) {
  treatmentDir <- file.path(outDir, treatment)
  de_file <- file.path(treatmentDir, paste0("edgeR_", treatment, "_vs_", control, ".csv"))
  
  # Check if the file exists
  if(file.exists(de_file)) {
    # Read the DE table
    de_table <- read.csv(de_file, row.names = 1)
    de_tables[[treatment]] <- de_table
    
    # Collect all repeats
    all_repeats[[treatment]] <- rownames(de_table)
  }
}

# Get unique repeats across all treatments
unique_repeats <- unique(unlist(all_repeats))

# Create an empty matrix to store the DE status
de_matrix <- matrix(0, nrow = length(unique_repeats), ncol = length(treatments))
rownames(de_matrix) <- unique_repeats
colnames(de_matrix) <- treatments

# Fill the matrix with DE status: 1 for up, -1 for down, 0 for not significant
for(i in seq_along(treatments)) {
  treatment <- treatments[i]
  
  if(treatment %in% names(de_tables)) {
    de_table <- de_tables[[treatment]]
    
    # Get indices of repeats in the matrix
    repeat_indices <- match(rownames(de_table), unique_repeats)
    
    # Set values based on DE status
    de_matrix[repeat_indices[!is.na(repeat_indices)], i] <- 
      ifelse(de_table$is_significant == "Up", 1,
             ifelse(de_table$is_significant == "Down", -1, 0))
  }
}

# Convert to data frame for easier handling
de_status_df <- as.data.frame(de_matrix)

# Write the matrix to a CSV file
write.csv(de_status_df, file = file.path(outDir, "repeat_DE_status_matrix.csv"))

# Create a heatmap of the DE status matrix
# Extract repeat family for annotation
extract_class <- function(x) {
  parts <- strsplit(x, "\\|")
  sapply(parts, function(p) ifelse(length(p) >= 2, p[2], NA))
}

# Create annotation for rows (repeat families)
row_anno <- data.frame(
  repeat_class = extract_class(rownames(de_status_df))
)
rownames(row_anno) <- rownames(de_status_df)
class_colors <- setNames(distinct_colors[1:length(unique(row_anno$repeat_class))], 
                            unique(row_anno$repeat_class))

# Create a custom color palette for the heatmap
de_colors <- c("-1" = "blue", "0" = "white", "1" = "red")

# Generate the heatmap
pdf(file.path(outDir, "repeat_DE_status_heatmap.pdf"), width = 10, height = 20)
pheatmap(
  de_matrix,
  annotation_row = row_anno,
  annotation_colors = list(
    repeat_class = class_colors
  ),
  color = de_colors,
  breaks = c(-1.5, -0.5, 0.5, 1.5),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,  # Too many rows to display names
  fontsize_col = 10,
  border = TRUE,
  border_color = "grey25",
  main = "Differential Expression Status Across Treatments"
)
dev.off()

# Create a filtered version with only repeats that are DE in at least one treatment
de_in_any <- apply(de_matrix, 1, function(x) any(x != 0))
de_filtered_matrix <- de_matrix[de_in_any, ]

# Create a heatmap of the filtered DE status matrix
pdf(file.path(outDir, "repeat_DE_status_filtered_heatmap.pdf"), width = 10, height = 20)
filtered_row_anno <- row_anno[rownames(de_filtered_matrix), , drop = FALSE]
de_filtered_class_colors <- class_colors[unique(filtered_row_anno$repeat_class)]
pheatmap(
  de_filtered_matrix,
  annotation_row = filtered_row_anno,
  annotation_colors = list(
    repeat_class = de_filtered_class_colors
  ),
  color = de_colors,
  breaks = c(-1.5, -0.5, 0.5, 1.5),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border = TRUE,
  border_color = "green4",
  show_rownames = ifelse(nrow(de_filtered_matrix) <= 100, TRUE, FALSE),
  fontsize_row = 6,
  fontsize_col = 10,
  main = "Differential Expression Status Across Treatments (DE in at least one treatment)"
)
dev.off()

# Create a heatmap of the filtered DE status matrix for pIR.3d and RNASEL.3d
pIR_RNASEL_matrix <- de_filtered_matrix[, c("pIR.3d", "RNASEL.3d")]
#remove all 0 rows
pIR_RNASEL_matrix <- pIR_RNASEL_matrix[rowSums(pIR_RNASEL_matrix != 0) > 0, ]

pdf(file.path(outDir, "pIR_RNASEL_DE_heatmap.pdf"), width = 10, height = 20)

# Choose class colors only relevant to pIR.3d and RNASEL.3d 
pIR_RNASEL_row_anno <- row_anno[rownames(pIR_RNASEL_matrix), , drop = FALSE]
pIR_RNASEL_class_colors <- class_colors[unique(pIR_RNASEL_row_anno$repeat_class)]
pheatmap(
  pIR_RNASEL_matrix,
  annotation_row = pIR_RNASEL_row_anno,
  annotation_colors = list(
    repeat_class = pIR_RNASEL_class_colors
  ),
  color = de_colors,
  breaks = c(-1.5, -0.5, 0.5, 1.5),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border = TRUE,
  border_color = "black",
  show_rownames = ifelse(nrow(pIR_RNASEL_matrix) <= 200, TRUE, FALSE),
  fontsize_row = 6,
  fontsize_col = 10,
  main = "Differential Expression Status Across Treatments (pIR.3d vs RNASEL.3d)"
)
dev.off()

# Create a heatmap of the filtered DE status matrix all treatments, but only for LTR repeats
LTR_anno <- row_anno[row_anno$repeat_class == "LTR", , drop = FALSE]
LTR_de_matrix <- de_filtered_matrix[rownames(de_filtered_matrix) %in% rownames(LTR_anno), ]
#remove all 0 rows
LTR_de_matrix <- LTR_de_matrix[rowSums(LTR_de_matrix != 0) > 0, ]
# write the LTR_de_matrix to a csv file
write.csv(LTR_de_matrix, file = file.path(outDir, "LTR_DE_matrix.csv"))

pdf(file.path(outDir, "LTR_DE_heatmap.pdf"), width = 10, height = 10)
pheatmap(
  LTR_de_matrix,
  annotation_row = LTR_anno,
  annotation_colors = list(
    repeat_class = class_colors
  ),
  color = de_colors,
  breaks = c(-1.5, -0.5, 0.5, 1.5),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border = TRUE,
  border_color = "black",
  show_rownames = ifelse(nrow(LTR_de_matrix) <= 200, TRUE, FALSE),
  fontsize_row = 6,
  fontsize_col = 10,
  main = "Differential Expression Status Across Treatments (LTR)"
)
dev.off()

