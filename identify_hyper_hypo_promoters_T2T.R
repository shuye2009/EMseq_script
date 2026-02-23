library(GenomicRanges)
library(data.table)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)
library(rtracklayer)

# ==============================================================================
# Settings and Paths
# ==============================================================================

# Input Paths
promoter_bed_file <- "C:/PROJECTS/resource/T2T_CHM13/T2T_CHM13_promoters.bed"
dmrseq_base_dir   <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/dmrseq"
targeted_base_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/Promoter_targeted"
rnaseq_base_dir   <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/RNAseq/deseq2/single_factor_analysis_gene_level"

# Output Path
output_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/Promoter_methylation_summary"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Targeted file suffix pattern
targeted_suffix <- "all_diff_T2T_CHM13_promoters.tab"

# Filtering Thresholds
dmrseq_pval_cutoff <- 0.05
dmrseq_beta_cutoff <- 0.25

targeted_pval_cutoff <- 0.05
targeted_diff_cutoff <- 10

rnaseq_padj_cutoff <- 0.05
rnaseq_lfc_cutoff <- 1

# ==============================================================================
# Load Promoter Annotations & Gene Links
# ==============================================================================
message("Loading Promoter Annotations...")
# Expecting standard 6-column BED: chrom, start, end, name, score, strand
promoters <- import(promoter_bed_file, format = "BED")
# Ensure names are unique or usable as IDs
if (any(duplicated(promoters$name))) {
  message("Found duplicate promoter names, making them unique...")
  promoters$name <- paste0(promoters$name, "_", seq_along(promoters$name))
}

message("Loading Promoter-Gene Links...")
# The promoter BED file contains gene info in the 4th column (name)
# Format: GeneName_ENSGID
# We need to extract the gene name and create a mapping
promoter_dt <- data.table(
  name = promoters$name,
  chrom = as.character(seqnames(promoters)),
  start = start(promoters),
  end = end(promoters),
  strand = as.character(strand(promoters))
)

# Extract gene name from the name column (split by "_")
promoter_dt$GeneName <- sub("_.*", "", promoter_dt$name)

# Create a gene map: promoter name -> gene name
gene_map <- unique(promoter_dt[, c("name", "GeneName")])
setkey(gene_map, name)

# ==============================================================================
# Initialize Summary Stats
# ==============================================================================
summary_stats <- list()

# ==============================================================================
# Get Comparison List
# ==============================================================================
message("Finding comparisons...")
comps <- list.dirs(dmrseq_base_dir, recursive = FALSE)
comps <- basename(comps)
# Filter to only include directories with "_vs_" in the name
comps <- comps[grepl("_vs_", comps)]
message("Found ", length(comps), " comparisons: ", paste(comps, collapse = ", "))

# ==============================================================================
# Helper Functions
# ==============================================================================
# Function to annotate regulated genes from RNA-seq data
annotate_rna <- function(target_genes_str, reg_map) {
  if (is.na(target_genes_str) || target_genes_str == "") return(NA_character_)
  genes <- unlist(strsplit(target_genes_str, ", "))
  genes <- trimws(genes)
  
  reg_status <- reg_map[genes]
  hits <- !is.na(reg_status)
  
  if (!any(hits)) return(NA_character_)
  
  return(paste(paste0(genes[hits], "(", reg_status[hits], ")"), collapse = ", "))
}

# Function to count unique regulated genes from a dataframe
count_reg_genes <- function(df) {
  if (is.null(df) || !"Regulated_TargetGenes" %in% colnames(df)) {
    return(c(Up = 0, Down = 0))
  }
  
  regulated_col <- na.omit(df$Regulated_TargetGenes)
  regulated_col <- regulated_col[regulated_col != ""]
  
  if (length(regulated_col) == 0) {
    return(c(Up = 0, Down = 0))
  }
  
  # Extract all gene names with regulation status
  all_entries <- unlist(strsplit(regulated_col, ", "))
  all_entries <- trimws(all_entries)
  
  # Count unique up and down regulated genes
  up_genes <- unique(sub("\\(.*\\)", "", all_entries[grepl("\\(Up\\)$", all_entries)]))
  down_genes <- unique(sub("\\(.*\\)", "", all_entries[grepl("\\(Down\\)$", all_entries)]))
  
  n_up <- length(up_genes)
  n_down <- length(down_genes)
  
  return(c(Up = n_up, Down = n_down))
}

# ==============================================================================
# Main Loop: Process Each Comparison
# ==============================================================================
for (comp in comps) {
  message("\nProcessing comparison: ", comp)
  
  # ==============================================================================
  # Load RNA-seq Data
  # ==============================================================================
  # Map EM-seq comparison (IR...) to RNA-seq comparison
  # Split by "_vs_", remove "IR" prefix from each part if present, then rejoin
  parts <- strsplit(comp, "_vs_")[[1]]
  clean_parts <- gsub("^IR", "", parts)
  rnaseq_comp <- paste(clean_parts, collapse = "_vs_")
  
  rnaseq_file <- file.path(rnaseq_base_dir, rnaseq_comp, paste0(rnaseq_comp, "_DESeq2_results_regular.csv"))
  
  regulated_genes_map <- c()
  if (file.exists(rnaseq_file)) {
    message("  Loading RNA-seq results: ", rnaseq_comp)
    rna_data <- fread(rnaseq_file)
    
    # Filter for significant genes
    sig_rna <- rna_data[padj < rnaseq_padj_cutoff & abs(log2FoldChange) > rnaseq_lfc_cutoff]
    
    if (nrow(sig_rna) > 0) {
      sig_rna$Regulation <- ifelse(sig_rna$log2FoldChange > 0, "Up", "Down")
      sig_rna <- unique(sig_rna[, c("gene_name", "Regulation")])
      regulated_genes_map <- setNames(sig_rna$Regulation, sig_rna$gene_name)
    }
  } else {
    message("  RNA-seq file not found: ", rnaseq_file)
  }
  
  # ==============================================================================
  # Load and Process dmrseq Results
  # ==============================================================================
  dmrseq_file <- file.path(dmrseq_base_dir, comp, "DMR", "DMR_table.tsv")
  if (!file.exists(dmrseq_file)) {
    message("  DMR file not found, skipping comparison: ", comp)
    next
  }
  
  message("  Loading dmrseq results from DMR_table.tsv")
  dmrs <- fread(dmrseq_file, sep = "\t")
  
  # Create GRanges for DMRs
  dmr_gr <- GRanges(
    seqnames = dmrs$seqnames,
    ranges = IRanges(dmrs$start, dmrs$end),
    strand = "*",
    beta = dmrs$beta,
    pval = dmrs$pval
  )
  
  # Find overlapping promoters
  promoter_gr <- GRanges(
    seqnames = seqnames(promoters),
    ranges = IRanges(start(promoters), end(promoters)),
    strand = strand(promoters),
    name = promoters$name
  )
  
  overlaps <- findOverlaps(dmr_gr, promoter_gr)
  
  if (length(overlaps) == 0) {
    message("    No overlapping promoters found for dmrseq")
    dmrseq_hyper_promoters <- data.table()
    dmrseq_hypo_promoters <- data.table()
  } else {
    # Create a data.table with overlapping results
    overlap_dt <- data.table(
      dmr_index = queryHits(overlaps),
      promoter_index = subjectHits(overlaps)
    )
    
    # Merge with DMR and promoter data
    dmr_subset <- as.data.frame(dmrs)[overlap_dt$dmr_index, c("seqnames", "start", "end", "beta", "pval")]
    prom_subset <- as.data.frame(promoters)[overlap_dt$promoter_index, "name", drop = FALSE]
    
    # Ensure all subsets have the same number of rows
    if (nrow(dmr_subset) == nrow(overlap_dt) && nrow(prom_subset) == nrow(overlap_dt)) {
      overlap_dt <- cbind(overlap_dt, dmr_subset, prom_subset)
    } else {
      message("    Error: Mismatch in row counts during merge. Skipping this comparison.")
      next
    }
    
    # Filter significant DMRs based on p-value and beta threshold
    sig_dmrs <- overlap_dt[pval < dmrseq_pval_cutoff & abs(beta) > dmrseq_beta_cutoff]
    
    # Classify as hyper or hypo methylated
    dmrseq_hyper_promoters <- sig_dmrs[sig_dmrs$beta > dmrseq_beta_cutoff, c("seqnames", "start", "end", "name")]
    colnames(dmrseq_hyper_promoters) <- c("chrom", "start", "end", "name")
    dmrseq_hyper_promoters$source <- "dmrseq"
    dmrseq_hypo_promoters <- sig_dmrs[sig_dmrs$beta < -dmrseq_beta_cutoff, c("seqnames", "start", "end", "name")]
    colnames(dmrseq_hypo_promoters) <- c("chrom", "start", "end", "name")
    dmrseq_hypo_promoters$source <- "dmrseq"
    
    message("    dmrseq Hyper Promoters: ", nrow(dmrseq_hyper_promoters))
    message("    dmrseq Hypo Promoters:  ", nrow(dmrseq_hypo_promoters))
  }
  
  # ==============================================================================
  # Load and Process Targeted Results
  # ==============================================================================
  message("  Loading targeted results...")
  targeted_files <- list.files(targeted_base_dir, pattern = paste0("*", targeted_suffix), recursive = TRUE)
  targeted_files <- targeted_files[grepl(comp, targeted_files)]
  
  targeted_hyper_promoters <- data.table()
  targeted_hypo_promoters <- data.table()
  
  if (length(targeted_files) > 0) {
    for (tf in targeted_files) {
      targeted_res <- fread(file.path(targeted_base_dir, tf), sep = "\t")
      
      # Check if required columns exist
      req_cols <- c("meth.diff", "pvalue")
      if (!all(req_cols %in% colnames(targeted_res))) {
        message("    Skipping file (missing columns): ", tf)
        next
      }
      
      # Filter based on thresholds
      sig_targeted <- targeted_res[pvalue < targeted_pval_cutoff & abs(meth.diff) > targeted_diff_cutoff]
      
      # Classify as hyper or hypo methylated
      hyper <- sig_targeted[sig_targeted$meth.diff > targeted_diff_cutoff, c("chr", "start", "end", "name")]
      colnames(hyper) <- c("chrom", "start", "end", "name")
      hyper$source <- "targeted"
      hypo <- sig_targeted[sig_targeted$meth.diff < -targeted_diff_cutoff, c("chr", "start", "end", "name")]
      colnames(hypo) <- c("chrom", "start", "end", "name")
      hypo$source <- "targeted"
      
      targeted_hyper_promoters <- rbind(targeted_hyper_promoters, hyper)
      targeted_hypo_promoters <- rbind(targeted_hypo_promoters, hypo)
    }
    
    # Remove duplicates
    targeted_hyper_promoters <- unique(targeted_hyper_promoters)
    targeted_hypo_promoters <- unique(targeted_hypo_promoters)
    
    message("    Targeted Hyper Promoters: ", nrow(targeted_hyper_promoters))
    message("    Targeted Hypo Promoters:  ", nrow(targeted_hypo_promoters))
  } else {
    message("    No targeted files found for comparison: ", comp)
  }
  
  # ==============================================================================
  # Create Unions and Annotate
  # ==============================================================================
  # Combine dmrseq and targeted results
  # Ensure both data frames have the same structure
  required_cols <- c("chrom", "start", "end", "name", "source")
  
  # Add missing columns to dmrseq results if needed
  if (nrow(dmrseq_hyper_promoters) > 0) {
    for (col in required_cols) {
      if (!col %in% colnames(dmrseq_hyper_promoters)) {
        dmrseq_hyper_promoters[[col]] <- NA_character_
      }
    }
    dmrseq_hyper_promoters <- dmrseq_hyper_promoters[, ..required_cols]
  }
  
  if (nrow(dmrseq_hypo_promoters) > 0) {
    for (col in required_cols) {
      if (!col %in% colnames(dmrseq_hypo_promoters)) {
        dmrseq_hypo_promoters[[col]] <- NA_character_
      }
    }
    dmrseq_hypo_promoters <- dmrseq_hypo_promoters[, ..required_cols]
  }
  
  # Add missing columns to targeted results if needed
  if (nrow(targeted_hyper_promoters) > 0) {
    for (col in required_cols) {
      if (!col %in% colnames(targeted_hyper_promoters)) {
        targeted_hyper_promoters[[col]] <- NA_character_
      }
    }
    targeted_hyper_promoters <- targeted_hyper_promoters[, ..required_cols]
  }
  
  if (nrow(targeted_hypo_promoters) > 0) {
    for (col in required_cols) {
      if (!col %in% colnames(targeted_hypo_promoters)) {
        targeted_hypo_promoters[[col]] <- NA_character_
      }
    }
    targeted_hypo_promoters <- targeted_hypo_promoters[, ..required_cols]
  }
  
  hyper_combined <- unique(rbind(dmrseq_hyper_promoters, targeted_hyper_promoters))
  hypo_combined <- unique(rbind(dmrseq_hypo_promoters, targeted_hypo_promoters))
  
  # Get promoter coordinates for merging
  prom_coords <- unique(promoter_dt[, c("chrom", "start", "end", "name")])
  setkey(prom_coords, name)
  
  # Annotate hypermethylated promoters
  if (nrow(hyper_combined) > 0) {
    setkey(hyper_combined, name)
    hyper_combined <- merge(hyper_combined, prom_coords, by = "name", all.x = TRUE)
    hyper_combined <- merge(hyper_combined, gene_map, by = "name", all.x = TRUE)
    if (!"GeneName" %in% colnames(hyper_combined)) hyper_combined$GeneName <- NA_character_
    hyper_combined$Regulated_TargetGenes <- sapply(hyper_combined$GeneName, annotate_rna, reg_map = regulated_genes_map)
    
    # Use the correct column names after merge (prom_coords columns override original ones)
    hyper_combined <- hyper_combined[, c("chrom.y", "start.y", "end.y", "name", "source", "GeneName", "Regulated_TargetGenes"), with = FALSE]
    # Rename columns back to expected names
    setnames(hyper_combined, c("chrom.y", "start.y", "end.y"), c("chrom", "start", "end"))
    
    # Save hypermethylated promoters
    out_file_hyper <- file.path(output_dir, paste0(comp, "_Hypermethylated_Promoters_Union.csv"))
    fwrite(hyper_combined, out_file_hyper, row.names = FALSE)
  }
  
  # Annotate hypomethylated promoters
  if (nrow(hypo_combined) > 0) {
    setkey(hypo_combined, name)
    hypo_combined <- merge(hypo_combined, prom_coords, by = "name", all.x = TRUE)
    hypo_combined <- merge(hypo_combined, gene_map, by = "name", all.x = TRUE)
    if (!"GeneName" %in% colnames(hypo_combined)) hypo_combined$GeneName <- NA_character_
    hypo_combined$Regulated_TargetGenes <- sapply(hypo_combined$GeneName, annotate_rna, reg_map = regulated_genes_map)
    
    # Use the correct column names after merge (prom_coords columns override original ones)
    hypo_combined <- hypo_combined[, c("chrom.y", "start.y", "end.y", "name", "source", "GeneName", "Regulated_TargetGenes"), with = FALSE]
    # Rename columns back to expected names
    setnames(hypo_combined, c("chrom.y", "start.y", "end.y"), c("chrom", "start", "end"))
    
    # Save hypomethylated promoters
    out_file_hypo <- file.path(output_dir, paste0(comp, "_Hypomethylated_Promoters_Union.csv"))
    fwrite(hypo_combined, out_file_hypo, row.names = FALSE)
  }
  
  # Count regulated genes
  hyper_counts <- count_reg_genes(hyper_combined)
  hypo_counts <- count_reg_genes(hypo_combined)

  # ============================================================================
  # Fisher's Exact Test for Gene Regulation & Promoter Methylation
  # ============================================================================
  # Only run if there are regulated genes in either promoter set
  if (hypo_counts["Up"] > 0 || hypo_counts["Down"] > 0 || hyper_counts["Up"] > 0 || hyper_counts["Down"] > 0) {
    
    # The counts from summary_stats are based on unique genes, which is what we need.
    # Use these counts directly for the contingency table.
    # Matrix is filled by row to match the desired table structure.
    a <- hypo_counts["Up"]    # Hypo & Up
    b <- hypo_counts["Down"]  # Hypo & Down
    c <- hyper_counts["Up"]   # Hyper & Up
    d <- hyper_counts["Down"] # Hyper & Down
    
    # Ensure all counts are non-negative
    counts <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE, # Fill by row
                     dimnames = list(Promoter_Methylation = c("Hypomethylated", "Hypermethylated"),
                                     Gene_Regulation = c("Upregulated", "Downregulated")))
    
    # Perform Fisher's exact test
    fisher_res <- fisher.test(counts)
    
    # Write results to a text file
    out_file_fisher <- file.path(output_dir, paste0(comp, "_Fisher_Exact_Methylation_vs_Regulation.txt"))
    sink(out_file_fisher)
    cat("Fisher's Exact Test 1: Association between Promoter Methylation and Gene Regulation\n")
    cat("Comparison: ", comp, "\n\n")
    cat("Contingency Table:\n")
    print(counts)
    cat("\n")
    cat("P-value: ", fisher_res$p.value, "\n")
    cat("Odds Ratio: ", fisher_res$estimate, "\n")
    cat("95% CI: [", fisher_res$conf.int[1], ", ", fisher_res$conf.int[2], "]\n")
    
    # ============================================================================
    # Second Fisher's Exact Test: Downregulated Genes & Hypermethylated Promoters
    # ============================================================================
    # Construct a 2x2 table focusing on the association between downregulation and hypermethylation
    # Rows: Hypermethylated vs. Not Hypermethylated (i.e., Hypomethylated)
    # Columns: Downregulated vs. Not Downregulated (i.e., Upregulated)
    
    a_hyper_down <- hyper_counts["Down"]   # Hyper & Down
    b_hyper_up   <- hyper_counts["Up"]     # Hyper & Up (Not Down)
    c_hypo_down  <- hypo_counts["Down"]    # Not Hyper & Down
    d_hypo_up    <- hypo_counts["Up"]      # Not Hyper & Not Down (Up)
    
    counts_hyper_down <- matrix(c(a_hyper_down, b_hyper_up, c_hypo_down, d_hypo_up), nrow = 2, byrow = TRUE,
                                dimnames = list(Hypermethylated_Promoter = c("Yes", "No"),
                                                Downregulated_Gene = c("Yes", "No")))
    
    fisher_res_hyper_down <- fisher.test(counts_hyper_down)
    
    cat("\n\n")
    cat("Fisher's Exact Test 2: Association between Downregulated Genes and Hypermethylated Promoters\n")
    cat("Comparison: ", comp, "\n\n")
    cat("Contingency Table:\n")
    print(counts_hyper_down)
    cat("\n")
    cat("P-value: ", fisher_res_hyper_down$p.value, "\n")
    cat("Odds Ratio: ", fisher_res_hyper_down$estimate, "\n")
    cat("95% CI: [", fisher_res_hyper_down$conf.int[1], ", ", fisher_res_hyper_down$conf.int[2], "]\n")
    
    sink()
    
    message("  Fisher's exact test results saved to: ", basename(out_file_fisher))
  }

  # Collect Stats
  summary_stats[[comp]] <- data.frame(
    Comparison = comp,
    Hyper_dmrseq = length(dmrseq_hyper_promoters$name),
    Hyper_Targeted = length(targeted_hyper_promoters$name),
    Hyper_Union = ifelse(is.null(hyper_combined), 0, nrow(hyper_combined)),
    Hyper_Genes_Up = hyper_counts["Up"],
    Hyper_Genes_Down = hyper_counts["Down"],
    
    Hypo_dmrseq = length(dmrseq_hypo_promoters$name),
    Hypo_Targeted = length(targeted_hypo_promoters$name),
    Hypo_Union = ifelse(is.null(hypo_combined), 0, nrow(hypo_combined)),
    Hypo_Genes_Up = hypo_counts["Up"],
    Hypo_Genes_Down = hypo_counts["Down"]
  )
}

# ==============================================================================
# Save Summary Stats
# ==============================================================================
if (length(summary_stats) > 0) {
  final_stats <- do.call(rbind, summary_stats)
  print(final_stats)
  fwrite(final_stats, file.path(output_dir, "Promoter_Methylation_Summary_Counts.csv"))
} else {
  message("No stats collected.")
}

message("\nAnalysis Complete. Results saved to: ", output_dir)
