
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
enhancer_bed_file <- "C:/PROJECTS/resource/T2T_CHM13/ENCFF912FUA_MCF10A_element_gene_links_thresholded_Engreitz_T2T_6col.bed"
enhancer_gene_file <- "C:/PROJECTS/resource/T2T_CHM13/ENCFF912FUA_MCF10A_element_gene_links_thresholded_Engreitz_T2T.bed"
dmrseq_base_dir   <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/dmrseq"
targeted_base_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/Enhancer_targeted"
rnaseq_base_dir   <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/RNAseq/deseq2/single_factor_analysis_gene_level"

# Output Path
output_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/Enhancer_methylation_summary"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Targeted file suffix pattern
targeted_suffix <- "significant_diff_ENCFF912FUA_MCF10A_element_gene_links_thresholded_Engreitz_T2T_6col.tab"

# Filtering Thresholds
dmrseq_pval_cutoff <- 0.05
dmrseq_beta_cutoff <- 0.5

targeted_pval_cutoff <- 0.05
targeted_diff_cutoff <- 10

rnaseq_padj_cutoff <- 0.05
rnaseq_lfc_cutoff <- 1

# ==============================================================================
# Load Enhancer Annotations & Gene Links
# ==============================================================================
message("Loading Enhancer Annotations...")
# Expecting standard 6-column BED: chrom, start, end, name, score, strand
enhancers <- import(enhancer_bed_file, format = "BED")
# Ensure names are unique or usable as IDs
if (any(duplicated(enhancers$name))) {
  warning("Duplicate enhancer names found. Using unique coordinates as fallback IDs if needed.")
}

# Convert to data.table for easier lookup
enhancers_dt <- as.data.table(enhancers)
setkey(enhancers_dt, name)

message("Loading Enhancer-Gene Links...")
# Load the full file to get TargetGene info
links_dt <- fread(enhancer_gene_file)
# Expecting columns 'name' and 'TargetGene'
if (all(c("name", "TargetGene") %in% colnames(links_dt))) {
  # Aggregate genes by enhancer name
  gene_map <- links_dt[, .(TargetGenes = paste(unique(TargetGene), collapse = ", ")), by = name]
  setkey(gene_map, name)
} else {
  warning("Could not find 'name' or 'TargetGene' columns in link file.")
  gene_map <- data.table(name = character(), TargetGenes = character())
}

# ==============================================================================
# Process Comparisons
# ==============================================================================

# Get list of comparisons from dmrseq directory
comparisons <- list.dirs(dmrseq_base_dir, full.names = FALSE, recursive = FALSE)
# Filter for actual comparison folders (contain "_vs_")
comparisons <- comparisons[grepl("_vs_", comparisons)]
# Filter out files or non-comparison folders if any (simple check: must be in targeted dir too usually)
comparisons <- comparisons[dir.exists(file.path(targeted_base_dir, comparisons))]

message("Found ", length(comparisons), " comparisons: ", paste(comparisons, collapse = ", "))

summary_stats <- list()

for (comp in comparisons) {
  message("\nProcessing comparison: ", comp)
  
  # ----------------------------------------------------------------------------
  # Prepare RNA-seq Map
  # ----------------------------------------------------------------------------
  # Map EM-seq comparison (IR...) to RNA-seq comparison
  # Split by "_vs_", remove "IR" prefix from each part if present, then rejoin
  parts <- strsplit(comp, "_vs_")[[1]]
  clean_parts <- gsub("^IR", "", parts)
  rnaseq_comp <- paste(clean_parts, collapse = "_vs_")
  
  rnaseq_file <- file.path(rnaseq_base_dir, rnaseq_comp, paste0(rnaseq_comp, "_DESeq2_results_regular.csv"))
  
  regulated_genes_map <- character()
  
  if (file.exists(rnaseq_file)) {
    message("  Loading RNA-seq results: ", rnaseq_comp)
    res_rna <- fread(rnaseq_file)
    # Check columns
    if (all(c("padj", "log2FoldChange", "gene_name") %in% colnames(res_rna))) {
        # Filter
        sig_rna <- res_rna[padj < rnaseq_padj_cutoff & abs(log2FoldChange) > rnaseq_lfc_cutoff]
        if (nrow(sig_rna) > 0) {
            # Create map: Gene -> "Up" or "Down"
            sig_rna[, Regulation := ifelse(log2FoldChange > 0, "Up", "Down")]
            # Handle duplicate gene names if any by taking first or distinct (usually unique)
            sig_rna <- unique(sig_rna[, .(gene_name, Regulation)])
            regulated_genes_map <- setNames(sig_rna$Regulation, sig_rna$gene_name)
        }
    } else {
        warning("  RNA-seq file missing required columns: ", rnaseq_file)
    }
  } else {
    warning("  RNA-seq file not found: ", rnaseq_file)
  }
  
  # Helper function to annotate RNA regulation
  annotate_rna <- function(target_genes_str, reg_map) {
    if (is.na(target_genes_str) || target_genes_str == "") return(NA_character_)
    genes <- unlist(strsplit(target_genes_str, ", "))
    genes <- trimws(genes)
    
    # Check if genes are in the map
    reg_status <- reg_map[genes]
    # Filter NAs
    hits <- !is.na(reg_status)
    
    if (!any(hits)) return(NA_character_)
    
    # Format: "GeneA(Up), GeneB(Down)"
    return(paste(paste0(genes[hits], "(", reg_status[hits], ")"), collapse = ", "))
  }

  # ----------------------------------------------------------------------------
  # 1. Process dmrseq Results
  # ----------------------------------------------------------------------------
  dmrseq_hyper_enhancers <- character()
  dmrseq_hypo_enhancers  <- character()
  
  # Try loading DMR table (tsv preferred) or annotated xlsx
  dmr_table_path <- file.path(dmrseq_base_dir, comp, "DMR", "DMR_table.tsv")
  
  if (file.exists(dmr_table_path)) {
    message("  Loading dmrseq results from DMR_table.tsv")
    dmrs <- fread(dmr_table_path)
    
    # Filter significant DMRs based on p-value and beta threshold
    sig_dmrs <- dmrs[pval < dmrseq_pval_cutoff & abs(beta) > dmrseq_beta_cutoff]
    
    if (nrow(sig_dmrs) > 0) {
      # Split into Hyper and Hypo
      hyper_dmrs <- sig_dmrs[beta > dmrseq_beta_cutoff]
      hypo_dmrs  <- sig_dmrs[beta < -dmrseq_beta_cutoff]
      
      # Convert to GRanges
      if (nrow(hyper_dmrs) > 0) {
        hyper_gr <- makeGRangesFromDataFrame(hyper_dmrs, keep.extra.columns = TRUE)
        overlaps <- findOverlaps(enhancers, hyper_gr)
        dmrseq_hyper_enhancers <- unique(enhancers$name[queryHits(overlaps)])
      }
      
      if (nrow(hypo_dmrs) > 0) {
        hypo_gr <- makeGRangesFromDataFrame(hypo_dmrs, keep.extra.columns = TRUE)
        overlaps <- findOverlaps(enhancers, hypo_gr)
        dmrseq_hypo_enhancers <- unique(enhancers$name[queryHits(overlaps)])
      }
    }
  } else {
    warning("  dmrseq table not found for ", comp)
  }
  
  message("    dmrseq Hyper Enhancers: ", length(dmrseq_hyper_enhancers))
  message("    dmrseq Hypo Enhancers:  ", length(dmrseq_hypo_enhancers))
  
  # ----------------------------------------------------------------------------
  # 2. Process Targeted Results
  # ----------------------------------------------------------------------------
  targeted_hyper_enhancers <- character()
  targeted_hypo_enhancers  <- character()
  
  targeted_file <- file.path(targeted_base_dir, comp, "Targeted", targeted_suffix)
  
  if (file.exists(targeted_file)) {
    message("  Loading targeted results...")
    # These files usually contain ONLY significant differences
    targeted_res <- fread(targeted_file)
    
    if (nrow(targeted_res) > 0) {
      # Filter based on thresholds
      sig_targeted <- targeted_res[pvalue < targeted_pval_cutoff & abs(meth.diff) > targeted_diff_cutoff]
      
      # Hyper: meth.diff > 0
      targeted_hyper_enhancers <- unique(sig_targeted[meth.diff > 0]$name)
      # Hypo: meth.diff < 0
      targeted_hypo_enhancers  <- unique(sig_targeted[meth.diff < 0]$name)
    }
  } else {
    warning("  Targeted results file not found for ", comp)
  }
  
  message("    Targeted Hyper Enhancers: ", length(targeted_hyper_enhancers))
  message("    Targeted Hypo Enhancers:  ", length(targeted_hypo_enhancers))
  
  # ----------------------------------------------------------------------------
  # 3. Union and Classify
  # ----------------------------------------------------------------------------
  
  # Function to combine and classify
  combine_results <- function(dmrseq_ids, targeted_ids) {
    all_ids <- unique(c(dmrseq_ids, targeted_ids))
    if (length(all_ids) == 0) return(NULL)
    
    res <- data.frame(name = all_ids, source = "NA", stringsAsFactors = FALSE)
    
    in_dmrseq <- all_ids %in% dmrseq_ids
    in_targeted <- all_ids %in% targeted_ids
    
    res$source[in_dmrseq & in_targeted] <- "Both"
    res$source[in_dmrseq & !in_targeted] <- "dmrseq_only"
    res$source[!in_dmrseq & in_targeted] <- "Targeted_only"
    
    return(res)
  }
  
  # Hypermethylated Union
  hyper_combined <- combine_results(dmrseq_hyper_enhancers, targeted_hyper_enhancers)
  if (!is.null(hyper_combined)) {
    # Add coordinates - Ensure we only take atomic columns from enhancers_dt
    # Also ensure name is character to match
    enh_coords <- enhancers_dt[, .(chrom = as.character(seqnames), start, end, name = as.character(name))]
    
    hyper_combined <- merge(hyper_combined, enh_coords, by = "name", all.x = TRUE)
    # Add Target Genes
    hyper_combined <- merge(hyper_combined, gene_map, by = "name", all.x = TRUE)
    
    # Reorder columns and ensure no list columns
    # Handle missing TargetGenes if any
    if (!"TargetGenes" %in% colnames(hyper_combined)) hyper_combined$TargetGenes <- NA_character_
    
    # Annotate Regulated Genes
    hyper_combined$Regulated_TargetGenes <- sapply(hyper_combined$TargetGenes, annotate_rna, reg_map = regulated_genes_map)
    
    hyper_combined <- hyper_combined[, c("chrom", "start", "end", "name", "source", "TargetGenes", "Regulated_TargetGenes")]
    
    # Save
    out_file_hyper <- file.path(output_dir, paste0(comp, "_Hypermethylated_Enhancers_Union.csv"))
    fwrite(hyper_combined, out_file_hyper)
    
    # Create BED file
    hyper_bed <- hyper_combined
    hyper_bed$name <- paste0(hyper_bed$name, "|", hyper_bed$source)
    hyper_bed$score <- 0
    hyper_bed$strand <- "*"
    
    # Ensure all columns are atomic before writing
    hyper_bed_df <- as.data.frame(hyper_bed)
    fwrite(hyper_bed_df[, c("chrom", "start", "end", "name", "score", "strand")], 
           file.path(output_dir, paste0(comp, "_Hypermethylated_Enhancers_Union.bed")), 
           sep = "\t", col.names = FALSE)
  }
  
  # Hypomethylated Union
  hypo_combined <- combine_results(dmrseq_hypo_enhancers, targeted_hypo_enhancers)
  if (!is.null(hypo_combined)) {
    # Add coordinates
    enh_coords <- enhancers_dt[, .(chrom = as.character(seqnames), start, end, name = as.character(name))]
    
    hypo_combined <- merge(hypo_combined, enh_coords, by = "name", all.x = TRUE)
    # Add Target Genes
    hypo_combined <- merge(hypo_combined, gene_map, by = "name", all.x = TRUE)
    
    # Reorder columns
    if (!"TargetGenes" %in% colnames(hypo_combined)) hypo_combined$TargetGenes <- NA_character_
    
    # Annotate Regulated Genes
    hypo_combined$Regulated_TargetGenes <- sapply(hypo_combined$TargetGenes, annotate_rna, reg_map = regulated_genes_map)
    
    hypo_combined <- hypo_combined[, c("chrom", "start", "end", "name", "source", "TargetGenes", "Regulated_TargetGenes")]
    
    # Save
    out_file_hypo <- file.path(output_dir, paste0(comp, "_Hypomethylated_Enhancers_Union.csv"))
    fwrite(hypo_combined, out_file_hypo)
    
    # BED
    hypo_bed <- hypo_combined
    hypo_bed$name <- paste0(hypo_bed$name, "|", hypo_bed$source)
    hypo_bed$score <- 0
    hypo_bed$strand <- "*"
    
    hypo_bed_df <- as.data.frame(hypo_bed)
    fwrite(hypo_bed_df[, c("chrom", "start", "end", "name", "score", "strand")], 
           file.path(output_dir, paste0(comp, "_Hypomethylated_Enhancers_Union.bed")), 
           sep = "\t", col.names = FALSE)
  }
  
  # Collect Stats
  summary_stats[[comp]] <- data.frame(
    Comparison = comp,
    Hyper_dmrseq = length(dmrseq_hyper_enhancers),
    Hyper_Targeted = length(targeted_hyper_enhancers),
    Hyper_Union = ifelse(is.null(hyper_combined), 0, nrow(hyper_combined)),
    Hypo_dmrseq = length(dmrseq_hypo_enhancers),
    Hypo_Targeted = length(targeted_hypo_enhancers),
    Hypo_Union = ifelse(is.null(hypo_combined), 0, nrow(hypo_combined))
  )
}

# ==============================================================================
# Save Summary Stats
# ==============================================================================
if (length(summary_stats) > 0) {
  final_stats <- do.call(rbind, summary_stats)
  print(final_stats)
  fwrite(final_stats, file.path(output_dir, "Enhancer_Methylation_Summary_Counts.csv"))
} else {
  message("No stats collected.")
}

message("\nAnalysis Complete. Results saved to: ", output_dir)
