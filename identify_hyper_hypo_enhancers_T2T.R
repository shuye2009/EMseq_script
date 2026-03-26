
# ==============================================================================
# Package Installation and Loading
# ==============================================================================

# Function to check and install CRAN packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message("Installing package: ", pkg)
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Function to check and install Bioconductor packages
install_bioc_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message("Installing Bioconductor package: ", pkg)
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Install and load CRAN packages
cran_packages <- c("data.table", "openxlsx", "dplyr", "tidyr", "stringr")
for (pkg in cran_packages) {
  install_if_missing(pkg)
}

# Install and load Bioconductor packages
bioc_packages <- c("GenomicRanges", "rtracklayer", "clusterProfiler", "org.Hs.eg.db", "mygene")
for (pkg in bioc_packages) {
  install_bioc_if_missing(pkg)
}

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
targeted_suffix <- "all_diff_ENCFF912FUA_MCF10A_element_gene_links_thresholded_Engreitz_T2T_6col.tab"

# Filtering Thresholds
dmrseq_pval_cutoff <- 0.05
dmrseq_beta_cutoff <- 0.25

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
message("  Total enhancers loaded: ", length(enhancers))
message("  Sample enhancer names from BED: ", paste(head(enhancers$name, 3), collapse = ", "))
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
  dmrseq_hyper_dmr_info <- list()
  dmrseq_hypo_dmr_info <- list()
  
  # Try loading DMR table (tsv preferred) or annotated xlsx
  dmr_table_path <- file.path(dmrseq_base_dir, comp, "DMRs", "DMRs_annotated.xlsx")
  
  if (file.exists(dmr_table_path)) {
    message("  Loading dmrseq results from DMRs_annotated.xlsx")
    dmrs <- as.data.table(read.xlsx(dmr_table_path))
    
    # Filter significant DMRs based on p-value and beta threshold
    sig_dmrs <- dmrs[p.value < dmrseq_pval_cutoff & abs(betaCoefficient) > dmrseq_beta_cutoff]
    
    if (nrow(sig_dmrs) > 0) {
      # Split into Hyper and Hypo
      hyper_dmrs <- sig_dmrs[betaCoefficient > dmrseq_beta_cutoff]
      hypo_dmrs  <- sig_dmrs[betaCoefficient < -dmrseq_beta_cutoff]
      
      # Convert to GRanges
      if (nrow(hyper_dmrs) > 0) {
        hyper_gr <- makeGRangesFromDataFrame(hyper_dmrs, keep.extra.columns = TRUE)
        overlaps <- findOverlaps(enhancers, hyper_gr, minoverlap = min(width(hyper_gr)))
        
        if (length(overlaps) > 0) {
          # Calculate overlap widths
          enh_hits <- enhancers[queryHits(overlaps)]
          dmr_hits <- hyper_gr[subjectHits(overlaps)]
          overlap_ranges <- pintersect(enh_hits, dmr_hits)
          overlap_widths <- width(overlap_ranges)
          
          # Store DMR info for each enhancer
          for (i in seq_along(overlaps)) {
            enh_name <- enhancers$name[queryHits(overlaps)[i]]
            dmr_chr <- as.character(seqnames(dmr_hits[i]))
            dmr_start <- start(dmr_hits[i])
            dmr_end <- end(dmr_hits[i])
            dmr_width <- width(dmr_hits[i])
            ovlp_width <- overlap_widths[i]
            
            dmr_string <- paste0(dmr_chr, ":", dmr_start, "-", dmr_end)
            
            if (is.null(dmrseq_hyper_dmr_info[[enh_name]])) {
              dmrseq_hyper_dmr_info[[enh_name]] <- list(
                dmr_coords = dmr_string,
                dmr_widths = dmr_width,
                overlap_widths = ovlp_width
              )
            } else {
              dmrseq_hyper_dmr_info[[enh_name]]$dmr_coords <- c(
                dmrseq_hyper_dmr_info[[enh_name]]$dmr_coords, dmr_string
              )
              dmrseq_hyper_dmr_info[[enh_name]]$dmr_widths <- c(
                dmrseq_hyper_dmr_info[[enh_name]]$dmr_widths, dmr_width
              )
              dmrseq_hyper_dmr_info[[enh_name]]$overlap_widths <- c(
                dmrseq_hyper_dmr_info[[enh_name]]$overlap_widths, ovlp_width
              )
            }
          }
          
          dmrseq_hyper_enhancers <- unique(enhancers$name[queryHits(overlaps)])
        }
      }
      
      if (nrow(hypo_dmrs) > 0) {
        hypo_gr <- makeGRangesFromDataFrame(hypo_dmrs, keep.extra.columns = TRUE)
        overlaps <- findOverlaps(enhancers, hypo_gr, minoverlap = min(width(hypo_gr)))
        
        if (length(overlaps) > 0) {
          # Calculate overlap widths
          enh_hits <- enhancers[queryHits(overlaps)]
          dmr_hits <- hypo_gr[subjectHits(overlaps)]
          overlap_ranges <- pintersect(enh_hits, dmr_hits)
          overlap_widths <- width(overlap_ranges)
          
          # Store DMR info for each enhancer
          for (i in seq_along(overlaps)) {
            enh_name <- enhancers$name[queryHits(overlaps)[i]]
            dmr_chr <- as.character(seqnames(dmr_hits[i]))
            dmr_start <- start(dmr_hits[i])
            dmr_end <- end(dmr_hits[i])
            dmr_width <- width(dmr_hits[i])
            ovlp_width <- overlap_widths[i]
            
            dmr_string <- paste0(dmr_chr, ":", dmr_start, "-", dmr_end)
            
            if (is.null(dmrseq_hypo_dmr_info[[enh_name]])) {
              dmrseq_hypo_dmr_info[[enh_name]] <- list(
                dmr_coords = dmr_string,
                dmr_widths = dmr_width,
                overlap_widths = ovlp_width
              )
            } else {
              dmrseq_hypo_dmr_info[[enh_name]]$dmr_coords <- c(
                dmrseq_hypo_dmr_info[[enh_name]]$dmr_coords, dmr_string
              )
              dmrseq_hypo_dmr_info[[enh_name]]$dmr_widths <- c(
                dmrseq_hypo_dmr_info[[enh_name]]$dmr_widths, dmr_width
              )
              dmrseq_hypo_dmr_info[[enh_name]]$overlap_widths <- c(
                dmrseq_hypo_dmr_info[[enh_name]]$overlap_widths, ovlp_width
              )
            }
          }
          
          dmrseq_hypo_enhancers <- unique(enhancers$name[queryHits(overlaps)])
        }
      }
    }
  } else {
    warning("  dmrseq table not found for ", comp)
  }
  
  message("    dmrseq Hyper Enhancers: ", length(dmrseq_hyper_enhancers))
  if (length(dmrseq_hyper_enhancers) > 0) {
    message("      Sample dmrseq Hyper names: ", paste(head(dmrseq_hyper_enhancers, 3), collapse = ", "))
    message("      Class: ", class(dmrseq_hyper_enhancers), ", Type: ", typeof(dmrseq_hyper_enhancers))
  }
  message("    dmrseq Hypo Enhancers:  ", length(dmrseq_hypo_enhancers))
  if (length(dmrseq_hypo_enhancers) > 0) {
    message("      Sample dmrseq Hypo names: ", paste(head(dmrseq_hypo_enhancers, 3), collapse = ", "))
  }
  
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
    message("    Total rows in targeted file: ", nrow(targeted_res))
    if (nrow(targeted_res) > 0) {
      message("    Columns in targeted file: ", paste(colnames(targeted_res), collapse = ", "))
      message("    Sample 'name' values from targeted file: ", paste(head(targeted_res$name, 3), collapse = ", "))
    }
    
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
  if (length(targeted_hyper_enhancers) > 0) {
    message("      Sample Targeted Hyper names: ", paste(head(targeted_hyper_enhancers, 3), collapse = ", "))
    message("      Class: ", class(targeted_hyper_enhancers), ", Type: ", typeof(targeted_hyper_enhancers))
    
    # Check if any dmrseq names exist in the targeted list
    if (length(dmrseq_hyper_enhancers) > 0) {
      check_names <- head(dmrseq_hyper_enhancers, 5)
      matches <- check_names %in% targeted_hyper_enhancers
      message("      Checking if dmrseq names exist in targeted: ", paste(check_names, "=", matches, collapse = "; "))
    }
  }
  message("    Targeted Hypo Enhancers:  ", length(targeted_hypo_enhancers))
  if (length(targeted_hypo_enhancers) > 0) {
    message("      Sample Targeted Hypo names: ", paste(head(targeted_hypo_enhancers, 3), collapse = ", "))
  }
  
  # ----------------------------------------------------------------------------
  # 3. Union and Classify
  # ----------------------------------------------------------------------------
  
  # Function to combine and classify using coordinate-based matching
  combine_results <- function(dmrseq_ids, targeted_ids) {
    all_ids <- unique(c(dmrseq_ids, targeted_ids))
    if (length(all_ids) == 0) return(NULL)
    
    # Parse coordinates from name strings (format: chr:start-end)
    parse_coords <- function(name_vec) {
      if (length(name_vec) == 0) return(NULL)
      parts <- strsplit(name_vec, "[:-]")
      data.frame(
        name = name_vec,
        chr = sapply(parts, `[`, 1),
        start = as.integer(sapply(parts, `[`, 2)),
        end = as.integer(sapply(parts, `[`, 3)),
        stringsAsFactors = FALSE
      )
    }
    
    dmrseq_coords <- parse_coords(dmrseq_ids)
    targeted_coords <- parse_coords(targeted_ids)
    
    # Create result dataframe
    res <- data.frame(name = all_ids, source = "NA", stringsAsFactors = FALSE)
    
    # For each enhancer, check if it's in dmrseq, targeted, or both
    # Use coordinate overlap with tolerance for off-by-one errors
    for (i in seq_len(nrow(res))) {
      enh_name <- res$name[i]
      
      # Check if in dmrseq
      in_dmrseq <- enh_name %in% dmrseq_ids
      
      # Check if in targeted (with coordinate matching)
      in_targeted <- enh_name %in% targeted_ids
      
      # If not exact match, try coordinate-based matching
      if (!in_dmrseq && !is.null(targeted_coords) && enh_name %in% targeted_coords$name) {
        # This is a targeted enhancer, check if dmrseq has same region
        tgt_row <- targeted_coords[targeted_coords$name == enh_name, ]
        if (!is.null(dmrseq_coords)) {
          # Check for coordinate overlap (allowing ±2bp tolerance)
          matches <- dmrseq_coords$chr == tgt_row$chr & 
                     abs(dmrseq_coords$start - tgt_row$start) <= 2 & 
                     abs(dmrseq_coords$end - tgt_row$end) <= 2
          if (any(matches)) {
            in_dmrseq <- TRUE
            # Use the dmrseq name as canonical
            res$name[i] <- dmrseq_coords$name[which(matches)[1]]
          }
        }
      } else if (!in_targeted && !is.null(dmrseq_coords) && enh_name %in% dmrseq_coords$name) {
        # This is a dmrseq enhancer, check if targeted has same region
        dmr_row <- dmrseq_coords[dmrseq_coords$name == enh_name, ]
        if (!is.null(targeted_coords)) {
          # Check for coordinate overlap (allowing ±2bp tolerance)
          matches <- targeted_coords$chr == dmr_row$chr & 
                     abs(targeted_coords$start - dmr_row$start) <= 2 & 
                     abs(targeted_coords$end - dmr_row$end) <= 2
          if (any(matches)) {
            in_targeted <- TRUE
          }
        }
      }
      
      # Classify
      if (in_dmrseq && in_targeted) {
        res$source[i] <- "Both"
      } else if (in_dmrseq) {
        res$source[i] <- "dmrseq_only"
      } else if (in_targeted) {
        res$source[i] <- "Targeted_only"
      }
    }
    
    return(res)
  }
  
  # Hypermethylated Union
  message("  Combining Hypermethylated results...")
  hyper_combined <- combine_results(dmrseq_hyper_enhancers, targeted_hyper_enhancers)
  if (!is.null(hyper_combined)) {
    message("    Source distribution: ", paste(names(table(hyper_combined$source)), "=", table(hyper_combined$source), collapse = ", "))
    
    # Add coordinates and target genes using coordinate-based matching
    hyper_combined$chrom <- NA_character_
    hyper_combined$start <- NA_integer_
    hyper_combined$end <- NA_integer_
    hyper_combined$TargetGenes <- NA_character_
    
    for (i in seq_len(nrow(hyper_combined))) {
      enh_name <- hyper_combined$name[i]
      
      # Parse coordinates from name string
      parts <- strsplit(enh_name, "[:-]")[[1]]
      if (length(parts) == 3) {
        query_chr <- parts[1]
        query_start <- as.integer(parts[2])
        query_end <- as.integer(parts[3])
        
        # Find matching enhancer in BED file (with ±2bp tolerance)
        matches <- which(
          as.character(enhancers_dt$seqnames) == query_chr &
          abs(enhancers_dt$start - query_start) <= 2 &
          abs(enhancers_dt$end - query_end) <= 2
        )
        
        if (length(matches) > 0) {
          # Use first match
          match_idx <- matches[1]
          hyper_combined$chrom[i] <- as.character(enhancers_dt$seqnames[match_idx])
          hyper_combined$start[i] <- enhancers_dt$start[match_idx]
          hyper_combined$end[i] <- enhancers_dt$end[match_idx]
          
          # Get target genes for this enhancer
          enh_name_bed <- as.character(enhancers_dt$name[match_idx])
          if (enh_name_bed %in% gene_map$name) {
            hyper_combined$TargetGenes[i] <- gene_map[name == enh_name_bed]$TargetGenes
          }
        }
      }
    }
    
    # Reorder columns and ensure no list columns
    # Handle missing TargetGenes if any
    if (!"TargetGenes" %in% colnames(hyper_combined)) hyper_combined$TargetGenes <- NA_character_
    
    # Annotate Regulated Genes
    hyper_combined$Regulated_TargetGenes <- sapply(hyper_combined$TargetGenes, annotate_rna, reg_map = regulated_genes_map)
    
    # Add enhancer width (1-based inclusive coordinates)
    hyper_combined$Enhancer_Width <- hyper_combined$end - hyper_combined$start + 1
    
    # Add DMR coordinates, widths, overlap widths, and fractions for dmrseq_only and Both
    hyper_combined$DMR_Coordinates <- NA_character_
    hyper_combined$DMR_Width <- NA_character_
    hyper_combined$DMR_Overlap_Width <- NA_character_
    hyper_combined$Overlap_Fraction <- NA_real_
    hyper_combined$DMR_Overlap_Fraction <- NA_real_
    
    for (i in seq_len(nrow(hyper_combined))) {
      enh_name <- hyper_combined$name[i]
      src <- hyper_combined$source[i]
      
      if (src %in% c("dmrseq_only", "Both") && enh_name %in% names(dmrseq_hyper_dmr_info)) {
        dmr_info <- dmrseq_hyper_dmr_info[[enh_name]]
        hyper_combined$DMR_Coordinates[i] <- paste(dmr_info$dmr_coords, collapse = "; ")
        hyper_combined$DMR_Width[i] <- paste(dmr_info$dmr_widths, collapse = "; ")
        hyper_combined$DMR_Overlap_Width[i] <- paste(dmr_info$overlap_widths, collapse = "; ")
        
        # Calculate overlap fraction (sum of all overlaps / enhancer width)
        total_overlap <- sum(dmr_info$overlap_widths)
        enh_width <- hyper_combined$Enhancer_Width[i]
        hyper_combined$Overlap_Fraction[i] <- if (enh_width > 0) total_overlap / enh_width else NA_real_
        
        # Calculate DMR overlap fraction (sum of overlap widths / sum of DMR widths)
        total_dmr_width <- sum(dmr_info$dmr_widths)
        hyper_combined$DMR_Overlap_Fraction[i] <- if (total_dmr_width > 0) total_overlap / total_dmr_width else NA_real_
      }
    }
    
    hyper_combined <- hyper_combined[, c("chrom", "start", "end", "name", "source", "TargetGenes", "Regulated_TargetGenes", "Enhancer_Width", "DMR_Coordinates", "DMR_Width", "DMR_Overlap_Width", "Overlap_Fraction", "DMR_Overlap_Fraction")]
    
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
  message("  Combining Hypomethylated results...")
  hypo_combined <- combine_results(dmrseq_hypo_enhancers, targeted_hypo_enhancers)
  if (!is.null(hypo_combined)) {
    message("    Source distribution: ", paste(names(table(hypo_combined$source)), "=", table(hypo_combined$source), collapse = ", "))
    
    # Add coordinates and target genes using coordinate-based matching
    hypo_combined$chrom <- NA_character_
    hypo_combined$start <- NA_integer_
    hypo_combined$end <- NA_integer_
    hypo_combined$TargetGenes <- NA_character_
    
    for (i in seq_len(nrow(hypo_combined))) {
      enh_name <- hypo_combined$name[i]
      
      # Parse coordinates from name string
      parts <- strsplit(enh_name, "[:-]")[[1]]
      if (length(parts) == 3) {
        query_chr <- parts[1]
        query_start <- as.integer(parts[2])
        query_end <- as.integer(parts[3])
        
        # Find matching enhancer in BED file (with ±2bp tolerance)
        matches <- which(
          as.character(enhancers_dt$seqnames) == query_chr &
          abs(enhancers_dt$start - query_start) <= 2 &
          abs(enhancers_dt$end - query_end) <= 2
        )
        
        if (length(matches) > 0) {
          # Use first match
          match_idx <- matches[1]
          hypo_combined$chrom[i] <- as.character(enhancers_dt$seqnames[match_idx])
          hypo_combined$start[i] <- enhancers_dt$start[match_idx]
          hypo_combined$end[i] <- enhancers_dt$end[match_idx]
          
          # Get target genes for this enhancer
          enh_name_bed <- as.character(enhancers_dt$name[match_idx])
          if (enh_name_bed %in% gene_map$name) {
            hypo_combined$TargetGenes[i] <- gene_map[name == enh_name_bed]$TargetGenes
          }
        }
      }
    }
    
    # Reorder columns
    if (!"TargetGenes" %in% colnames(hypo_combined)) hypo_combined$TargetGenes <- NA_character_
    
    # Annotate Regulated Genes
    hypo_combined$Regulated_TargetGenes <- sapply(hypo_combined$TargetGenes, annotate_rna, reg_map = regulated_genes_map)
    
    # Add enhancer width (1-based inclusive coordinates)
    hypo_combined$Enhancer_Width <- hypo_combined$end - hypo_combined$start + 1
    
    # Add DMR coordinates, widths, overlap widths, and fractions for dmrseq_only and Both
    hypo_combined$DMR_Coordinates <- NA_character_
    hypo_combined$DMR_Width <- NA_character_
    hypo_combined$DMR_Overlap_Width <- NA_character_
    hypo_combined$Overlap_Fraction <- NA_real_
    hypo_combined$DMR_Overlap_Fraction <- NA_real_
    
    for (i in seq_len(nrow(hypo_combined))) {
      enh_name <- hypo_combined$name[i]
      src <- hypo_combined$source[i]
      
      if (src %in% c("dmrseq_only", "Both") && enh_name %in% names(dmrseq_hypo_dmr_info)) {
        dmr_info <- dmrseq_hypo_dmr_info[[enh_name]]
        hypo_combined$DMR_Coordinates[i] <- paste(dmr_info$dmr_coords, collapse = "; ")
        hypo_combined$DMR_Width[i] <- paste(dmr_info$dmr_widths, collapse = "; ")
        hypo_combined$DMR_Overlap_Width[i] <- paste(dmr_info$overlap_widths, collapse = "; ")
        
        # Calculate overlap fraction (sum of all overlaps / enhancer width)
        total_overlap <- sum(dmr_info$overlap_widths)
        enh_width <- hypo_combined$Enhancer_Width[i]
        hypo_combined$Overlap_Fraction[i] <- if (enh_width > 0) total_overlap / enh_width else NA_real_
        
        # Calculate DMR overlap fraction (sum of overlap widths / sum of DMR widths)
        total_dmr_width <- sum(dmr_info$dmr_widths)
        hypo_combined$DMR_Overlap_Fraction[i] <- if (total_dmr_width > 0) total_overlap / total_dmr_width else NA_real_
      }
    }
    
    hypo_combined <- hypo_combined[, c("chrom", "start", "end", "name", "source", "TargetGenes", "Regulated_TargetGenes", "Enhancer_Width", "DMR_Coordinates", "DMR_Width", "DMR_Overlap_Width", "Overlap_Fraction", "DMR_Overlap_Fraction")]
    
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
  
  # ----------------------------------------------------------------------------
  # Count Regulated Genes
  # ----------------------------------------------------------------------------
  count_reg_genes <- function(df) {
    if (is.null(df) || !"Regulated_TargetGenes" %in% colnames(df)) return(c(Up = 0, Down = 0))
    
    # Get column, remove NAs and empty strings
    vals <- na.omit(df$Regulated_TargetGenes)
    vals <- vals[vals != ""]
    
    if (length(vals) == 0) return(c(Up = 0, Down = 0))
    
    # Split multiple genes: "GeneA(Up), GeneB(Down)"
    all_entries <- unlist(strsplit(vals, ", "))
    # Unique genes per enhancer set
    unique_entries <- unique(all_entries)
    
    n_up <- sum(grepl("\\(Up\\)$", unique_entries))
    n_down <- sum(grepl("\\(Down\\)$", unique_entries))
    
    return(c(Up = n_up, Down = n_down))
  }
  
  hyper_counts <- count_reg_genes(hyper_combined)
  hypo_counts <- count_reg_genes(hypo_combined)

  # ----------------------------------------------------------------------------
  # Fisher's Exact Test for Gene Regulation & Enhancer Methylation
  # ----------------------------------------------------------------------------
  # Only run if there are regulated genes in either enhancer set
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
                     dimnames = list(Enhancer_Methylation = c("Hypomethylated", "Hypermethylated"),
                                     Gene_Regulation = c("Upregulated", "Downregulated")))
    
    # Perform Fisher's exact test
    fisher_res <- fisher.test(counts)
    
    # Write results to a text file
    out_file_fisher <- file.path(output_dir, paste0(comp, "_Fisher_Exact_Methylation_vs_Regulation.txt"))
    sink(out_file_fisher)
    cat("Fisher's Exact Test 1: Association between Enhancer Methylation and Gene Regulation\n")
    cat("Comparison: ", comp, "\n\n")
    cat("Contingency Table:\n")
    print(counts)
    cat("\n")
    cat("P-value: ", fisher_res$p.value, "\n")
    cat("Odds Ratio: ", fisher_res$estimate, "\n")
    cat("95% CI: [", fisher_res$conf.int[1], ", ", fisher_res$conf.int[2], "]\n")
    
    # ==============================================================================
    # Second Fisher's Exact Test: Downregulated Genes & Hypermethylated Enhancers
    # ==============================================================================
    # Construct a 2x2 table focusing on the association between downregulation and hypermethylation
    # Rows: Hypermethylated vs. Not Hypermethylated (i.e., Hypomethylated)
    # Columns: Downregulated vs. Not Downregulated (i.e., Upregulated)
    
    a_hyper_down <- hyper_counts["Down"]   # Hyper & Down
    b_hyper_up   <- hyper_counts["Up"]     # Hyper & Up (Not Down)
    c_hypo_down  <- hypo_counts["Down"]    # Not Hyper & Down
    d_hypo_up    <- hypo_counts["Up"]      # Not Hyper & Not Down (Up)
    
    counts_hyper_down <- matrix(c(a_hyper_down, b_hyper_up, c_hypo_down, d_hypo_up), nrow = 2, byrow = TRUE,
                                dimnames = list(Enhancer_Methylation = c("Hypermethylated", "Hypomethylated"),
                                                Gene_Regulation = c("Downregulated", "Upregulated")))
    
    fisher_res_hyper_down <- fisher.test(counts_hyper_down)
    
    cat("\n\n")
    cat("Fisher's Exact Test 2: Association between Downregulated Genes and Hypermethylated Enhancers\n")
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

  # ----------------------------------------------------------------------------
  # GO:BP Enrichment Analysis
  # ----------------------------------------------------------------------------
  
  # Extract all target genes from hypermethylated enhancers
  hyper_target_genes <- c()
  hyper_upreg_genes <- c()
  
  if (!is.null(hyper_combined) && "TargetGenes" %in% colnames(hyper_combined)) {
    # Get all target genes
    target_vals <- na.omit(hyper_combined$TargetGenes)
    target_vals <- target_vals[target_vals != ""]
    if (length(target_vals) > 0) {
      all_genes <- unlist(strsplit(target_vals, ","))
      all_genes <- trimws(all_genes)
      hyper_target_genes <- unique(all_genes[all_genes != ""])
    }
    
    # Get upregulated target genes
    if ("Regulated_TargetGenes" %in% colnames(hyper_combined)) {
      reg_vals <- na.omit(hyper_combined$Regulated_TargetGenes)
      reg_vals <- reg_vals[reg_vals != ""]
      if (length(reg_vals) > 0) {
        all_reg <- unlist(strsplit(reg_vals, ", "))
        up_reg <- all_reg[grepl("\\(Up\\)$", all_reg)]
        hyper_upreg_genes <- unique(gsub("\\(Up\\)$", "", up_reg))
      }
    }
  }
  
  # Extract all target genes from hypomethylated enhancers
  hypo_target_genes <- c()
  hypo_upreg_genes <- c()
  
  if (!is.null(hypo_combined) && "TargetGenes" %in% colnames(hypo_combined)) {
    # Get all target genes
    target_vals <- na.omit(hypo_combined$TargetGenes)
    target_vals <- target_vals[target_vals != ""]
    if (length(target_vals) > 0) {
      all_genes <- unlist(strsplit(target_vals, ","))
      all_genes <- trimws(all_genes)
      hypo_target_genes <- unique(all_genes[all_genes != ""])
    }
    
    # Get upregulated target genes
    if ("Regulated_TargetGenes" %in% colnames(hypo_combined)) {
      reg_vals <- na.omit(hypo_combined$Regulated_TargetGenes)
      reg_vals <- reg_vals[reg_vals != ""]
      if (length(reg_vals) > 0) {
        all_reg <- unlist(strsplit(reg_vals, ", "))
        up_reg <- all_reg[grepl("\\(Up\\)$", all_reg)]
        hypo_upreg_genes <- unique(gsub("\\(Up\\)$", "", up_reg))
      }
    }
  }
  
  # Perform GO:BP enrichment for hypermethylated enhancer target genes
  if (length(hyper_target_genes) >= 3) {
    message("  Running GO:BP enrichment for hypermethylated enhancer target genes...")
    tryCatch({
      go_hyper_all <- enrichGO(
        gene = hyper_target_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )
      
      if (!is.null(go_hyper_all) && nrow(as.data.frame(go_hyper_all)) > 0) {
        out_file_go <- file.path(output_dir, paste0(comp, "_Hypermethylated_TargetGenes_GO_BP.csv"))
        fwrite(as.data.frame(go_hyper_all), out_file_go)
        message("    Saved to: ", basename(out_file_go))
      } else {
        message("    No significant GO:BP terms found")
      }
    }, error = function(e) {
      message("    GO enrichment failed: ", e$message)
    })
  } else {
    message("  Insufficient genes for GO:BP enrichment (hypermethylated target genes)")
  }
  
  # Perform GO:BP enrichment for hypermethylated enhancer upregulated target genes
  if (length(hyper_upreg_genes) >= 3) {
    message("  Running GO:BP enrichment for hypermethylated enhancer upregulated target genes...")
    tryCatch({
      go_hyper_up <- enrichGO(
        gene = hyper_upreg_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )
      
      if (!is.null(go_hyper_up) && nrow(as.data.frame(go_hyper_up)) > 0) {
        out_file_go <- file.path(output_dir, paste0(comp, "_Hypermethylated_Upregulated_TargetGenes_GO_BP.csv"))
        fwrite(as.data.frame(go_hyper_up), out_file_go)
        message("    Saved to: ", basename(out_file_go))
      } else {
        message("    No significant GO:BP terms found")
      }
    }, error = function(e) {
      message("    GO enrichment failed: ", e$message)
    })
  } else {
    message("  Insufficient genes for GO:BP enrichment (hypermethylated upregulated genes)")
  }
  
  # Perform GO:BP enrichment for hypomethylated enhancer target genes
  if (length(hypo_target_genes) >= 3) {
    message("  Running GO:BP enrichment for hypomethylated enhancer target genes...")
    tryCatch({
      go_hypo_all <- enrichGO(
        gene = hypo_target_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )
      
      if (!is.null(go_hypo_all) && nrow(as.data.frame(go_hypo_all)) > 0) {
        out_file_go <- file.path(output_dir, paste0(comp, "_Hypomethylated_TargetGenes_GO_BP.csv"))
        fwrite(as.data.frame(go_hypo_all), out_file_go)
        message("    Saved to: ", basename(out_file_go))
      } else {
        message("    No significant GO:BP terms found")
      }
    }, error = function(e) {
      message("    GO enrichment failed: ", e$message)
    })
  } else {
    message("  Insufficient genes for GO:BP enrichment (hypomethylated target genes)")
  }
  
  # Perform GO:BP enrichment for hypomethylated enhancer upregulated target genes
  if (length(hypo_upreg_genes) >= 3) {
    message("  Running GO:BP enrichment for hypomethylated enhancer upregulated target genes...")
    tryCatch({
      go_hypo_up <- enrichGO(
        gene = hypo_upreg_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )
      
      if (!is.null(go_hypo_up) && nrow(as.data.frame(go_hypo_up)) > 0) {
        out_file_go <- file.path(output_dir, paste0(comp, "_Hypomethylated_Upregulated_TargetGenes_GO_BP.csv"))
        fwrite(as.data.frame(go_hypo_up), out_file_go)
        message("    Saved to: ", basename(out_file_go))
      } else {
        message("    No significant GO:BP terms found")
      }
    }, error = function(e) {
      message("    GO enrichment failed: ", e$message)
    })
  } else {
    message("  Insufficient genes for GO:BP enrichment (hypomethylated upregulated genes)")
  }

  # ----------------------------------------------------------------------------
  # Gene Annotation with mygene.info
  # ----------------------------------------------------------------------------
  
  # Extract upregulated and downregulated genes from both hyper and hypo enhancers
  hyper_down_genes <- c()
  hypo_down_genes <- c()
  
  if (!is.null(hyper_combined) && "Regulated_TargetGenes" %in% colnames(hyper_combined)) {
    reg_vals <- na.omit(hyper_combined$Regulated_TargetGenes)
    reg_vals <- reg_vals[reg_vals != ""]
    if (length(reg_vals) > 0) {
      all_reg <- unlist(strsplit(reg_vals, ", "))
      down_reg <- all_reg[grepl("\\(Down\\)$", all_reg)]
      hyper_down_genes <- unique(gsub("\\(Down\\)$", "", down_reg))
    }
  }
  
  if (!is.null(hypo_combined) && "Regulated_TargetGenes" %in% colnames(hypo_combined)) {
    reg_vals <- na.omit(hypo_combined$Regulated_TargetGenes)
    reg_vals <- reg_vals[reg_vals != ""]
    if (length(reg_vals) > 0) {
      all_reg <- unlist(strsplit(reg_vals, ", "))
      down_reg <- all_reg[grepl("\\(Down\\)$", all_reg)]
      hypo_down_genes <- unique(gsub("\\(Down\\)$", "", down_reg))
    }
  }
  
  # Annotate upregulated genes from hypermethylated enhancers
  if (length(hyper_upreg_genes) > 0) {
    message("  Annotating upregulated genes from hypermethylated enhancers with mygene.info...")
    tryCatch({
      mg <- mygene::MyGene()
      hyper_up_annot <- mygene::queryMany(hyper_upreg_genes, scopes="symbol", 
                                          fields=c("name", "summary", "entrezgene", "ensembl.gene", "go.BP"),
                                          species="human", return.as="DataFrame")
      hyper_up_annot <- as.data.frame(hyper_up_annot)
      
      # Flatten list columns for CSV export
      for (col in colnames(hyper_up_annot)) {
        if (is.list(hyper_up_annot[[col]])) {
          hyper_up_annot[[col]] <- sapply(hyper_up_annot[[col]], function(x) {
            if (is.null(x) || length(x) == 0) return(NA_character_)
            if (is.list(x)) return(paste(unlist(x), collapse="; "))
            return(paste(x, collapse="; "))
          })
        }
      }
      
      if (!is.null(hyper_up_annot) && nrow(hyper_up_annot) > 0) {
        out_file <- file.path(output_dir, paste0(comp, "_Hypermethylated_Upregulated_Genes_Annotated.csv"))
        fwrite(hyper_up_annot, out_file)
        message("    Saved to: ", basename(out_file))
      }
    }, error = function(e) {
      message("    Gene annotation failed: ", e$message)
    })
  }
  
  # Annotate downregulated genes from hypermethylated enhancers
  if (length(hyper_down_genes) > 0) {
    message("  Annotating downregulated genes from hypermethylated enhancers with mygene.info...")
    tryCatch({
      mg <- mygene::MyGene()
      hyper_down_annot <- mygene::queryMany(hyper_down_genes, scopes="symbol", 
                                            fields=c("name", "summary", "entrezgene", "ensembl.gene", "go.BP"),
                                            species="human", return.as="DataFrame")
      hyper_down_annot <- as.data.frame(hyper_down_annot)
      
      # Flatten list columns for CSV export
      for (col in colnames(hyper_down_annot)) {
        if (is.list(hyper_down_annot[[col]])) {
          hyper_down_annot[[col]] <- sapply(hyper_down_annot[[col]], function(x) {
            if (is.null(x) || length(x) == 0) return(NA_character_)
            if (is.list(x)) return(paste(unlist(x), collapse="; "))
            return(paste(x, collapse="; "))
          })
        }
      }
      
      if (!is.null(hyper_down_annot) && nrow(hyper_down_annot) > 0) {
        out_file <- file.path(output_dir, paste0(comp, "_Hypermethylated_Downregulated_Genes_Annotated.csv"))
        fwrite(hyper_down_annot, out_file)
        message("    Saved to: ", basename(out_file))
      }
    }, error = function(e) {
      message("    Gene annotation failed: ", e$message)
    })
  }
  
  # Annotate upregulated genes from hypomethylated enhancers
  if (length(hypo_upreg_genes) > 0) {
    message("  Annotating upregulated genes from hypomethylated enhancers with mygene.info...")
    tryCatch({
      mg <- mygene::MyGene()
      hypo_up_annot <- mygene::queryMany(hypo_upreg_genes, scopes="symbol", 
                                         fields=c("name", "summary", "entrezgene", "ensembl.gene", "go.BP"),
                                         species="human", return.as="DataFrame")
      hypo_up_annot <- as.data.frame(hypo_up_annot)
      
      # Flatten list columns for CSV export
      for (col in colnames(hypo_up_annot)) {
        if (is.list(hypo_up_annot[[col]])) {
          hypo_up_annot[[col]] <- sapply(hypo_up_annot[[col]], function(x) {
            if (is.null(x) || length(x) == 0) return(NA_character_)
            if (is.list(x)) return(paste(unlist(x), collapse="; "))
            return(paste(x, collapse="; "))
          })
        }
      }
      
      if (!is.null(hypo_up_annot) && nrow(hypo_up_annot) > 0) {
        out_file <- file.path(output_dir, paste0(comp, "_Hypomethylated_Upregulated_Genes_Annotated.csv"))
        fwrite(hypo_up_annot, out_file)
        message("    Saved to: ", basename(out_file))
      }
    }, error = function(e) {
      message("    Gene annotation failed: ", e$message)
    })
  }
  
  # Annotate downregulated genes from hypomethylated enhancers
  if (length(hypo_down_genes) > 0) {
    message("  Annotating downregulated genes from hypomethylated enhancers with mygene.info...")
    tryCatch({
      mg <- mygene::MyGene()
      hypo_down_annot <- mygene::queryMany(hypo_down_genes, scopes="symbol", 
                                           fields=c("name", "summary", "entrezgene", "ensembl.gene", "go.BP"),
                                           species="human", return.as="DataFrame")
      hypo_down_annot <- as.data.frame(hypo_down_annot)
      
      # Flatten list columns for CSV export
      for (col in colnames(hypo_down_annot)) {
        if (is.list(hypo_down_annot[[col]])) {
          hypo_down_annot[[col]] <- sapply(hypo_down_annot[[col]], function(x) {
            if (is.null(x) || length(x) == 0) return(NA_character_)
            if (is.list(x)) return(paste(unlist(x), collapse="; "))
            return(paste(x, collapse="; "))
          })
        }
      }
      
      if (!is.null(hypo_down_annot) && nrow(hypo_down_annot) > 0) {
        out_file <- file.path(output_dir, paste0(comp, "_Hypomethylated_Downregulated_Genes_Annotated.csv"))
        fwrite(hypo_down_annot, out_file)
        message("    Saved to: ", basename(out_file))
      }
    }, error = function(e) {
      message("    Gene annotation failed: ", e$message)
    })
  }

  # Collect Stats
  summary_stats[[comp]] <- data.frame(
    Comparison = comp,
    Hyper_dmrseq = length(dmrseq_hyper_enhancers),
    Hyper_Targeted = length(targeted_hyper_enhancers),
    Hyper_Union = ifelse(is.null(hyper_combined), 0, nrow(hyper_combined)),
    Hyper_Genes_Up = hyper_counts["Up"],
    Hyper_Genes_Down = hyper_counts["Down"],
    
    Hypo_dmrseq = length(dmrseq_hypo_enhancers),
    Hypo_Targeted = length(targeted_hypo_enhancers),
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
  fwrite(final_stats, file.path(output_dir, "Enhancer_Methylation_Summary_Counts.csv"))
} else {
  message("No stats collected.")
}

message("\nAnalysis Complete. Results saved to: ", output_dir)
