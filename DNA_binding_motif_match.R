# ==============================================================================
# Package Installation and Loading
# ==============================================================================

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

# Function to check and install CRAN packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message("Installing package: ", pkg)
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Install and load CRAN packages
cran_packages <- c("devtools", "data.table", "stringr")
for (pkg in cran_packages) {
  install_if_missing(pkg)
}

# Install and load Bioconductor packages
bioc_packages <- c("JASPAR2022", "TFBSTools", "motifmatchr", "rtracklayer", "SummarizedExperiment")
for (pkg in bioc_packages) {
  install_bioc_if_missing(pkg)
}


# ==============================================================================
# Motif Matching Example
# ==============================================================================

# Load the customr-built hs1 genome
# See build_hs1_BSgenome.R for the details

library(BSgenome.Hsapiens.UCSC.hs1)
print(packageVersion("BSgenome.Hsapiens.UCSC.hs1"))

# Get multiple TF motifs relevant to 10Gy IR at 10 days post-treatment
# Focus on DNA damage response, senescence, SASP, and stress response
motif_names <- c("RELA", "FOS", "JUN", "CEBPB", "STAT3", "TP53", "MYC")

message("\nRetrieving motifs from JASPAR2022...")
pwm_list <- list()
for (motif_name in motif_names) {
  tryCatch({
    # Get motif using getMatrixSet with name filter
    opts <- list(name = motif_name)
    pfm_list_obj <- getMatrixSet(JASPAR2022, opts)
    
    # Check if any motifs were found
    if (length(pfm_list_obj) == 0) {
      message("  WARNING: No motif found for ", motif_name)
      next
    }
    
    # Check how many motifs were found
    n_motifs <- length(pfm_list_obj)
    if (n_motifs > 1) {
      message("  Found: ", motif_name, " (using first of ", n_motifs, " versions)")
    } else {
      message("  Found: ", motif_name)
    }
    
    # Extract first motif from PFMatrixList and convert to PWMatrix
    pfm <- pfm_list_obj[[1]]
    pwm_list[[motif_name]] <- toPWM(pfm)
    
  }, error = function(e) {
    message("  ERROR retrieving ", motif_name, ": ", e$message)
  })
}

# Check if we got any motifs
if (length(pwm_list) == 0) {
  stop("No motifs were successfully retrieved. Check JASPAR2022 database connection.")
}

message("\nSuccessfully retrieved ", length(pwm_list), " motifs")

# Convert to PWMatrixList for batch processing
pwm_list <- do.call(PWMatrixList, pwm_list)

message("\nStarting motif matching with hs1 genome...")

# Load hypomethylated enhancers
targeted_base_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/Enhancer_methylation_summary"
enhancer_gr <- import.bed(file.path(targeted_base_dir, "IR10Gy6d_vs_NIR_Hypomethylated_Enhancers_Union.bed"))

enhancer_targetGene <- read.csv(file.path(targeted_base_dir, "IR10Gy6d_vs_NIR_Hypomethylated_Enhancers_Union.csv"))

# Remove duplicates from CSV if any exist (based on name + source)
n_before <- nrow(enhancer_targetGene)
enhancer_targetGene <- enhancer_targetGene[!duplicated(enhancer_targetGene[, c("name", "source")]), ]
n_after <- nrow(enhancer_targetGene)
if (n_before > n_after) {
  message("Removed ", n_before - n_after, " duplicate rows from CSV")
}

# Match all motifs at once
motif_matches <- matchMotifs(pwm_list, enhancer_gr, genome = "hs1", out = "scores")

message("Motif matching complete!")

# ==============================================================================
# Export Motif Match Results to CSV
# ==============================================================================

# Extract scores for all motifs
scores <- assays(motif_matches)$motifScores

# Create data frame with enhancer positions
# Split the name column from BED file (format: "enhancerName|source")
bed_names <- mcols(enhancer_gr)$name
name_parts <- str_split(bed_names, "\\|", n = 2, simplify = TRUE)

results_df <- data.frame(
  chrom = as.character(seqnames(enhancer_gr)),
  start = start(enhancer_gr),
  end = end(enhancer_gr),
  name = name_parts[, 1],  # First part is enhancer name
  source = name_parts[, 2]  # Last part is source
)

# Add motif scores as separate columns
motif_idx <- 1
for (motif_name in motif_names) {
  if (motif_name %in% names(pwm_list)) {
    col_name <- paste0(motif_name, "_score")
    results_df[[col_name]] <- scores[, motif_idx]
    motif_idx <- motif_idx + 1
  }
}

# Merge with target gene information from CSV
# Use name and source columns for matching
results_df <- merge(results_df, enhancer_targetGene, 
                    by = c("name", "source"), 
                    all.x = TRUE, 
                    sort = FALSE,
                    suffixes = c("", ".csv"))

# Remove any remaining duplicates based on enhancer coordinates
# Check for column name suffixes from merge
coord_cols <- if ("chrom" %in% colnames(results_df)) {
  c("chrom", "start", "end")
} else if ("chrom.csv" %in% colnames(results_df)) {
  c("chrom.csv", "start.csv", "end.csv")
} else {
  c("name", "source")  # Fallback to name+source
}

n_before_dedup <- nrow(results_df)
results_df <- results_df[!duplicated(results_df[, coord_cols]), ]
n_after_dedup <- nrow(results_df)
if (n_before_dedup > n_after_dedup) {
  message("Removed ", n_before_dedup - n_after_dedup, " duplicate rows after merge")
}

# Calculate max score across all motifs for sorting
score_cols <- grep("_score$", colnames(results_df), value = TRUE)
results_df$max_motif_score <- apply(results_df[, score_cols], 1, max, na.rm = TRUE)

# Identify which motif has the maximum score for each enhancer
# Only show motif name if max score > 0
results_df$max_motif_name <- apply(results_df[, score_cols], 1, function(x) {
  if (all(is.na(x))) return(NA_character_)
  max_score <- max(x, na.rm = TRUE)
  if (max_score <= 0) return(NA_character_)
  max_idx <- which.max(x)
  motif_name <- sub("_score$", "", score_cols[max_idx])
  return(motif_name)
})

# Sort by max motif score (highest first)
results_df <- results_df[order(-results_df$max_motif_score), ]

# Export to CSV
output_file <- file.path(targeted_base_dir, "IR10Gy6d_vs_NIR_Hypomethylated_Enhancers_Multi_Motif_Scores.csv")
fwrite(results_df, output_file)

message("\n=== Motif Match Results ===")
message("Output file: ", basename(output_file))
message("Total enhancers analyzed: ", nrow(results_df))
message("\nMotif-specific statistics:")
for (motif_name in motif_names) {
  col_name <- paste0(motif_name, "_score")
  if (col_name %in% colnames(results_df)) {
    n_positive <- sum(results_df[[col_name]] > 0, na.rm = TRUE)
    mean_score <- mean(results_df[[col_name]], na.rm = TRUE)
    max_score <- max(results_df[[col_name]], na.rm = TRUE)
    message(sprintf("  %s: %d positive (%.1f%%), mean=%.2f, max=%.2f", 
                    motif_name, n_positive, 100*n_positive/nrow(results_df), 
                    mean_score, max_score))
  }
}


library(clusterProfiler)
library(msigdbr)
library(dplyr)

# Get senescence gene sets
senescence_sets <- msigdbr(species = "Homo sapiens", 
                           collection = "C5", 
                           subcollection = "GO:BP") %>%
  filter(grepl("SENESCENCE", gs_name))

# Convert to GMT format
# GMT format: gene_set_name\tdescription\tgene1\tgene2\tgene3...
gmt_lines <- senescence_sets %>%
  group_by(gs_name) %>%
  summarize(
    genes = paste(gene_symbol, collapse = "\t"),
    .groups = "drop"
  ) %>%
  mutate(gmt_line = paste(gs_name, gs_name, genes, sep = "\t")) %>%
  dplyr::pull(gmt_line)

# Write to GMT file
gmt_output <- file.path(targeted_base_dir, "senescence_gene_sets.gmt")
writeLines(gmt_lines, gmt_output)

message("\n=== Senescence Gene Sets ===")
message("GMT file: ", basename(gmt_output))
message("Number of gene sets: ", length(gmt_lines))
message("Gene sets included:")
for (set_name in unique(senescence_sets$gs_name)) {
  n_genes <- sum(senescence_sets$gs_name == set_name)
  message(sprintf("  %s: %d genes", set_name, n_genes))
}
