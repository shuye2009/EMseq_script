library(bsseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(R.utils)
library(preprocessCore)


# Direct liftOver: hg38 -> hg19
message("Setting up liftOver: hg38 -> hg19...")
chain_dir <- "C:/PROJECTS/resource/T2T_CHM13"

# hg38 to hg19 chain from UCSC
chain_gz <- file.path(chain_dir, "hg38ToHg19.over.chain.gz")
chain_file <- file.path(chain_dir, "hg38ToHg19.over.chain")

if (!file.exists(chain_file)) {
    if (!file.exists(chain_gz)) {
        message("Downloading hg38 to hg19 chain file...")
        download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
                      chain_gz, mode = "wb")
    }
    message("Decompressing chain file...")
    gunzip(chain_gz, destname = chain_file, remove = FALSE)
}

# Load chain file
message("Loading chain file...")
chain <- import.chain(chain_file)
message("Loaded hg38->hg19 chain: ", length(chain), " chains")

# Install annotation packages if needed
if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}
if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
}
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Extract MAPLE probe IDs from example data
maple_example <- read.csv("C:/PROJECTS/Shane/Harding_250611/wo_chrY/DSS/Beta_values.csv", row.names = 1)
maple_probe_ids <- rownames(maple_example)
message("MAPLE uses ", length(maple_probe_ids), " probes")

# Get both EPIC and 450K annotations (hg19)
annotation_epic <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annotation_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
message("EPIC annotation: ", nrow(annotation_epic), " probes")
message("450K annotation: ", nrow(annotation_450k), " probes")

# Combine annotations (prioritize EPIC, fill in missing with 450K)
# First try EPIC

epic_found <- maple_probe_ids[maple_probe_ids %in% rownames(annotation_epic)]
fourfivezeroK_found <- maple_probe_ids[maple_probe_ids %in% rownames(annotation_450k)]
fourfivezeroK_only <- setdiff(fourfivezeroK_found, epic_found)

common_cols <- intersect(colnames(annotation_epic), colnames(annotation_450k))

annotation_maple_epic <- annotation_epic[epic_found, common_cols]
annotation_maple_450k <- annotation_450k[fourfivezeroK_only, common_cols]

annotation_maple <- rbind(annotation_maple_epic, annotation_maple_450k)

message("MAPLE probes found in EPIC: ", length(epic_found))
message("MAPLE probes filled from 450K: ", length(fourfivezeroK_only))
message("Total MAPLE probes in annotation: ", nrow(annotation_maple), " probes")

# Remove probes with NA coordinates
valid_coords <- !is.na(annotation_maple$chr) & !is.na(annotation_maple$pos)
annotation_maple_clean <- annotation_maple[valid_coords, ]
message("Probes with valid coordinates: ", nrow(annotation_maple_clean), " probes")
message("Probes removed due to NA coordinates: ", sum(!valid_coords))

# Convert to GRanges for coordinate matching
maple_gr <- GRanges(
  seqnames = annotation_maple_clean$chr,
  ranges = IRanges(start = annotation_maple_clean$pos, width = 1),
  strand = "*"
)
names(maple_gr) <- rownames(annotation_maple_clean)

# load data
message("Loading bsseq object...")
bs <- readRDS("C:/PROJECTS/Shane/Harding_250611/wo_chrY/DSS/bs.RDS")

# get methylation data
message("Extracting methylation data...")
meth <- bsseq::getMeth(bs, type="raw", what="perBase", withDimnames=TRUE)
head(meth)

# get coverage data
message("Extracting coverage data...")
coverage <- bsseq::getCoverage(bs, type="Cov", what="perBase", withDimnames=TRUE)
head(coverage)

pos_data <- GenomicRanges::granges(bs)

# remove bs object to free memory
remove(bs)
gc()

# liftOver: T2T-CHM13 -> hg19
message("Lifting over positions...")
pos_data_hg19_list <- liftOver(pos_data, chain)



# Track which positions successfully lifted (some may fail)
keep_idx <- lengths(pos_data_hg19_list) > 0

# Unlist and keep only successful liftOvers
pos_data_hg19 <- unlist(pos_data_hg19_list)

# Fix seqlevels ordering
seqlevels(pos_data_hg19) <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

# After liftover, check for duplicate hg19 positions
pos_hg19_coords <- paste0(seqnames(pos_data_hg19), ":", start(pos_data_hg19))
dup_in_hg19 <- duplicated(pos_hg19_coords)
message("T2T CpGs that map to duplicate hg19 positions: ", sum(dup_in_hg19))

# See examples
dup_coords <- pos_hg19_coords[duplicated(pos_hg19_coords)]
message("Sample duplicate hg19 coordinates: ", paste(head(dup_coords, 3), collapse = ", "))

pos_df <- data.frame(
          chr = as.character(seqnames(pos_data_hg19)),
          start = start(pos_data_hg19),
          end = end(pos_data_hg19),
          strand = strand(pos_data_hg19),
          stringsAsFactors = FALSE) %>%
          mutate(CpG = paste0(chr, ":", start, "-", end))

head(pos_df)

meth <- meth[keep_idx, ]
coverage <- coverage[keep_idx, ]

rownames(meth) <- pos_df$CpG
rownames(coverage) <- pos_df$CpG

head(meth)
head(coverage)

# filter meth_df by coverage
message("Filtering by coverage...")

# Parameters
#coverage_threshold <- 10
sample_frac <- 0.5

for (coverage_threshold in c(1, 5, 10)) {
  keep <- rowSums(coverage > coverage_threshold) > as.integer(sample_frac * ncol(coverage))
  filtered_meth <- as.data.frame(meth[keep, ])
  filtered_cov <- as.data.frame(coverage[keep, ])


  # Extract methylation data for MAPLE probes, and translate rownames to maple probe id
  # Create GRanges from filtered EM-seq positions
  message("Identifying coveraged-filtered EM-seq positions...")
  filtered_pos <- pos_data_hg19[keep]
  message("Filtered EM-seq CpGs: ", length(filtered_pos))

  # Find overlaps between filtered EM-seq data and MAPLE probes
  message("Finding overlaps between filtered EM-seq data and MAPLE probes...")
  overlaps <- findOverlaps(filtered_pos, maple_gr, type = "equal")
  message("EM-seq CpGs matching MAPLE probes: ", length(overlaps))

  # Extract matching methylation data
  emseq_idx <- queryHits(overlaps)
  maple_idx <- subjectHits(overlaps)

  # Get probe IDs for matched positions
  probe_ids <- names(maple_gr)[maple_idx]

  # Check for duplicates (multiple EM-seq CpGs mapping to same MAPLE probe)
  dup_probes <- duplicated(probe_ids)
  message("Duplicate probe mappings: ", sum(dup_probes))

  # Create data frame with probe IDs and methylation values
  meth_with_probes <- cbind(probe_id = probe_ids, filtered_meth[emseq_idx, ])

  # Average methylation for duplicate probe IDs
  filtered_meth_maple <- meth_with_probes %>%
    group_by(probe_id) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
    as.data.frame()

  # Set probe IDs as rownames
  rownames(filtered_meth_maple) <- filtered_meth_maple$probe_id
  filtered_meth_maple <- filtered_meth_maple[, -1]  # Remove probe_id column
  write.csv(filtered_meth_maple, paste0("C:/PROJECTS/Shane/Harding_250611/wo_chrY/DSS/filtered_meth_maple_cov", coverage_threshold, ".csv"))
  message("MAPLE format beta matrix dimensions: ", nrow(filtered_meth_maple), " x ", ncol(filtered_meth_maple))
  message("Sample probe IDs: ", paste(head(rownames(filtered_meth_maple), 3), collapse = ", "))

  # Missing value imputation using KNN
  message("Checking for missing values...")
  filtered_meth_maple_matrix <- as.matrix(filtered_meth_maple)

  na_count <- sum(is.na(filtered_meth_maple_matrix))
  if (na_count > 0) {
    message("Found ", na_count, " missing values (", 
            round(100 * na_count / length(filtered_meth_maple_matrix), 2), "% of total)")
    
    # Increase R expression limit to handle large clusters
    old_expr <- getOption("expressions")
    options(expressions = 500000)
    
    message("Applying KNN imputation (k=5) with chunking to avoid stack overflow...")
    
    if (!requireNamespace("impute", quietly = TRUE)) {
      BiocManager::install("impute")
    }
    library(impute)
    
    # Process in chunks to avoid node stack overflow for large matrices
    chunk_size <- 5000
    n_probes <- nrow(filtered_meth_maple_matrix)
    
    if (n_probes > chunk_size) {
      message("Processing ", n_probes, " probes in chunks of ", chunk_size)
      imputed_matrix <- matrix(nrow = n_probes, ncol = ncol(filtered_meth_maple_matrix))
      
      for (i in seq(1, n_probes, by = chunk_size)) {
        end_idx <- min(i + chunk_size - 1, n_probes)
        message("  Imputing chunk ", ceiling(i/chunk_size), " (rows ", i, " to ", end_idx, ")...")
        
        chunk <- filtered_meth_maple_matrix[i:end_idx, , drop = FALSE]
        imputed_chunk <- impute.knn(chunk, k = 5, rowmax = 0.7, colmax = 0.8)$data
        imputed_matrix[i:end_idx, ] <- imputed_chunk
      }
      
      rownames(imputed_matrix) <- rownames(filtered_meth_maple_matrix)
      colnames(imputed_matrix) <- colnames(filtered_meth_maple_matrix)
      filtered_meth_maple_matrix <- imputed_matrix
    } else {
      # Small matrix, process all at once
      imputed_result <- impute.knn(filtered_meth_maple_matrix, k = 5, rowmax = 0.7, colmax = 0.8)
      filtered_meth_maple_matrix <- imputed_result$data
    }
    # Restore expression limit
    options(expressions = old_expr)
    
    message("Imputation complete. Remaining NAs: ", sum(is.na(filtered_meth_maple_matrix)))
  } else {
    message("No missing values found. Skipping imputation.")
  }

  # Use quantile normalization instead of ChAMP, which produces too many NA values

  message("Applying quantile normalization ...")
  filtered_meth_maple_normalized <- normalize.quantiles(filtered_meth_maple_matrix, keep.names = TRUE)

  # Convert back to data frame
  filtered_meth_maple_normalized_df <- as.data.frame(filtered_meth_maple_normalized)

  message("Normalization complete. Normalized dimensions: ", nrow(filtered_meth_maple_normalized_df), " x ", ncol(filtered_meth_maple_normalized_df))


  # Round to 3 decimal places to reduce file size
  filtered_meth_maple_normalized_df <- round(filtered_meth_maple_normalized_df, 3)
  
  # Save MAPLE-compatible beta values
  message("Saving MAPLE-compatible beta values...")
  write.csv(filtered_meth_maple_normalized_df, 
            paste0("C:/PROJECTS/Shane/Harding_250611/wo_chrY/DSS/normalized_Beta_values_MAPLE_cov", coverage_threshold, ".csv"), 
            row.names = TRUE)
  # add rownames as the first column of meth
  meth_df <- data.frame(CpG = rownames(filtered_meth_maple_normalized_df), filtered_meth_maple_normalized_df)
  head(meth_df)
  dim(meth_df)
  write.table(meth_df, 
              paste0("C:/PROJECTS/Shane/Harding_250611/wo_chrY/DSS/normalized_methylation_data_MAPLE_cov", coverage_threshold, ".tsv"), 
              row.names = FALSE, 
              sep = "\t", 
              quote = FALSE)

  # Also save the sample metadata

  meta_data <- data.frame(sample_id = colnames(filtered_meth_maple_normalized_df))
  rownames(meta_data) <- colnames(filtered_meth_maple_normalized_df)
  write.csv(meta_data, "C:/PROJECTS/Shane/Harding_250611/wo_chrY/DSS/meta_data_MAPLE.csv", row.names = TRUE)

}