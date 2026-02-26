#!/usr/bin/env Rscript

# Function to create chromInfo dataframe from FASTA file
create_chrominfo_from_fasta <- function(fasta_file) {
  # Required packages
  if (!require("Biostrings", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("Biostrings")
    library(Biostrings)
  }
  
  # Check if file exists
  if (!file.exists(fasta_file)) {
    stop("FASTA file not found: ", fasta_file)
  }
  
  # Check file extension and report compression status
  if (grepl("\\.gz$", fasta_file)) {
    message("Reading compressed FASTA file: ", fasta_file)
  } else {
    message("Reading FASTA file: ", fasta_file)
  }
  
  # Read the FASTA file (Biostrings automatically handles .gz files)
  fasta_sequences <- readDNAStringSet(fasta_file)
  
  # Extract chromosome names and lengths
  chrom_names <- names(fasta_sequences)
  chrom_lengths <- width(fasta_sequences)
  
  # Create chromInfo dataframe
  chromInfo <- data.frame(
    chrom = chrom_names,
    length = chrom_lengths,
    stringsAsFactors = FALSE
  )
  
  # Sort by chromosome name (you can modify this sorting as needed)
  # For human chromosomes, you might want to sort numerically for chr1-22, X, Y
  if (all(grepl("^chr", chrom_names))) {
    # Extract chromosome numbers for proper sorting
    chrom_nums <- gsub("^chr", "", chrom_names)
    # Handle X, Y, M differently
    is_numeric <- grepl("^[0-9]+$", chrom_nums)
    numeric_chroms <- chrom_names[is_numeric]
    numeric_nums <- as.numeric(chrom_nums[is_numeric])
    non_numeric_chroms <- chrom_names[!is_numeric]
    
    # Sort numeric chromosomes
    sorted_numeric <- numeric_chroms[order(numeric_nums)]
    # Sort non-numeric chromosomes alphabetically
    sorted_non_numeric <- sort(non_numeric_chroms)
    
    # Combine and reorder chromInfo
    final_order <- c(sorted_numeric, sorted_non_numeric)
    chromInfo <- chromInfo[match(final_order, chromInfo$chrom), ]
  }
  
  message("Created chromInfo with ", nrow(chromInfo), " chromosomes")
  return(chromInfo)
}

# Example usage
if (interactive() || !exists("sourced")) {
  # Example with uncompressed FASTA file
  # fasta_file <- "path/to/your/genome.fa"
  # chromInfo <- create_chrominfo_from_fasta(fasta_file)
  # print(chromInfo)
  
  # Example with compressed FASTA file (.fa.gz)
  # fasta_file_gz <- "path/to/your/genome.fa.gz"
  # chromInfo_gz <- create_chrominfo_from_fasta(fasta_file_gz)
  # print(chromInfo_gz)
  
  # Save to file
  # write.table(chromInfo, "chromInfo.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("Function loaded. Use create_chrominfo_from_fasta('your_file.fa' or 'your_file.fa.gz') to create chromInfo.")
}

# Mark as sourced when called from another script
sourced <- TRUE
