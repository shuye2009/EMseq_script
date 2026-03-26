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
cran_packages <- c("devtools")
for (pkg in cran_packages) {
  install_if_missing(pkg)
}

# Install and load Bioconductor packages
bioc_packages <- c("BSgenome", "BSgenomeForge", "Biostrings", "rtracklayer")
for (pkg in bioc_packages) {
  install_bioc_if_missing(pkg)
}


# ==============================================================================
# Create BSgenome Package for hs1 (T2T-CHM13v2.0)
# ==============================================================================

# Paths
genome_fa <- "C:/PROJECTS/resource/T2T_CHM13/chm13v2.0.fa"
work_dir <- "C:/PROJECTS/resource/T2T_CHM13/BSgenome_hs1"
seqs_dir <- file.path(work_dir, "seqs")

# Create directories
if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
if (!dir.exists(seqs_dir)) dir.create(seqs_dir, recursive = TRUE)

message("Step 1: Splitting genome FASTA into individual chromosome files...")

# Read and split genome FASTA by chromosome
library(Biostrings)
genome_seqs <- readDNAStringSet(genome_fa)

# Get chromosome names
chr_names <- names(genome_seqs)
message("Found ", length(chr_names), " sequences in genome")

# Write individual chromosome FASTA files
for (i in seq_along(genome_seqs)) {
  chr_name <- chr_names[i]
  # Clean up chromosome name (remove everything after first space)
  chr_id <- sub(" .*", "", chr_name)
  
  message("  Writing ", chr_id, ".fa")
  writeXStringSet(genome_seqs[i], file.path(seqs_dir, paste0(chr_id, ".fa")))
}

message("\nStep 2: Creating seed file...")

# Create seed file for BSgenome package
# Format seqnames as R vector syntax: c("chr1", "chr2", ...)
chr_ids <- sub(" .*", "", chr_names)
seqnames_vector <- paste0('c("', paste(chr_ids, collapse='", "'), '")')

seed_content <- c(
  "Package: BSgenome.Hsapiens.UCSC.hs1",
  "Title: Full genome sequences for Homo sapiens (T2T-CHM13v2.0/hs1)",
  "Description: Full genome sequences for Homo sapiens as provided by the Telomere-to-Telomere Consortium (T2T-CHM13v2.0, UCSC version hs1).",
  "Version: 1.0.0",
  "organism: Homo sapiens",
  "common_name: Human",
  "provider: T2T",
  "genome: hs1",
  "release_date: 2022-04",
  "source_url: https://github.com/marbl/CHM13",
  "organism_biocview: Homo_sapiens",
  "BSgenomeObjname: Hsapiens",
  paste0("seqs_srcdir: ", seqs_dir),
  paste0("seqnames: ", seqnames_vector),
  "seqfiles_suffix: .fa"
)

seed_file <- file.path(work_dir, "BSgenome.Hsapiens.UCSC.hs1-seed")
writeLines(seed_content, seed_file)
message("Seed file created: ", seed_file)

message("\nStep 3: Forging BSgenome package...")

# Capture output and errors from forgeBSgenomeDataPkg
tryCatch({
  forgeBSgenomeDataPkg(seed_file, destdir = work_dir, replace = TRUE)
  message("forgeBSgenomeDataPkg completed")
}, error = function(e) {
  message("Error in forgeBSgenomeDataPkg: ", e$message)
  stop(e)
})

# List what was actually created
message("\nChecking work directory contents:")
work_contents <- list.dirs(work_dir, recursive = FALSE)
message("Directories found: ", paste(basename(work_contents), collapse = ", "))

message("\nStep 4: Building and installing package...")
pkg_dir <- file.path(work_dir, "BSgenome.Hsapiens.UCSC.hs1")

if (dir.exists(pkg_dir)) {
  # Build the package
  message("Building package...")
  pkg_tar <- devtools::build(pkg_dir, path = work_dir)
  
  # Install the package
  message("Installing package...")
  install.packages(pkg_tar, repos = NULL, type = "source")
  
  message("\nBSgenome package successfully created and installed!")
  message("Load with: library(BSgenome.Hsapiens.UCSC.hs1)")
} else {
  message("\nExpected package directory not found: ", pkg_dir)
  message("Available directories in work_dir:")
  print(work_contents)
  stop("Package directory not created. Check the directory listing above.")
}
