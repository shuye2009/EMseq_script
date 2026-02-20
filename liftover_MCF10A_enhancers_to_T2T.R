# LiftOver MCF10A enhancers from hg38 to T2T-CHM13

rm(list = ls())

library(rtracklayer)
library(GenomicRanges)
library(R.utils)
library(data.table)

# Direct liftOver: hg38 -> T2T-CHM13 (Hs1)
message("Setting up liftOver: hg38 -> T2T-CHM13...")
output_dir <- "C:/PROJECTS/resource/T2T_CHM13"

# hg38 to T2T chain from UCSC
chain_gz <- file.path(output_dir, "hg38ToHs1.over.chain.gz")
chain_file <- file.path(output_dir, "hg38ToHs1.over.chain")

if (!file.exists(chain_file)) {
    if (!file.exists(chain_gz)) {
        message("Downloading hg38 to T2T-CHM13 chain file...")
        download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz",
                      chain_gz, mode = "wb")
    }
    message("Decompressing chain file...")
    gunzip(chain_gz, destname = chain_file, remove = FALSE)
}

# Load chain file
message("Loading chain file...")
chain <- import.chain(chain_file)
message("Loaded hg38->T2T chain: ", length(chain), " chains")

#######################################################################################
# LiftOver enhancer-gene associations file
#######################################################################################

assoc_input <- "C:/PROJECTS/resource/ENCODE_enhancers/ENCFF912FUA_MCF10A_element_gene_links_thresholded_Engreitz.bed"
assoc_output <- file.path(output_dir, "ENCFF912FUA_MCF10A_element_gene_links_thresholded_Engreitz_T2T.bed")

if (file.exists(assoc_input)) {
    message("\n=== LiftOver Enhancer-Gene Associations ===")
    message("Loading enhancer-gene associations...")
    
    assoc_dt <- fread(assoc_input)
    message("Loaded ", nrow(assoc_dt), " associations")
    
    # Create GRanges from enhancer coordinates
    assoc_gr <- GRanges(
        seqnames = assoc_dt$`#chr`,
        ranges = IRanges(start = assoc_dt$start, end = assoc_dt$end)
    )
    
    # Perform liftOver
    message("Performing liftOver on enhancer coordinates...")
    assoc_lifted <- liftOver(assoc_gr, chain)
    
    # Check mapping success
    assoc_mapped_idx <- lengths(assoc_lifted) > 0
    message("Successfully mapped: ", sum(assoc_mapped_idx), " / ", nrow(assoc_dt), " associations")
    
    
    # Handle multiple mappings - keep first
    if (any(lengths(assoc_lifted[assoc_mapped_idx]) > 1)) {
        message("Note: Some regions mapped to multiple locations. Keeping first mapping only.")
    }
    
    # Filter to mapped rows and update enhancer coordinates
    assoc_t2t <- assoc_dt[assoc_mapped_idx]
    
    assoc_lifted_gr <- unlist(assoc_lifted[assoc_mapped_idx])
    row_mapping <- rep(seq_along(which(assoc_mapped_idx)), lengths(assoc_lifted[assoc_mapped_idx]))
    first_occurrence <- !duplicated(row_mapping)
    assoc_lifted_first <- assoc_lifted_gr[first_occurrence]
    
    assoc_t2t$`#chr` <- as.character(seqnames(assoc_lifted_first))
    assoc_t2t$start <- start(assoc_lifted_first)
    assoc_t2t$end <- end(assoc_lifted_first)
    assoc_t2t$name <- paste0(assoc_t2t$`#chr`, ":", assoc_t2t$start, "-", assoc_t2t$end)
    
    # Also liftOver TargetGeneTSS (single position on same chromosome as gene)
    message("Performing liftOver on TargetGeneTSS coordinates...")
    tss_gr <- GRanges(
        seqnames = assoc_dt$`#chr`[assoc_mapped_idx],
        ranges = IRanges(start = assoc_dt$TargetGeneTSS[assoc_mapped_idx],
                         end = assoc_dt$TargetGeneTSS[assoc_mapped_idx])
    )
    
    tss_lifted <- liftOver(tss_gr, chain)
    tss_mapped_idx <- lengths(tss_lifted) > 0
    message("TargetGeneTSS mapped: ", sum(tss_mapped_idx), " / ", length(tss_gr))
    
    # Extract lifted TSS positions
    tss_lifted_gr <- unlist(tss_lifted[tss_mapped_idx])
    tss_row_mapping <- rep(seq_along(which(tss_mapped_idx)), lengths(tss_lifted[tss_mapped_idx]))
    tss_first <- !duplicated(tss_row_mapping)
    tss_lifted_first <- tss_lifted_gr[tss_first]
    
    assoc_t2t$TargetGeneTSS[tss_mapped_idx] <- start(tss_lifted_first)
    assoc_t2t$TargetGeneTSS[!tss_mapped_idx] <- NA_integer_
    message("TSS positions unmapped: ", sum(!tss_mapped_idx))
    
    # Save full lifted associations
    fwrite(assoc_t2t, assoc_output, sep = "\t")
    message("Saved lifted associations to: ", assoc_output)
    message("Lifted associations: ", nrow(assoc_t2t))
    message("Unique genes: ", length(unique(assoc_t2t$TargetGene)))
    
    # Save as standard 6-column BED file
    bed6_output <- sub("\\.bed$", "_6col.bed", assoc_output)
    bed6 <- data.table(
        chrom = assoc_t2t$`#chr`,
        start = assoc_t2t$start,
        end = assoc_t2t$end,
        name = assoc_t2t$name,
        score = 0,
        strand = "*"
    )
    bed6 <- unique(bed6)
    fwrite(bed6, bed6_output, sep = "\t", col.names = FALSE)
    message("Saved 6-column BED to: ", bed6_output)
} else {
    message("Enhancer-gene associations file not found: ", assoc_input)
}

message("\n=== All Done ===")
