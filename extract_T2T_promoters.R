rm(list = ls())

library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

gtf_file <- "C:/PROJECTS/resource/T2T_CHM13/Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf"
txdb <- txdbmaker::makeTxDbFromGFF(gtf_file, format = "gtf", organism = "Homo sapiens")
# Extract promoter coordinates
genes <- genes(txdb)
promoters_gr <- promoters(genes, upstream = 2000, downstream = 500)
# Export promoters with gene names
gtf_data <- rtracklayer::import(gtf_file)
gene_mapping <- as.data.frame(gtf_data) %>%
    dplyr::filter(type == "gene") %>%
    dplyr::select(gene_id, gene_name) %>%
    dplyr::distinct()

# only keep distinct gene names
promoters_df <- as.data.frame(promoters_gr) %>%
    dplyr::left_join(gene_mapping, by = c("gene_id" = "gene_id")) %>%
    dplyr::filter(!is.na(gene_name)) %>%
    dplyr::distinct(gene_name, .keep_all = TRUE)

# Save as BED for genome browser
promoters_bed <- promoters_df %>%
    dplyr::mutate(score = 0) %>%
    dplyr::mutate(name = paste0(gene_name, "_", gene_id)) %>%
    dplyr::select(seqnames, start, end, name, score, strand)
    
write.table(promoters_bed, "C:/PROJECTS/resource/T2T_CHM13/T2T_CHM13_promoters.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)