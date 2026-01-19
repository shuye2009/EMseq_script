rm(list = ls())

library(GenomicFeatures)
library(txdbmaker)
library(rtracklayer)
library(readr)
library(dplyr)

wd <- "C:/PROJECTS/resource/T2T_CHM13/"
setwd(wd)

# Create TxDb from GFF3
txdb <- makeTxDbFromGFF("Homo_sapiens-GCA_009914755.4-2022_07-genes.gff3",
                        format = "gff3",
                        organism = "Homo sapiens")

# Also read GFF3 to get Name and biotype fields
gff3_orig <- import("Homo_sapiens-GCA_009914755.4-2022_07-genes.gff3", format = "gff3")
unique(gff3_orig$type)

# Create lookup tables for gene names and biotypes
gene_info <- gff3_orig[gff3_orig$type %in% c("gene", "ncRNA_gene", "pseudogene")]
gene_lookup <- data.frame(
  gene_id = gene_info$gene_id,
  gene_name = gene_info$Name,
  gene_type = gene_info$biotype,
  gene_description = gene_info$description,
  stringsAsFactors = FALSE
)

tx_info <- gff3_orig[gff3_orig$type %in% c("pseudogenic_transcript", "mRNA", "lnc_RNA", "ncRNA", "rRNA", "snRNA", "snoRNA", "miRNA", "scRNA")]
tx_lookup <- data.frame(
  transcript_id = tx_info$transcript_id,
  gene_id = sapply(gsub("gene:", "", tx_info$Parent), function(x)x[1]),
  transcript_type = tx_info$biotype,
  stringsAsFactors = FALSE
)

tx_lookup <- dplyr::left_join(tx_lookup, gene_lookup)

exon_info <- gff3_orig[gff3_orig$type == "exon"]
sum(duplicated(exon_info$exon_id))
head(exon_info[duplicated(exon_info$exon_id),])

exon_lookup <- data.frame(
  exon_id = exon_info$exon_id,
  exon_rank = exon_info$rank,
  transcript_id = sapply(gsub("transcript:", "", exon_info$Parent), function(x)x[1])
)

exon_lookup <- dplyr::left_join(exon_lookup, tx_lookup)

# add gene_name to gene_info
genes <- gene_info
genes$gene_name <- gene_lookup$gene_name

# add gene_name to tx_info
transcripts <- tx_info
transcripts$gene_name <- tx_lookup$gene_name
transcripts$gene_id <- tx_lookup$gene_id

# find transcripts without gene_id
na_transcripts <- transcripts[is.na(transcripts$gene_id),]
length(na_transcripts$transcript_id)
transcripts <- transcripts[!is.na(transcripts$gene_id),]

# add gene_name to tx_info
exons <- exon_info
exons$gene_name <- exon_lookup$gene_name
exons$gene_id <- exon_lookup$gene_id
exons$transcript_id <- exon_lookup$transcript_id

# find exons without gene_id, to make sure each exon has a gene_id which is required by RSEM index
na_exons <- exons[is.na(exons$gene_id),]
length(na_exons$exon_id)
exons <- exons[!is.na(exons$gene_id),]

# Combine all features
all_features <- c(genes, transcripts, exons)

# Convert seqnames to UCSC style
seqlevelsStyle(all_features) <- "UCSC"

# Remove unwanted mcols
mcols(all_features)$ID <- NULL
mcols(all_features)$Name <- NULL
mcols(all_features)$Parent <- NULL
mcols(all_features)$constitutive <- NULL
mcols(all_features)$protein_id <- NULL

# change chromosome style to UCSC style


# Export to GTF
export(all_features, "Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf", 
       format = "gtf")

# convert bed to gtf
bed_repeat <- read_tsv("chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed", col_names = FALSE)

# Modify duplicated entries in X4 column
bed_repeat <- bed_repeat %>%
  group_by(X4) %>%
  mutate(dup_count = row_number() - 1) %>%
  mutate(X4dup = if_else(dup_count > 0, 
                      paste0(X4, "_dup", dup_count), 
                      X4)) %>%
  ungroup() %>%
  select(-dup_count) %>%
  mutate(X4 = X4dup) %>%
  select(-X4dup)
write_tsv(bed_repeat, "chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14_dup.bed", col_names = FALSE)
# Convert to GRanges for GTF export
bed_repeat_gr <- GRanges(
  seqnames = bed_repeat$X1,
  ranges = IRanges(start = bed_repeat$X2 + 1, end = bed_repeat$X3),  # BED is 0-based, GTF is 1-based
  strand = bed_repeat$X6
)

# Add required GTF attributes
bed_repeat_gr$source <- "RepeatMasker"
bed_repeat_gr$type <- "exon"
bed_repeat_gr$biotype <- "repeat"
bed_repeat_gr$gene_id <- bed_repeat$X4
bed_repeat_gr$transcript_id <- bed_repeat$X4dup
bed_repeat_gr$family_id <- bed_repeat$X8
bed_repeat_gr$class_id <- bed_repeat$X7
bed_repeat_gr$gene_name <- paste0(bed_repeat$X4, ":TE")


# Export to GTF
export(bed_repeat_gr, "chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.gtf", format = "gtf")

censat <- read_tsv("chm13v2.0_censat_v2.1.bed", col_names = FALSE)

# Replace ; with : only inside parentheses which causes parsing error for STAR and RSEM
replace_semicolon_in_brackets <- function(x) {
  chars <- strsplit(x, "", fixed = TRUE)[[1]]
  in_bracket <- FALSE
  for (i in seq_along(chars)) {
    if (chars[i] == "(") in_bracket <- TRUE
    if (chars[i] == ")") in_bracket <- FALSE
    if (in_bracket && chars[i] %in% c(";", ",")) chars[i] <- "|"
  }
  paste(chars, collapse = "")
}

censat$X4_clean <- sapply(censat$X4, replace_semicolon_in_brackets)

censat_gr <- GRanges(
  seqnames = censat$X1,
  ranges = IRanges(start = censat$X2 + 1, end = censat$X3),
  strand = "+"
)

# Add required GTF attributes using cleaned names
censat_gr$source <- "T2T_CHM13"
censat_gr$type <- "exon"
censat_gr$biotype <- "censat"
censat_gr$gene_id <- censat$X4_clean
censat_gr$transcript_id <- censat$X4_clean
censat_gr$gene_name <- paste0(censat$X4_clean, ":T2T")

export(censat_gr, "chm13v2.0_censat_v2.1.gtf", format = "gtf")


# concatenate all three gtf files into one
gtf1 <- import("Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf", format = "gtf")
gtf2 <- import("chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.gtf", format = "gtf")
gtf3 <- import("chm13v2.0_censat_v2.1.gtf", format = "gtf")
all_gtf <- c(gtf1, gtf2, gtf3)

# export to gtf
export(all_gtf, "Homo_sapiens-GCA_009914755.4-2022_07-genes_repeat_censat.gtf", format = "gtf")
