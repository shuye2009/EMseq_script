# Use fgsea to perform GSEA to identify differentially methylated promoters.
# First step is to create a ranked list of promoters based on difference in methylation, 
# 1. Load differential methylation results from DML table of DSS analysis
# 2. Filter DMLs by intersecting with promoter regions
# 3. Create ranked list of promoters based on difference in methylation using the stat column
# 4. Create a gmt file to associate promoters with DMLs
# Next step is to perform GSEA using fgsea to identify hypomethylated and hypermethylated Promoter sets
# Last step is correlate GSEA results with gene expression

rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(openxlsx)
library(fgsea)
library(GenomicRanges)
library(IRanges)

# Load promoter bed file
promoter_bed <- read.table("C:/PROJECTS/resource/T2T_CHM13/T2T_CHM13_promoters.bed", header = FALSE, sep = "\t")
colnames(promoter_bed) <- c("chr", "start", "end", "name", "score", "strand")

# Function to filter DMLs by promoter regions
filter_dml_by_promoter <- function(dml_df, promoter_bed_df) {
    if (nrow(dml_df) == 0) {
        return(dml_df)
    }
    if (!all(c("chr", "start", "end") %in% colnames(promoter_bed_df))) {
        stop("promoter_bed must contain columns: chr, start, end")
    }
    promoter_gr <- GRanges(
        seqnames = promoter_bed_df$chr,
        ranges = IRanges(start = promoter_bed_df$start + 1, end = promoter_bed_df$end),
        strand = promoter_bed_df$strand
    )
    mcols(promoter_gr)$name <- promoter_bed_df$name

    if (!"chr" %in% colnames(dml_df)) {
        stop("DML table must contain a 'chr' column")
    }

    if (all(c("start", "end") %in% colnames(dml_df))) {
        dml_gr <- GRanges(
            seqnames = dml_df$chr,
            ranges = IRanges(start = dml_df$start, end = dml_df$end),
            strand = "*"
        )
    } else if ("pos" %in% colnames(dml_df)) {
        dml_gr <- GRanges(
            seqnames = dml_df$chr,
            ranges = IRanges(start = dml_df$pos, end = dml_df$pos),
            strand = "*"
        )
    } else if ("position" %in% colnames(dml_df)) {
        dml_gr <- GRanges(
            seqnames = dml_df$chr,
            ranges = IRanges(start = dml_df$position, end = dml_df$position),
            strand = "*"
        )
    } else {
        stop("DML table must contain either (start,end) or pos/position columns")
    }

    hits <- findOverlaps(dml_gr, promoter_gr, ignore.strand = TRUE)
    if (length(hits) == 0) {
        return(dml_df[0, , drop = FALSE])
    }

    kept_rows <- S4Vectors::queryHits(hits)
    promoter_name <- mcols(promoter_gr)$name[S4Vectors::subjectHits(hits)]
    out <- dml_df[kept_rows, , drop = FALSE]
    out$promoter_name <- promoter_name
    out
}

# Function to create ranked list of promoters based on difference in methylation
create_ranked_list <- function(dml_in_promoter_df) {
    stat_df <- dml_in_promoter_df %>%
        mutate(cpg = paste0(chr, "_", pos)) %>%
        distinct(cpg, .keep_all = TRUE) %>%
        arrange(stat)
    ranked_list <- stat_df$stat
    names(ranked_list) <- stat_df$cpg
    return(ranked_list)
}

# Function to make gmt file that associates promoter name with cpg site
make_gmt_file <- function(dml_in_promoter_df, output_file) {
    gmt_df <- dml_in_promoter_df %>%
        mutate(cpg = paste0(chr, "_", pos)) %>%
        select(promoter_name, cpg) %>%
        group_by(promoter_name) %>%
        summarise(cpgs = list(cpg), .groups = "drop")

    con <- file(output_file, "w")
    for (i in seq_len(nrow(gmt_df))) {
        line <- paste(c(gmt_df$promoter_name[i], "NA", unlist(gmt_df$cpgs[i])), collapse = "\t")
        writeLines(line, con)
    }
    close(con)
    invisible(gmt_df)
}

expression_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/RNAseq/deseq2/single_factor_analysis_gene_level"

emseq_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/DSS"

sample_dirs <- list.dirs(emseq_dir, full.names = FALSE, recursive = FALSE)

names(sample_dirs) <- sample_dirs

# redefine order of samples
sample_names <- c("IR2Gy24h_vs_NIR", "IR10Gy24h_vs_NIR", "IR10Gy24h_vs_IR2Gy24h", "IR2Gy6d_vs_NIR", "IR10Gy6d_vs_NIR", "IR10Gy6d_vs_IR2Gy6d", "IR2Gy6d_vs_IR2Gy24h", "IR10Gy6d_vs_IR10Gy24h")

sample_dirs <- sample_dirs[sample_names] 

for (sample_name in sample_names) {
    cat("Processing", sample_name, "\n")
    output_dir <- file.path(emseq_dir, sample_name, "GSEA")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    dml_file <- file.path(emseq_dir, sample_name, paste0("DML_", sample_name, ".tsv"))
    dml <- read.table(dml_file, header = TRUE, sep = "\t")
    dml_promoter <- filter_dml_by_promoter(dml, promoter_bed)
    ranked_list <- create_ranked_list(dml_promoter)
    make_gmt_file(dml_promoter, file.path(output_dir, "DML_in_promoter.gmt"))
    out_file <- file.path(output_dir, paste0("DML_in_promoter_", sample_name, ".tsv"))
    write.table(dml_promoter, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    # run GSEA
    gmt_list <- gmtPathways(file.path(output_dir, "DML_in_promoter.gmt"))
    max_size <- max(lengths(gmt_list))
    cat("Max set size:", max_size, "\n")
    gsea_result <- fgsea(pathways = gmt_list, stats = ranked_list, minSize = 5, maxSize = max_size)
    
    # sort by NES then padj
    gsea_result <- gsea_result %>%
        arrange(desc(NES), padj)
    
    # save GSEA result
    gsea_result$leadingEdge <- sapply(gsea_result$leadingEdge, paste, collapse = ",")
    write.table(gsea_result, file.path(output_dir, paste0("GSEA_result_", sample_name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    
    # plot top significant pathways for positive and negative NES
    n_top <- 3
    top_pos <- gsea_result %>% filter(NES > 0) %>% arrange(padj) %>% head(n_top)
    top_neg <- gsea_result %>% filter(NES < 0) %>% arrange(padj) %>% head(n_top)
    top_pathways <- bind_rows(top_pos, top_neg) %>%
        mutate(pathway = factor(pathway, levels = pathway[order(NES)]))
    
    if (nrow(top_pathways) > 0) {
        p <- ggplot(top_pathways, aes(x = NES, y = pathway, fill = padj)) +
            geom_col() +
            scale_fill_gradient(low = "red", high = "blue", name = "padj") +
            labs(title = paste("Top GSEA Pathways -", sample_name),
                 x = "Normalized Enrichment Score (NES)",
                 y = "Promoter") +
            theme_minimal() +
            theme(axis.text.y = element_text(size = 8))
        ggsave(file.path(output_dir, paste0("GSEA_top_promoters_", sample_name, ".pdf")), p, width = 10, height = 8)
    }
    
    # plotEnrichment for top positive and negative NES promoters
    pdf(file.path(output_dir, paste0("GSEA_enrichment_plots_", sample_name, ".pdf")), width = 8, height = 6)
    for (prom in top_pos$pathway) {
        print(plotEnrichment(gmt_list[[prom]], ranked_list) + 
              labs(title = paste("Hypermethylated:", prom)))
    }
    for (prom in top_neg$pathway) {
        print(plotEnrichment(gmt_list[[prom]], ranked_list) + 
              labs(title = paste("Hypomethylated:", prom)))
    }
    dev.off()
}
