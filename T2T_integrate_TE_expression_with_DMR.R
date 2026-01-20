rm(list = ls())

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(openxlsx)
library(ggplot2)
library(cowplot)
library(GenomicRanges)
library(plyranges)
library(GenomicFeatures)
library(ChIPseeker)
# Create TxDb from GTF
txdb <- makeTxDbFromGFF("C:/PROJECTS/resource/T2T_CHM13/Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf",
                        format = "gtf",
                        organism = "Homo sapiens")

TE_gtf <- "C:/PROJECTS/resource/T2T_CHM13/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.gtf"

TE_gr <- plyranges::read_gff(TE_gtf)
TE_df <- data.frame(gene_id = TE_gr$gene_id, 
                    gene_name = TE_gr$gene_name, 
                    transcript_id = TE_gr$transcript_id, 
                    family_id = TE_gr$family_id,
                    class_id = TE_gr$class_id)
n_family <- length(unique(TE_df$family_id))
n_class <- length(unique(TE_df$class_id))
n_gene <- length(unique(TE_df$gene_id))
n_gene_name <- length(unique(TE_df$gene_name))
n_transcript <- length(unique(TE_df$transcript_id))
Categories <- c("Class", "Family", "Gene", "Transcript")
Numbers <- c(n_class, n_family, n_gene, n_transcript)

TE_stats <- data.frame(Categories, Numbers)
print(TE_stats)

# get TE transcript count at gene level and keep class_id and family_id
TE_transcript_count_gene_level <- TE_df %>%
    group_by(class_id, family_id, gene_id) %>%
    summarise(transcript_count = n_distinct(transcript_id), .groups = "drop")

TE_gene_count_family_level <- TE_transcript_count_gene_level %>%
    group_by(class_id, family_id) %>%
    summarise(gene_count = n_distinct(gene_id), .groups = "drop")

TE_family_count_class_level <- TE_gene_count_family_level %>%
    group_by(class_id) %>%
    summarise(family_count = n_distinct(family_id), .groups = "drop")


# write TE stats
base_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13"
write_tsv(TE_stats, file.path(base_dir, "TE_stats.tsv"))
write_tsv(TE_transcript_count_gene_level, file.path(base_dir, "TE_transcript_count_gene_level.tsv"))
write_tsv(TE_gene_count_family_level, file.path(base_dir, "TE_gene_count_family_level.tsv"))
write_tsv(TE_family_count_class_level, file.path(base_dir, "TE_family_count_class_level.tsv"))

# count number of transcripts per TE family
TE_transcript_count <- TE_df %>%
    group_by(family_id, class_id) %>%
    summarise(transcript_count = n_distinct(transcript_id), .groups = "drop") %>%
    arrange(desc(transcript_count)) %>%
    group_by(family_id) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    top_n(30, transcript_count)

TE_class_transcript_count <- TE_df %>%
    group_by(class_id) %>%
    summarise(transcript_count = n_distinct(transcript_id), .groups = "drop") %>%
    arrange(desc(transcript_count))

# write TE class transcript count
write_tsv(TE_class_transcript_count, file.path(base_dir, "TE_class_transcript_count.tsv"))

# Define consistent color palette for TE classes based on actual class_id values
unique_classes <- unique(TE_class_transcript_count$class_id)
n_classes <- length(unique_classes)

# Predefined colors for common TE classes
predefined_colors <- c(
    "LINE" = "#E74C3C",        # Red
    "SINE" = "#3498DB",        # Blue
    "LTR" = "#2ECC71",         # Green
    "DNA" = "#F39C12",         # Orange
    "Retroposon" = "#9B59B6",  # Purple
    "RC" = "#1ABC9C",          # Teal
    "DNA?" = "#E67E22",        # Dark Orange
    "LTR?" = "#16A085",        # Dark Teal
    "Satellite" = "#e751bd",   # Magenta
    "RNA" = "#8E44AD",         # Dark Purple
    "SINE?" = "#27AE60",       # Dark Green
    "RC?" = "#2980B9",         # Dark Blue
    "Unknown" = "#95A5A6",     # Gray
    "Other" = "#34495E"        # Dark Gray
)

# Build base_colors with only classes present in data
base_colors <- c()
for(class_name in unique_classes) {
    if(class_name %in% names(predefined_colors)) {
        base_colors[class_name] <- predefined_colors[class_name]
    }
}

# For any classes not in predefined_colors, generate additional colors
missing_classes <- setdiff(unique_classes, names(base_colors))
if(length(missing_classes) > 0) {
    n_missing <- length(missing_classes)
    if(n_missing <= 8) {
        additional_colors <- RColorBrewer::brewer.pal(max(3, n_missing), "Set2")[1:n_missing]
    } else if(n_missing <= 16) {
        additional_colors <- c(RColorBrewer::brewer.pal(8, "Set2"), 
                               RColorBrewer::brewer.pal(min(8, n_missing - 8), "Set3"))
    } else {
        additional_colors <- c(RColorBrewer::brewer.pal(8, "Set2"),
                               RColorBrewer::brewer.pal(8, "Set3"),
                               rainbow(n_missing - 16))
    }
    names(additional_colors) <- missing_classes
    base_colors <- c(base_colors, additional_colors)
}

cat("TE classes in data (", n_classes, "):", paste(unique_classes, collapse=", "), "\n")

# Use base_colors as te_class_colors
te_class_colors <- base_colors

transcript_plot <- ggplot(TE_transcript_count, aes(x = reorder(family_id, -transcript_count), y = transcript_count, fill = class_id)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = transcript_count), vjust = -0.5, size = 2) +
    scale_fill_manual(values = te_class_colors, na.value = "#BDC3C7") +
    labs(title = "Number of Transcripts per Repeat Family (Top 30)", 
         x = "Repeat Family", 
         y = "Number of Transcripts", 
         fill = "Repeat Class") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, max(TE_transcript_count$transcript_count) * 1.1)
ggsave(transcript_plot, file = file.path(base_dir, "TE_transcript_count.pdf"), width = 10, height = 8)
ggsave(transcript_plot, file = file.path(base_dir, "TE_transcript_count.png"), width = 10, height = 8, dpi = 300)

class_plot <- ggplot(TE_class_transcript_count, aes(x = reorder(class_id, -transcript_count), y = transcript_count, fill = class_id)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = transcript_count), vjust = -0.5, size = 3) +
    scale_fill_manual(values = te_class_colors, na.value = "#BDC3C7") +
    labs(title = "Number of Transcripts per Repeat Class", 
         x = "Repeat Class", 
         y = "Number of Transcripts", 
         fill = "Repeat Class") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, max(TE_class_transcript_count$transcript_count) * 1.1)

ggsave(class_plot, file = file.path(base_dir, "TE_class_transcript_count.pdf"), width = 10, height = 8)
ggsave(class_plot, file = file.path(base_dir, "TE_class_transcript_count.png"), width = 10, height = 8, dpi = 300)

# get TE annotation
annotation_dir <- "C:/PROJECTS/resource/T2T_CHM13/"
te_annotation <- read_tsv(file.path(annotation_dir, "T2T_CHM13_rmsk_TE_annotated_ChIPseeker.tsv"))
head(te_annotation)

# modify TE annotation to remove contents in () in the annotation column
te_annotation$annotation <- gsub(" \\(.+\\)", "", te_annotation$annotation)
head(te_annotation)

# left join the differential TE with the TE annotation
annotated_TE_df <- dplyr::left_join(TE_df, te_annotation, by = c("transcript_id" = "V4"))
head(annotated_TE_df)
unique(annotated_TE_df$annotation)
# count number of transcripts per TE family - get top 20 in each annotation group
TE_transcript_count_annotated <- annotated_TE_df %>%
    filter(!is.na(annotation)) %>%
    group_by(family_id, class_id, annotation) %>%
    summarise(transcript_count = n_distinct(transcript_id), .groups = "drop") %>%
    group_by(annotation) %>%
    arrange(desc(transcript_count), .by_group = TRUE) %>%
    slice_head(n = 20) %>%
    ungroup()
transcript_plot_annotated <- ggplot(TE_transcript_count_annotated, aes(x = reorder(family_id, -transcript_count), y = transcript_count, fill = class_id)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = te_class_colors, na.value = "#BDC3C7") +
    labs(title = "Number of Transcripts per Repeat Family (Top 20)", 
         x = "Repeat Family", 
         y = "Number of Transcripts", 
         fill = "Repeat Class") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 18)) +
    facet_wrap(~ annotation, scale = "free")
ggsave(transcript_plot_annotated, file = file.path(base_dir, "TE_transcript_count_annotated.pdf"), width = 20, height = 16)
ggsave(transcript_plot_annotated, file = file.path(base_dir, "TE_transcript_count_annotated.png"), width = 20, height = 16, dpi = 300)

## process samples
expression_dir <- "C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/RNAseq/deseq2/single_factor_analysis_transcript_level"
# Read BED file and create GRanges object
centromere_bed <- read.table("C:/PROJECTS/resource/T2T_CHM13/chm13v2.0_centromere.bed", 
                             header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(centromere_bed)[1:3] <- c("chr", "start", "end")
if(ncol(centromere_bed) >= 4) colnames(centromere_bed)[4] <- "name"
centromere_gr <- makeGRangesFromDataFrame(centromere_bed, keep.extra.columns = TRUE)

dmr_method <- "dmrseq" #"dmrseq", "TE_targeted" "DSS"

emseq_dir <- file.path("C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq", dmr_method)

dmr_folder <- "DMR"
sample_dirs <- list.dirs(emseq_dir, full.names = FALSE, recursive = FALSE)

names(sample_dirs) <- sample_dirs

# redefine order of samples
sample_names <- c("IR2Gy24h_vs_NIR", "IR10Gy24h_vs_NIR", "IR10Gy24h_vs_IR2Gy24h", "IR2Gy6d_vs_NIR", "IR10Gy6d_vs_NIR", "IR10Gy6d_vs_IR2Gy6d", "IR2Gy6d_vs_IR2Gy24h", "IR10Gy6d_vs_IR10Gy24h")

sample_dirs <- sample_dirs[sample_names] 

TE_expression_list <- list() # TE expression regardless of DMR
TE_expression_dmr_plot_list <- list() # all regulation in DMR
TE_expression_in_dmr_plot_list <- list() # only up-regulated and down-regulated in DMR
TE_transcripts_dmr_list <- list() # number of transcripts per TE family overlapping with DMR
TE_dmr_overlap_list <- list() # number of DMRs overlapping with TE family
dmr_annotation_count_list <- list() # number of DMRs in genomic regions

for(sample_name in sample_names){
    message("Processing ", sample_name, " ##############################################################################################")
    # remove leading "IR" from sample name
    SA <- gsub("^IR", "", sample_name)
    SA <- gsub("_IR", "_", SA)
    
    working_dir <- file.path(emseq_dir, sample_dirs[sample_name], dmr_folder)
    output_dir <- file.path(working_dir, "Integrate_TE_expression_with_DMR")
    if(!dir.exists(output_dir)) dir.create(output_dir)
    # import differential transcript expression results
    dte_TE <- read.csv(paste0(expression_dir, "/", SA, "/", SA, "_DESeq2_results_TE.csv")) %>%
        dplyr::left_join(annotated_TE_df)
    head(dte_TE[!is.na(dte_TE$family_id), ])
    
    dmr <- read_tsv(file.path(working_dir, "DMR_table.tsv")) %>%
        mutate(direction = ifelse(stat > 0, "Hypermethylated", "Hypomethylated"))

    # annotate DMRs
    dmr_gr <- makeGRangesFromDataFrame(dmr, keep.extra.columns = TRUE)
    dmrAnno <- annotatePeak(dmr_gr, 
                         tssRegion=c(-1000, 0),
                         TxDb=txdb, 
                         annoDb="org.Hs.eg.db")

    # View annotation summary
    print(dmrAnno)
    annotated_dmr <- as.data.frame(dmrAnno) %>%
        mutate(annotation = gsub(" \\(.+\\)", "", annotation))

    # count number of DMRs annotateed to Intron etc., grouped by direction and annotation
    dmr_annotation_count <- annotated_dmr %>%
        group_by(direction, annotation) %>%
        summarise(count = dplyr::n(), .groups = "drop") %>%
        mutate(comparison = sample_name)
    
    dmr_annotation_count_list[[sample_name]] <- dmr_annotation_count
    
    # find overlap between DMRs and TE
    overlap <- findOverlaps(dmr_gr, TE_gr, maxgap = 1000)

    # extract each DMR that overlaps with a TE
    overlap_dmr <- dmr_gr[queryHits(overlap)]
    
    # extract each TE that overlaps with a DMR
    overlap_TE <- TE_gr[subjectHits(overlap)]

    head(overlap_dmr)
    head(overlap_TE)

    # convert to data frame
    overlap_dmr_df <- data.frame(overlap_dmr)
    overlap_TE_df <- data.frame(overlap_TE) %>%
        dplyr::left_join(annotated_TE_df)

    head(overlap_dmr_df)
    head(overlap_TE_df)

    # merge DMRs and TE annotations
    colnames(overlap_TE_df) <- paste0("TE_", colnames(overlap_TE_df))
    colnames(overlap_dmr_df) <- paste0("DMR_", colnames(overlap_dmr_df))
    overlap_dmr_TE <- cbind(overlap_dmr_df, overlap_TE_df)
    
    # Calculate actual overlap width (intersection of DMR and TE ranges)
    overlap_dmr_TE$overlap_width <- pmin(overlap_dmr_TE$DMR_end, overlap_dmr_TE$TE_end) - 
                                     pmax(overlap_dmr_TE$DMR_start, overlap_dmr_TE$TE_start) + 1

    head(overlap_dmr_TE)

    # count DMRs that are not overlapping with TE
    overlapping_dmr_indices <- unique(queryHits(overlap))
    non_overlap_dmr <- dmr_gr[-overlapping_dmr_indices]
    non_overlap_dmr_df <- data.frame(non_overlap_dmr) %>%
        group_by(direction) %>%
        summarise(count = dplyr::n(), .groups = "drop") 
    non_overlap_df <- data.frame(TE_family_id = "Non-TE", 
                                DMR_direction = non_overlap_dmr_df$direction, 
                                count = non_overlap_dmr_df$count)
    
    message("Total DMRs: ", length(dmr_gr))
    message("DMRs overlapping with TEs: ", length(overlapping_dmr_indices))
    message("DMRs NOT overlapping with TEs: ", length(non_overlap_dmr))

    # plot overlapping DMRs count by direction and TE family
    overlap_dmr_count <- overlap_dmr_TE %>%
        group_by(TE_family_id, DMR_direction) %>%
        summarise(count = n_distinct(c(DMR_seqnames, DMR_start, DMR_end)), .groups = "drop")


    overlap_dmr_count_hyper <- overlap_dmr_TE %>%
        filter(DMR_direction == "Hypermethylated") %>%
        group_by(TE_family_id) %>%
        summarise(hyper_count = n_distinct(c(DMR_seqnames, DMR_start, DMR_end)), .groups = "drop") 
    
    overlap_dmr_count_hypo <- overlap_dmr_TE %>%
        filter(DMR_direction == "Hypomethylated") %>%
        group_by(TE_family_id) %>%
        summarise(hypo_count = n_distinct(c(DMR_seqnames, DMR_start, DMR_end)), .groups = "drop")
    # full join to get all TE families

    
    # Combine hyper and hypo counts
    overlap_dmr_count <- full_join(
        overlap_dmr_count_hyper,
        overlap_dmr_count_hypo,
        by = "TE_family_id"
    ) %>%
        replace_na(list(hyper_count = 0, hypo_count = 0)) %>%
        arrange(desc(hyper_count + hypo_count)) %>%
        slice_head(n = 10)

    # pivot longer for plotting
    overlap_dmr_long <- overlap_dmr_count %>%
        pivot_longer(cols = c(hyper_count, hypo_count), names_to = "DMR_direction", values_to = "count") %>%
        mutate(DMR_direction = case_when(
            DMR_direction == "hyper_count" ~ "Hypermethylated",
            DMR_direction == "hypo_count" ~ "Hypomethylated",
            TRUE ~ DMR_direction
        )) %>%
        mutate(DMR_direction = factor(DMR_direction, levels = c("Hypermethylated", "Hypomethylated")))

        
    overlap_dmr_count_long <- rbind(overlap_dmr_long, non_overlap_df) %>%
        mutate(comparison = sample_name) %>%
        mutate(TE_family_id = factor(TE_family_id, levels = c("Non-TE", unique(TE_family_id)[!unique(TE_family_id) %in% "Non-TE"])))
    
    overlap_dmr_plot <- ggplot(overlap_dmr_count_long, aes(x = reorder(TE_family_id, -count), y = count, fill = DMR_direction)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = sample_name, 
             x = "Repeat Family", 
             y = "Number of DMRs Overlapping with Repeats") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(overlap_dmr_plot, file = file.path(output_dir, "overlapping_dmr_count.pdf"), width = 10, height = 8)
    ggsave(overlap_dmr_plot, file = file.path(output_dir, "overlapping_dmr_count.png"), width = 10, height = 8, dpi = 300)
    TE_dmr_overlap_list[[sample_name]] <- overlap_dmr_count_long
    
    
    
    # plot overlapping distinct transcripts count by TE family
    overlap_transcripts_count_hyper <- overlap_dmr_TE %>%
        filter(DMR_direction == "Hypermethylated") %>%
        group_by(TE_family_id) %>%
        summarise(hyper_count = n_distinct(TE_transcript_id), .groups = "drop") 
    
    overlap_transcripts_count_hypo <- overlap_dmr_TE %>%
        filter(DMR_direction == "Hypomethylated") %>%
        group_by(TE_family_id) %>%
        summarise(hypo_count = n_distinct(TE_transcript_id), .groups = "drop")
    # full join to get all TE families

    
    # Combine hyper and hypo counts
    overlap_transcripts_count <- full_join(
        overlap_transcripts_count_hyper,
        overlap_transcripts_count_hypo,
        by = "TE_family_id"
    ) %>%
        replace_na(list(hyper_count = 0, hypo_count = 0)) %>%
        arrange(desc(hyper_count + hypo_count)) %>%
        slice_head(n = 10)

    # pivot longer for plotting
    overlap_transcripts_long <- overlap_transcripts_count %>%
        pivot_longer(cols = c(hyper_count, hypo_count), names_to = "DMR_direction", values_to = "count") %>%
        mutate(DMR_direction = case_when(
            DMR_direction == "hyper_count" ~ "Hypermethylated",
            DMR_direction == "hypo_count" ~ "Hypomethylated",
            TRUE ~ DMR_direction
        )) %>%
        mutate(DMR_direction = factor(DMR_direction, levels = c("Hypermethylated", "Hypomethylated"))) %>%
        mutate(comparison = sample_name)
    
    
    overlap_transcripts_plot <- ggplot(overlap_transcripts_long, aes(x = reorder(TE_family_id, -count), y = count, fill = DMR_direction)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = sample_name, 
             x = "Repeat Family", 
             y = "Number of Repeat Transcripts Overlapping with DMRs") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(overlap_transcripts_plot, file = file.path(output_dir, "overlapping_TE_transcripts_count.pdf"), width = 10, height = 8)
    ggsave(overlap_transcripts_plot, file = file.path(output_dir, "overlapping_TE_transcripts_count.png"), width = 10, height = 8, dpi = 300)
    TE_transcripts_dmr_list[[sample_name]] <- overlap_transcripts_long

    # overlap TE with centromere
    overlap_TE_centromere <- findOverlaps(TE_gr, centromere_gr, maxgap = 1000)
    overlap_TEc_df <- data.frame(TE_gr[queryHits(overlap_TE_centromere)])
  
    head(overlap_TEc_df)
    

    # merge DMRs and TE annotations with differential transcript expression results
    TE_dmr_overlap <- dplyr::left_join(dte_TE, overlap_dmr_TE, by = c("transcript_id"="TE_transcript_id", "gene_id"="TE_gene_id", "family_id"="TE_family_id", "class_id"="TE_class_id"))
    TE_dmr_overlap <- TE_dmr_overlap %>%
        dplyr::mutate(inCentromere = ifelse(transcript_id %in% overlap_TEc_df$transcript_id, "yes", "no"))
    # save to excel
    write.xlsx(TE_dmr_overlap, paste0(output_dir, "/TE_annotated_with_DMR.xlsx"))

    #####################################################################################################
    ######### Count TE_transcript_id by family for up and down regulated TEs ############################
    #####################################################################################################
    # Create up and down regulated dataframes
    up_df <- dte_TE %>% filter(log2FoldChange > 0 & padj < 0.05)
    down_df <- dte_TE %>% filter(log2FoldChange < 0 & padj < 0.05)
    

    # simplified 
    up_counts <- up_df %>%
        group_by(family_id) %>%
        summarise(up_count = dplyr::n()) %>%
        ungroup()
    
    down_counts <- down_df %>%
        group_by(family_id) %>%
        summarise(down_count = dplyr::n()) %>%
        ungroup()
    
    # Merge counts
    family_counts <- full_join(up_counts, down_counts, by = c("family_id")) %>%
        replace_na(list(up_count = 0, down_count = 0)) %>%
        arrange(desc(up_count + down_count)) %>% 
        slice_head(n = 10)
    
    # Reshape for plotting
    family_counts_long <- family_counts %>%
        pivot_longer(cols = c(up_count, down_count), 
                    names_to = "regulation", 
                    values_to = "count") %>%
        mutate(regulation = factor(regulation, 
                                   levels = c("up_count", "down_count"),
                                   labels = c("Up-regulated", "Down-regulated"))) %>%
        mutate(comparison = sample_name) %>%
        arrange(desc(count))
    
    # Create bar chart
    p <- ggplot(family_counts_long, aes(x = reorder(family_id, -count), y = count, fill = regulation)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = c("Up-regulated" = "#E74C3C", "Down-regulated" = "#3498DB")) +
        labs(title = sample_name,
             x = "Repeat Family",
             y = "Number of Repeat Transcripts",
             fill = "Regulation") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "top")
    
    # Save plot
    ggsave(paste0(output_dir, "/TE_family_expression_barplot.pdf"), 
           plot = p, width = 12, height = 6)
    ggsave(paste0(output_dir, "/TE_family_expression_barplot.png"), 
           plot = p, width = 12, height = 6, dpi = 300)
    
    # Save counts table
    write.xlsx(family_counts, paste0(output_dir, "/TE_family_expression_counts.xlsx"))
    TE_expression_list[[sample_name]] <- family_counts_long
    message("Created bar chart and saved counts for ", sample_name)

     # Count TE_transcript_id by family for up and down regulated TEs, to be faceted by annotation
    up_counts <- up_df %>%
        group_by(family_id, annotation) %>%
        summarise(up_count = dplyr::n()) %>%
        ungroup()
    
    down_counts <- down_df %>%
        group_by(family_id, annotation) %>%
        summarise(down_count = dplyr::n()) %>%
        ungroup()
    
    # Merge counts and get top 10 families by total count per annotation using slice
    family_counts <- full_join(up_counts, down_counts, by = c("family_id", "annotation")) %>%
        replace_na(list(up_count = 0, down_count = 0)) %>%
        group_by(annotation) %>%
        arrange(desc(up_count + down_count)) %>%
        slice_head(n = 10) %>%
        ungroup()

    # Save the combined counts
    write.xlsx(family_counts, paste0(output_dir, "/TE_family_expression_counts_annotated.xlsx"))
     # Reshape for plotting
    family_counts_long <- family_counts %>%
        pivot_longer(cols = c(up_count, down_count), 
                    names_to = "regulation", 
                    values_to = "count") %>%
        mutate(regulation = factor(regulation, 
                                   levels = c("up_count", "down_count"),
                                   labels = c("Up-regulated", "Down-regulated"))) %>%
        arrange(desc(count)) %>% 
        filter(!is.na(annotation))

    
    p_annotated <- ggplot(family_counts_long, aes(x = reorder(family_id, -count), y = count, fill = regulation)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = c("Up-regulated" = "#E74C3C", "Down-regulated" = "#3498DB")) +
        labs(title = sample_name,
             x = "Repeat Family",
             y = "Number of Repeat Transcripts",
             fill = "Regulation") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "top") +
        facet_wrap(~ annotation, scales = "free")
    
    # Save plot
    ggsave(paste0(output_dir, "/TE_family_expression_annotated_barplot.pdf"), 
           plot = p_annotated, width = 12, height = 6)
    ggsave(paste0(output_dir, "/TE_family_expression_annotated_barplot.png"), 
           plot = p_annotated, width = 12, height = 6, dpi = 300)

    #####################################################################################################
    ######### Count TE_transcript_id by family for up and down regulated TEs that overlap with DMRs #####
    #####################################################################################################
    # Debug: Check data before filtering
    message("Total rows in TE_dmr_overlap: ", nrow(TE_dmr_overlap))
    message("Rows with DMR overlap: ", sum(!is.na(TE_dmr_overlap$DMR_seqnames)))
    message("Rows with padj < 0.05: ", sum(TE_dmr_overlap$padj < 0.05, na.rm = TRUE))
    message("Rows with |log2FC| > 1: ", sum(abs(TE_dmr_overlap$log2FoldChange) > 1, na.rm = TRUE))
    
    # Filter for TEs that overlap with DMRs (non-NA DMR columns indicate overlap)
    TE_dmr_overlap_filtered <- TE_dmr_overlap %>%
        filter(!is.na(DMR_seqnames)) %>%
        mutate(regulation = case_when(
            log2FoldChange > 1  & padj < 0.05 ~ "Up-regulated",
            log2FoldChange < -1 & padj < 0.05 ~ "Down-regulated",
            TRUE ~ "No-change"
        )) 
    
    message("Rows after filtering: ", nrow(TE_dmr_overlap_filtered))
    
    if(nrow(TE_dmr_overlap_filtered) == 0){
        message("No significant TE overlap with DMRs for ", sample_name)
        next
    }
    # Count distinct TE_transcript_id by family, regulation, and methylation status
    dmr_family_counts <- TE_dmr_overlap_filtered %>%
        group_by(family_id, regulation, DMR_direction) %>%
        summarise(count = n_distinct(transcript_id), .groups = "drop") %>%
        mutate(regulation = as.character(regulation),
               DMR_direction = as.character(DMR_direction))
    
    # Create complete combinations to ensure all facets show all families
    all_combinations <- expand.grid(
        family_id = unique(dmr_family_counts$family_id),
        regulation = c("Up-regulated", "Down-regulated", "No-change"),
        DMR_direction = c("Hypermethylated", "Hypomethylated"),
        stringsAsFactors = FALSE
    )
    
    dmr_family_counts_complete <- all_combinations %>%
        left_join(dmr_family_counts, by = c("family_id", "regulation", "DMR_direction")) %>%
        replace_na(list(count = 0)) %>%
        arrange(desc(count))

    # find the top 10 unique families
    top_10_families <- dmr_family_counts_complete %>%
        filter(DMR_direction == "Hypomethylated", regulation == "No-change") %>%
        top_n(10, count) %>%
        pull(family_id)
    
    # Keep only top 10 families for plotting
    dmr_family_counts_complete <- dmr_family_counts_complete %>%
        filter(family_id %in% top_10_families)
    
    
    # Create faceted bar chart
    p_dmr <- ggplot(dmr_family_counts_complete, 
                    aes(x = reorder(family_id, -count), y = count, fill = regulation)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~DMR_direction, ncol = 2) +
        scale_fill_manual(values = c("Up-regulated" = "#E74C3C", "Down-regulated" = "#3498DB", "No-change" = "#95A5A6")) +
        labs(title = sample_name,
             x = "Repeat Family",
             y = "Number of regulated Repeat Transcripts Overlapping with DMRs",
             fill = "Regulation") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "top",
              strip.background = element_rect(fill = "gray90"),
              strip.text = element_text(face = "bold")) +
        ylim(0, 800)
    
    # Save plot
    ggsave(paste0(output_dir, "/TE_family_expression_in_DMR_barplot_faceted.pdf"), 
           plot = p_dmr, width = 14, height = 6)
    ggsave(paste0(output_dir, "/TE_family_expression_in_DMR_barplot_faceted.png"), 
           plot = p_dmr, width = 14, height = 6, dpi = 300)
    
    # filter out the no change
    TE_dmr_overlap_filtered <- TE_dmr_overlap_filtered %>%
        dplyr::filter(regulation != "No-change") %>%
        dplyr::select(c("transcript_id", "seqnames", "start", "end", "gene_id", "gene_name", "family_id", "class_id", "log2FoldChange", "padj", "regulation",  
                        "DMR_direction")) %>%
        dplyr::arrange(regulation, DMR_direction) %>%
        dplyr::distinct()
        
    write.xlsx(TE_dmr_overlap_filtered, paste0(output_dir, "/regulated_TE_family_expression_in_DMR.xlsx"))
    TE_expression_dmr_plot_list[[sample_name]] <- p_dmr
    
    message("Created DMR overlap bar chart and saved counts for ", sample_name)

    # only plot regulated TE expression in DMR
    regulated_dmr <- dmr_family_counts_complete %>%
        filter(regulation != "No-change")
    
    # Create faceted bar chart
    p_dmr <- ggplot(regulated_dmr, 
                    aes(x = reorder(family_id, -count), y = count, fill = regulation)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~DMR_direction, ncol = 2) +
        scale_fill_manual(values = c("Up-regulated" = "#E74C3C", "Down-regulated" = "#3498DB")) +
        labs(title = sample_name,
             x = "Repeat Family",
             y = "Number of regulated Repeat Transcripts Overlapping with DMRs",
             fill = "Regulation") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "top",
              strip.background = element_rect(fill = "gray90"),
              strip.text = element_text(face = "bold")) +
        ylim(0, 50)
    
    # Save plot
    ggsave(paste0(output_dir, "/TE_family_regulated_in_DMR_barplot_faceted.pdf"), 
           plot = p_dmr, width = 14, height = 6)
    ggsave(paste0(output_dir, "/TE_family_regulated_in_DMR_barplot_faceted.png"), 
           plot = p_dmr, width = 14, height = 6, dpi = 300)
    TE_expression_in_dmr_plot_list[[sample_name]] <- p_dmr
    message("Created regulated DMR overlap bar chart and saved counts for ", sample_name)
}

# Save all plots use cowplot as 2 x 4 multiplot
if(length(TE_expression_list) > 0){
    TE_expression_df <- dplyr::bind_rows(TE_expression_list) %>%
        mutate(comparison = factor(comparison, levels = sample_names))
    multiplot_te <- ggplot(TE_expression_df, aes(x = reorder(family_id, -count), y = count, fill = regulation)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~comparison, scales = "free_x") +
        scale_fill_manual(values = c("Up-regulated" = "#E74C3C", "Down-regulated" = "#3498DB")) +
        labs(title = "", x = "", y = "Number of transcripts", fill = "Regulation") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
            axis.text.y = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            plot.title = element_text(size = 16),
            legend.text = element_text(size = 16),
            strip.text = element_text(size = 20)) 

    ggsave(paste0(base_dir, "/TE_expression_counts_barplot_multiplot.pdf"), 
           plot = multiplot_te, width = 20, height = 15)
    ggsave(paste0(base_dir, "/TE_expression_counts_barplot_multiplot.png"), 
           plot = multiplot_te, width = 20, height = 15, dpi = 300)
    message("Saved combined TE expression counts plot with ", length(TE_expression_list), " panels")
}

if(length(TE_expression_dmr_plot_list) > 0){
    multiplot_dmr <- plot_grid(plotlist = TE_expression_dmr_plot_list, ncol = 3, nrow = 3, labels = "AUTO")
    ggsave(paste0(base_dir, "/TE_expression_in_DMR_barplot_faceted_multiplot.pdf"), 
           plot = multiplot_dmr, width = 20, height = 15)
    ggsave(paste0(base_dir, "/TE_expression_in_DMR_barplot_faceted_multiplot.png"), 
           plot = multiplot_dmr, width = 20, height = 15, dpi = 300)
    message("Saved combined TE-DMR overlap plot with ", length(TE_expression_dmr_plot_list), " panels")
}

if(length(TE_transcripts_dmr_list) > 0){
    TE_transcripts_dmr_df <- dplyr::bind_rows(TE_transcripts_dmr_list) %>%
        mutate(comparison = factor(comparison, levels = sample_names))
    multiplot_dmr <- ggplot(TE_transcripts_dmr_df, aes(x = reorder(TE_family_id, -count), y = count, fill = DMR_direction)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~comparison, scales = "free_x") +
        labs(title = "", x = "", 
             y = "Number of Repeat Transcripts Overlapping with DMRs") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          plot.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          strip.text = element_text(size = 20)) 

    ggsave(paste0(base_dir, "/TE_DMR_overlap_counts_barplot_faceted_multiplot.pdf"), 
           plot = multiplot_dmr, width = 20, height = 15)
    ggsave(paste0(base_dir, "/TE_DMR_overlap_counts_barplot_faceted_multiplot.png"), 
           plot = multiplot_dmr, width = 20, height = 15, dpi = 300)
    message("Saved combined TE-DMR overlap plot with ", nrow(TE_transcripts_dmr_df), " panels")
}

if(length(TE_dmr_overlap_list) > 0){
    TE_dmr_overlap_df <- dplyr::bind_rows(TE_dmr_overlap_list) %>%
        mutate(comparison = factor(comparison, levels = sample_names))
    multiplot_dmr <- ggplot(TE_dmr_overlap_df, aes(x = reorder(TE_family_id, -count), y = count, fill = DMR_direction)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~comparison, scales = "free_x") +
        labs(title = "", x = "", 
             y = "Number of DMRs Overlapping with Repeats") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          plot.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          strip.text = element_text(size = 20)) 

    ggsave(paste0(base_dir, "/DMR_TE_overlap_counts_barplot_faceted_multiplot.pdf"), 
           plot = multiplot_dmr, width = 20, height = 15)
    ggsave(paste0(base_dir, "/DMR_TE_overlap_counts_barplot_faceted_multiplot.png"), 
           plot = multiplot_dmr, width = 20, height = 15, dpi = 300)
    message("Saved combined TE-DMR overlap plot with ", nrow(TE_dmr_overlap_df), " panels")
}

if(length(TE_expression_in_dmr_plot_list) > 0){
    multiplot_dmr <- plot_grid(plotlist = TE_expression_in_dmr_plot_list, ncol = 3, nrow = 3, labels = "AUTO")
    ggsave(paste0(base_dir, "/regulated_TE_expression_in_DMR_barplot_faceted_multiplot.pdf"), 
           plot = multiplot_dmr, width = 20, height = 15)
    ggsave(paste0(base_dir, "/regulated_TE_expression_in_DMR_barplot_faceted_multiplot.png"), 
           plot = multiplot_dmr, width = 20, height = 15, dpi = 300)
    message("Saved combined regulated TE-DMR overlap plot with ", length(TE_expression_in_dmr_plot_list), " panels")
}

# combine dmr_annotation_count_list into a single data frame
dmr_annotation_count_df <- dplyr::bind_rows(dmr_annotation_count_list) %>%
    mutate(comparison = factor(comparison, levels = sample_names))

# plot dmr_annotation_count_df, facet by comparison
p_dmr_annotation <- ggplot(dmr_annotation_count_df, aes(x = annotation, y = count, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~comparison, scales = "free_x") +
    labs(title = "", x = "", y = "Number of DMRs", fill = "Direction") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          plot.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          strip.text = element_text(size = 20)) 

ggsave(paste0(base_dir, "/DMR_annotation_barplot_faceted_multiplot.pdf"), 
       plot = p_dmr_annotation, width = 20, height = 15)
ggsave(paste0(base_dir, "/DMR_annotation_barplot_faceted_multiplot.png"), 
       plot = p_dmr_annotation, width = 20, height = 15, dpi = 300)
message("Saved combined DMR annotation plot with ", length(dmr_annotation_count_list), " panels")
           
    