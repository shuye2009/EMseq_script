# import DMRs and TE annotations
library(openxlsx)
library(GenomicRanges)
library(plyranges)
library(dplyr)
library(tidyr)

# Function to perform TE enrichment analysis (works with family_id or gene_id)
perform_te_enrichment <- function(dmr_te_data, analysis_name, output_dir, grouping_var = "family_id") {
    cat(paste("\n=== TE", stringr::str_to_title(grouping_var), "Enrichment Analysis for", analysis_name, "===\n"))
    
    if(nrow(dmr_te_data) == 0) {
        cat("No DMR-TE overlaps found for", analysis_name, "\n")
        return(NULL)
    }
    
    # Determine column names based on grouping variable
    dmrTE_col <- paste0("TE_", grouping_var)
    genomeTE_col <- grouping_var
    
    # Check if required columns exist
    if(!dmrTE_col %in% colnames(dmr_te_data)) {
        cat("Error: Column", dmrTE_col, "not found in DMR data\n")
        cat("Available columns:", paste(colnames(dmr_te_data), collapse = ", "), "\n")
        return(NULL)
    }
    
    if(!genomeTE_col %in% colnames(mcols(TE_gr))) {
        cat("Error: Column", genomeTE_col, "not found in TE genome data\n")
        cat("Available columns:", paste(colnames(mcols(TE_gr)), collapse = ", "), "\n")
        return(NULL)
    }
    
    # Count TE groups in DMRs
    group_counts_dmrTE <- table(dmr_te_data[[dmrTE_col]])
    cat(paste("TE", grouping_var, "found in DMRs:\n"))
    print(group_counts_dmrTE)
    
    # Calculate width-based statistics for DMR overlaps
    # Sum the width of actual overlaps between DMRs and TEs by group
    group_width_dmrTE <- tapply(dmr_te_data$overlap_width, dmr_te_data[[dmrTE_col]], sum)
    cat(paste("\nTotal genomic width (bp) of TE", grouping_var, "overlapping DMRs:\n"))
    print(group_width_dmrTE)
    
    # Calculate width-based statistics for genome background
    # Sum the width of all TEs in genome by group
    group_width_genomeTE <- tapply(width(TE_gr), mcols(TE_gr)[[genomeTE_col]], sum)
    
    # Count all TE groups in the genome (background)
    group_counts_genomeTE <- table(mcols(TE_gr)[[genomeTE_col]])
    
    # Perform enrichment analysis
    group_enrichment <- data.frame(
        TE_id = character(),
        dmrTE_count = integer(),
        genomeTE_count = integer(),
        total_dmrTE_count = integer(),
        total_genomeTE_count = integer(),
        dmrTE_width_bp = numeric(),
        genomeTE_width_bp = numeric(),
        total_dmrTE_width_bp = numeric(),
        total_genomeTE_width_bp = numeric(),
        expected_width_bp = numeric(),
        fold_enrichment_width = numeric(),
        fold_enrichment_count = numeric(),
        p_value = numeric(),
        stringsAsFactors = FALSE
    )
    
    total_dmrTE_overlaps <- nrow(dmr_te_data)
    total_genomeTE <- length(TE_gr)
    # Sum width of distinct DMRs only (to avoid counting DMRs multiple times when they overlap multiple TEs)
    total_dmrTE_width <- dmr_te_data %>% 
        distinct(DMR_seqnames, DMR_start, DMR_end, .keep_all = TRUE) %>% 
        pull(DMR_width) %>% 
        sum()
    total_genomeTE_width <- sum(as.numeric(width(TE_gr)))
    
    # Test each group found in DMRs
    for(group in names(group_counts_dmrTE)){
        dmrTE_count <- as.numeric(group_counts_dmrTE[group])
        genomeTE_count <- as.numeric(group_counts_genomeTE[group])
        dmrTE_width <- as.numeric(group_width_dmrTE[group])
        genomeTE_width <- as.numeric(group_width_genomeTE[group])
        
        # Expected width under null hypothesis
        expected_width <- (genomeTE_width * total_dmrTE_width) / total_genomeTE_width
        
        # Fold enrichment based on width
        fold_enrichment_width <- dmrTE_width / expected_width
        
        # Fold enrichment based on count (for comparison)
        expected_count <- (genomeTE_count * total_dmrTE_overlaps) / total_genomeTE
        fold_enrichment_count <- dmrTE_count / expected_count
        
        # Fisher's exact test using width (binned to nearest 100bp for computational efficiency)
        # This tests if the observed width is significantly different from expected
        dmrTE_width_bin <- round(dmrTE_width / 100)
        genomeTE_width_bin <- round(genomeTE_width / 100)
        total_dmrTE_width_bin <- round(total_dmrTE_width / 100)
        total_genomeTE_width_bin <- round(total_genomeTE_width / 100)
        
        contingency_matrix <- matrix(c(
            dmrTE_width_bin,
            total_dmrTE_width_bin - dmrTE_width_bin,
            genomeTE_width_bin - dmrTE_width_bin,
            total_genomeTE_width_bin - genomeTE_width_bin - total_dmrTE_width_bin + dmrTE_width_bin
        ), nrow = 2)
        
        # Only perform test if we have sufficient values
        if(all(contingency_matrix >= 0) && sum(contingency_matrix) > 0){
            fisher_test <- fisher.test(contingency_matrix, alternative = "greater")
            p_value <- fisher_test$p.value
        } else {
            p_value <- 1
        }
        
        # Add to results
        group_enrichment <- rbind(group_enrichment, data.frame(
            TE_id = group,
            dmrTE_count = dmrTE_count,
            genomeTE_count = genomeTE_count,
            total_dmrTE_count = total_dmrTE_overlaps,
            total_genomeTE_count = total_genomeTE,
            dmrTE_width_bp = dmrTE_width,
            genomeTE_width_bp = genomeTE_width,
            total_dmrTE_width_bp = total_dmrTE_width,
            total_genomeTE_width_bp = total_genomeTE_width,
            expected_width_bp = expected_width,
            fold_enrichment_width = fold_enrichment_width,
            fold_enrichment_count = fold_enrichment_count,
            p_value = p_value,
            stringsAsFactors = FALSE
        ))
    }
    
    # Adjust p-values for multiple testing
    if(nrow(group_enrichment) > 0) {
        group_enrichment$padj <- p.adjust(group_enrichment$p_value, method = "BH")
        
        # Sort by significance
        group_enrichment <- group_enrichment[order(group_enrichment$padj, group_enrichment$p_value), ]

        # join with dge_TE
        if(grouping_var == "gene_id") {
            group_enrichment <- dplyr::left_join(group_enrichment, dge_TE, by = c("TE_id"="gene_id"))

            # select only significant enrichments
            group_enrichment_sig <- group_enrichment[group_enrichment$padj.x < 0.05, ]
            # remove rows with NA in family_id
            group_enrichment_sig <- group_enrichment_sig[!is.na(group_enrichment_sig$family_id), ]  
            
            # Count gene regulation status by family_id from group_enrichment_sig table
            if(nrow(group_enrichment_sig) > 0){
                cat("\n=== Gene Regulation Count Analysis ===\n")
                
                # Create regulation status based on padj.y
                group_enrichment_sig$regulation_status <- ifelse(is.na(group_enrichment_sig$padj.y) | group_enrichment_sig$padj.y >= 0.05, 
                                                               "unchanged",
                                                               ifelse(group_enrichment_sig$log2FoldChange > 0, 
                                                                      "up_regulated", 
                                                                      "down_regulated"))
                
                # Count genes by family_id and regulation status
                regulation_counts <- group_enrichment_sig %>%
                    group_by(family_id, regulation_status) %>%
                    summarise(gene_count = dplyr::n(), .groups = 'drop') %>%
                    pivot_wider(names_from = regulation_status, 
                               values_from = gene_count, 
                               values_fill = 0)
                
                # Ensure all columns exist
                if(!"up_regulated" %in% colnames(regulation_counts)) regulation_counts$up_regulated <- 0
                if(!"down_regulated" %in% colnames(regulation_counts)) regulation_counts$down_regulated <- 0
                if(!"unchanged" %in% colnames(regulation_counts)) regulation_counts$unchanged <- 0
                
                # Reorder columns
                regulation_counts <- regulation_counts[, c("family_id", "up_regulated", "down_regulated", "unchanged")]
                
                # Add total count
                regulation_counts$total_genes <- regulation_counts$up_regulated + regulation_counts$down_regulated + regulation_counts$unchanged
                
                # Display the table
                cat("\nGene regulation counts by family_id:\n")
                print(regulation_counts)
                
                # Save the results
                regulation_output_file <- file.path(output_dir, paste0("gene_regulation_counts_by_family_", analysis_name, ".xlsx"))
                write.xlsx(regulation_counts, regulation_output_file)
                cat(paste("\nSaved gene regulation counts to:", basename(regulation_output_file), "\n"))
                
                # Print summary statistics
                cat("\nSummary statistics:\n")
                cat(paste("Total families analyzed:", nrow(regulation_counts), "\n"))
                cat(paste("Total up-regulated genes:", sum(regulation_counts$up_regulated), "\n"))
                cat(paste("Total down-regulated genes:", sum(regulation_counts$down_regulated), "\n"))
                cat(paste("Total unchanged genes:", sum(regulation_counts$unchanged), "\n"))
                cat(paste("Total genes:", sum(regulation_counts$total_genes), "\n"))
                
            } else {
                cat("\nNo significant enrichments found for gene regulation analysis.\n")
            }
        } 
        
        # Display results
        print(head(group_enrichment))
        
        # Save enrichment results
        output_file <- file.path(output_dir, paste0("TE_", grouping_var, "_enrichment_", analysis_name, ".xlsx"))
        write.xlsx(group_enrichment, output_file)
        cat(paste("Saved enrichment results to:", basename(output_file), "\n"))
        
        # Print summary of significant enrichments (based on width)
        sig_enrichments <- group_enrichment[group_enrichment$padj < 0.05 & group_enrichment$fold_enrichment_width > 1, ]
        if(nrow(sig_enrichments) > 0){
            cat(paste("Found", nrow(sig_enrichments), "significantly enriched TE", grouping_var, "(padj < 0.05, based on width):\n"))
            for(i in seq_len(nrow(sig_enrichments))){
                cat(paste(" -", sig_enrichments$group_id[i], 
                        ": Fold enrichment (width) =", round(sig_enrichments$fold_enrichment_width[i], 2),
                        ", Fold enrichment (count) =", round(sig_enrichments$fold_enrichment_count[i], 2),
                        ", padj =", format(sig_enrichments$padj[i], scientific = TRUE, digits = 3), "\n"))
            }
        } else {
            cat(paste("No significantly enriched TE", grouping_var, "found (padj < 0.05)\n"))
        }
        
        return(group_enrichment)
    } else {
        cat("No enrichment results generated\n")
        return(NULL)
    }
}

# Apply enrichment analysis to different subsets

# Function to run both family_id and gene_id enrichment analyses
run_enrichment_analyses <- function(data, analysis_name, output_dir) {
    # Family-level enrichment
    cat("\n" , rep("=", 80), "\n")
    cat(paste("ANALYZING", toupper(analysis_name), "- FAMILY LEVEL\n"))
    cat(rep("=", 80), "\n")
    family_enrichment <- perform_te_enrichment(data, analysis_name, output_dir, "family_id")
    
    # Gene-level enrichment
    cat("\n" , rep("=", 80), "\n")
    cat(paste("ANALYZING", toupper(analysis_name), "- GENE LEVEL\n"))
    cat(rep("=", 80), "\n")
    gene_enrichment <- perform_te_enrichment(data, analysis_name, output_dir, "gene_id")
    
    return(list(family = family_enrichment, gene = gene_enrichment))
}

#### MAIN analysis 

# import TE annotations
TE_gtf <- "C:/PROJECTS/Shane/RNAseq/GRCh38_GENCODE_rmsk_TE.gtf/GRCh38_GENCODE_rmsk_TE.gtf"
TE_gr <- plyranges::read_gff(TE_gtf)
 
dmr_method <- "TE_targeted" #"dmrseq", "TE_targeted" "DSS"
subfolder <- "cutoff0.05" #"cutoff0.1", "cutoff0.05" only for dmrseq

base_dir <- file.path("C:/PROJECTS/Shane/Harding_250611/wo_chrY/", dmr_method)
if(dmr_method == "dmrseq") base_dir <- file.path(base_dir, subfolder)

dmr_folder <- "DMRs"
if(dmr_method == "TE_targeted") dmr_folder <- "Targeted"


sample_dirs <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)

if(dmr_method == "dmrseq"){
    sample_names <- gsub("dmrseq_", "", sample_dirs)
    sample_names <- gsub("_", "_vs_", sample_names)
    
}else{
    sample_names <- sample_dirs
}

names(sample_dirs) <- sample_names
#sample_names <- c("IR2Gy6d_vs_NIR", "IR2Gy24h_vs_NIR", "IR10Gy6d_vs_NIR", "IR10Gy24h_vs_NIR", "IR2Gy6d_vs_IR2Gy24h", "IR10Gy6d_vs_IR10Gy24h", "IR10Gy6d_vs_IR2Gy6d", "IR10Gy24h_vs_IR2Gy24h")

for(sample_name in sample_names){
    message("Processing ", sample_name, " ##############################################################################################")
    # remove leading "IR" from sample name
    SA <- gsub("^IR", "", sample_name)
    SA <- gsub("_IR", "_", SA)
    
    # import differential transcript expression results
    dte_TE <- read.csv(paste0("C:/PROJECTS/Shane/RNAseq/20250604_LH00244_0317_A22WFMGLT3_Harding_Shreya_RNAseq/results_TE/deseq2/single_factor_analysis_TE_transcript/", SA, "/", SA, "_DESeq2_results_TE.csv"))
    head(dte_TE)

    # import differential gene expression results
    dge_TE <- read.csv(paste0("C:/PROJECTS/Shane/RNAseq/20250604_LH00244_0317_A22WFMGLT3_Harding_Shreya_RNAseq/results_TE/deseq2/single_factor_analysis_TE_gene_level/", SA, "/", SA, "_DESeq2_results_TE.csv"))
    head(dge_TE)

    # import DMRs and TE annotations
    working_dir <- file.path(base_dir, sample_dirs[sample_name], dmr_folder)

    if(dmr_method == "TE_targeted"){
        dmr <- read.table(file.path(working_dir, "/significant_diff_GRCh38_GENCODE_rmsk_TE.tab"), header = TRUE, sep = "\t")
    } else {
        dmr <- read.xlsx(file.path(working_dir, "/DMRs_annotated.xlsx"))
    }
    # change column names "difference" to "diff.Methy" to keep consistent with DSS
    if(dmr_method == "dmrseq") colnames(dmr) <- gsub("difference", "diff.Methy", colnames(dmr)) 
    if(dmr_method == "TE_targeted") colnames(dmr) <- gsub("meth.diff", "diff.Methy", colnames(dmr)) 
    
    dmr_gr <- makeGRangesFromDataFrame(dmr, keep.extra.columns = TRUE)
    
    # find overlap between DMRs and TE
    overlap <- findOverlaps(dmr_gr, TE_gr)

    # extract each DMR that overlaps with a TE
    overlap_dmr <- dmr_gr[queryHits(overlap)]

    # extract each TE that overlaps with a DMR
    overlap_TE <- TE_gr[subjectHits(overlap)]

    head(overlap_dmr)
    head(overlap_TE)

    # convert to data frame
    overlap_dmr_df <- data.frame(overlap_dmr)
    overlap_TE_df <- data.frame(overlap_TE)

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


    # merge DMRs and TE annotations with differential transcript expression results
    overlap_dmr_TE <- dplyr::left_join(overlap_dmr_TE, dte_TE, by = c("TE_transcript_id"="transcript_id", "TE_gene_id"="gene_id", "TE_family_id"="family_id", "TE_class_id"="class_id"))

    # save to excel
    write.xlsx(overlap_dmr_TE, paste0(working_dir, "/DMRs_annotated_with_TE.xlsx"))


    # 1. All DMRs
    enrichment_all <- run_enrichment_analyses(overlap_dmr_TE, "all_DMRs", working_dir)

    # 2. Hypermethylated DMRs only
    if("DMR_diff.Methy" %in% colnames(overlap_dmr_TE)) {
        hyper_dmrs <- overlap_dmr_TE[overlap_dmr_TE$DMR_diff.Methy > 0, ]
        enrichment_hyper <- run_enrichment_analyses(hyper_dmrs, "hypermethylated_DMRs", working_dir)
        
        # 3. Hypomethylated DMRs only
        hypo_dmrs <- overlap_dmr_TE[overlap_dmr_TE$DMR_diff.Methy < 0, ]
        enrichment_hypo <- run_enrichment_analyses(hypo_dmrs, "hypomethylated_DMRs", working_dir)
    } else {
        cat("\nNote: diff.Methy column not found. Skipping hyper/hypo analysis.\n")
        cat("Available columns:", paste(colnames(overlap_dmr_TE), collapse = ", "), "\n")
    }

    # Select TE_transcripts that are not NA in baseMean
    dmr_te_dge <- overlap_dmr_TE[!is.na(overlap_dmr_TE$baseMean), ] %>%
        dplyr::arrange(padj)
    dim(dmr_te_dge)
    head(dmr_te_dge)

    # Output to excel
    write.xlsx(dmr_te_dge, paste0(working_dir, "/DMRs_annotated_with_TE_dge.xlsx"))

}
 
# create DMR and DML circos plot 
source("create_DMR_circos_plots.R")
source("create_DML_circos_plots.R")

generate_dmr_circos_plots(base_dir, sample_dirs, dmr_folder)
if(dmr_method == "DSS") generate_dml_circos_plots(base_dir, sample_dirs)
