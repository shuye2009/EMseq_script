rm(list=ls())

library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)
library(circlize)

hd <- "C:/PROJECTS/Shane/Harding_250611/wo_chrY"
setwd(hd)

outdir <- file.path(hd, "downstream_analysis")
if(!dir.exists(outdir)) dir.create(outdir)

col_list <- list(Count=list(colid=2, ranges=c(0, 100, 1000), type="#DMRs"),
                 Enrichment=list(colid=5, ranges=c(-100, 0, 100), type="-LogPval"))

for(dmr_method in c("dmrseq", "DSS")){
    homer_folder <- "HOMER"
    subfolder <- "cutoff0.05"
    base_dir <- file.path(hd, dmr_method)
    if(dmr_method == "dmrseq") base_dir <- file.path(base_dir, subfolder)
    if(dmr_method == "DSS") subfolder <- ""

    sample_dirs <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)

    sample_names <- sample_dirs

    names(sample_dirs) <- sample_names

    for(aname in names(col_list)){
        enrichment_list <- list()
        for (comp in names(sample_dirs)){
            
            dir <- file.path(base_dir, sample_dirs[comp], homer_folder)
            for (direction in c("hyper", "hypo")){
            
                stat_file <- file.path(dir, direction, 
                                    paste0("DMRs_", direction, "_annot_stat.tab"))
                annot_stat <- read.delim2(stat_file, header = TRUE)
                enrichment <- annot_stat[15:46, col_list[[aname]][["colid"]], drop = FALSE]
                
                rownames(enrichment) <- annot_stat[15:46, 1]
                colnames(enrichment) <- paste(comp, direction, sep=":")
            
                enrichment_list[[comp]][[direction]] <- enrichment
            }
            
        }
        
        enrichment_df <- dplyr::bind_cols(enrichment_list)
        enrichment_df <- type.convert(as.matrix(enrichment_df), as.is=TRUE)
        if(aname == "Enrichment"){
            enrichment_df <- -1 * enrichment_df
        }
        
        print(head(enrichment_df))
        
        annot_df <- as.data.frame(t(bind_cols(strsplit(colnames(enrichment_df), split=":"))))
        colnames(annot_df) <- c("Comparison", "Direction")
        annot_df <- annot_df %>% 
            mutate(test = sapply(strsplit(Comparison, "_"), function(x) x[1])) %>%
            mutate(dose = sapply(strsplit(test, "Gy"), function(x) x[1])) %>%
            mutate(dose = sub("IR", "", dose)) %>%
            mutate(dose = paste0(dose, "Gy")) %>%
            mutate(time = sapply(strsplit(test, "Gy"), function(x) x[2])) 

        vannot <- HeatmapAnnotation(df = annot_df[, c("dose", "time")], which = "column")
        
        write.table(enrichment_df, 
            file.path(outdir, paste0("Homer_annotation_heatmap_", aname, "_", dmr_method, subfolder, ".tab")), 
            row.names = TRUE, 
            col.names = NA, 
            sep = "\t", 
            quote = FALSE)

        pdf(file.path(outdir, paste0("Homer_annotation_heatmap_", aname, "_", dmr_method, subfolder, ".pdf")), width = 12, height = 12)
        p <- ComplexHeatmap::Heatmap(enrichment_df, 
                                name = col_list[[aname]][["type"]],
                                top_annotation = vannot, 
                                show_column_names = TRUE,
                                cluster_rows = TRUE,
                                cluster_columns = FALSE,
                                column_names_rot = 90,
                                column_title = paste0(dmr_method, " DMRs ", aname),
                                column_names_gp = gpar(fontsize = 8),
                                col = colorRamp2(col_list[[aname]][["ranges"]], 
                                                colors=c("black", "grey", "orange")))
        print(p)
        dev.off()

    }
}