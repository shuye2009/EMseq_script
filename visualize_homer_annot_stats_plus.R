library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)
library(circlize)

hd <- "C:/PROJECTS/Shane/Harding_combined/h4h/cgmaptools"
setwd(hd)

datadir <- file.path(hd, "cgmaptools_dmr")
outdir <- file.path(hd, "results")
if(!dir.exists(outdir)) dir.create(outdir)

ERV_df <- read.delim(file.path(outdir, "ERV_enrichment_test.tab"))

col_list <- list(Ndmr=list(colid=2, ranges=c(0, 1000, 10000), type="DMRs", ERVcid=5),
                 Enrich=list(colid=5, ranges=c(-100, 0, 100), type="LogPval", ERVcid=8))

comparisons <- c("Irradiation", "Irradiation+M", "Mutation", "Mutation+I")
names(comparisons) <- c("WT-NIR_vs_WT-IR", "R172K-NIR_vs_R172K-IR", "WT-NIR_vs_R172K-NIR", "WT-IR_vs_R172K-IR")

for(aname in names(col_list)){
  enrichment_list <- list()
  for (comp in names(comparisons)){
    for(context in c("CG", "CHG", "CHH")){
      dir <- file.path(datadir, paste0(comp,"_dmr.", context))
      for (direction in c("hyper", "hypo")){
      
        stat_file <- file.path(dir, "HOMER", direction, 
                               paste0("DMRs_", direction, "_annot_stat.tab"))
        annot_stat <- read.delim2(stat_file, header = TRUE)
        enrichment <- annot_stat[15:46, col_list[[aname]][["colid"]], drop = FALSE]
    
        rownames(enrichment) <- annot_stat[15:46, 1]
        colnames(enrichment) <- paste(comp, context, direction, sep=":")
        
        
        ERV_enrichment <- ERV_df %>% 
          filter(Comparison == comp, Context == context, Direction == direction,
                 Annotation != "CpGi") 
        erv_enrichment <- ERV_enrichment[, col_list[[aname]][["ERVcid"]], drop=FALSE]
        rownames(erv_enrichment) <- ERV_enrichment[, "Annotation", drop=TRUE]
        colnames(erv_enrichment) <- paste(comp, context, direction, sep=":")
        
        enrichment <- rbind(enrichment, erv_enrichment)
        
        enrichment_list[[comp]][[context]][[direction]] <- enrichment
      }
    }
  }
  
  enrichment_df <- dplyr::bind_cols(enrichment_list)
  enrichment_df <- type.convert(as.matrix(enrichment_df), as.is=TRUE)
  
  print(head(enrichment_df))
  
  annot_df <- as.data.frame(t(bind_cols(strsplit(colnames(enrichment_df), split=":"))))
  colnames(annot_df) <- c("Comparison", "Context", "Direction")
  annot_df$Comparison <- comparisons[annot_df[,"Comparison", drop=TRUE]]
  vannot <- HeatmapAnnotation(df = annot_df, which = "column")
  
  pdf(file.path(outdir, paste0("Homer_annotation_heatmap", aname, ".pdf")), width = 8, height = 10)
  p <- ComplexHeatmap::Heatmap(enrichment_df, 
                          name = col_list[[aname]][["type"]],
                          top_annotation = vannot, 
                          show_column_names = TRUE,
                          cluster_rows = TRUE,
                          cluster_columns = FALSE,
                          column_names_rot = 90,
                          column_names_gp = gpar(fontsize = 8),
                          col = colorRamp2(col_list[[aname]][["ranges"]], 
                                           colors=c("orange", "white", "black")))
  print(p)
  dev.off()

}
