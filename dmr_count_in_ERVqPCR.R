library(GenomicPlot)
library(bit64)
library(ComplexHeatmap)

hd <- "C:/PROJECTS/Shane/Harding_combined/h4h/cgmaptools"
setwd(hd)

dataDir <- file.path(hd, "cgmaptools_dmr")
outDir <- file.path(hd, "results")
if(!dir.exists(outDir)) dir.create(outDir)

dmrDirs <- list.dirs(dataDir, recursive = FALSE)

backgroundBed <- c(ERVqPCR="C:/PROJECTS/Shane/resource/ERV_primers_sequence_matched_loci.bed")
ERVnames <- read.delim(backgroundBed, header = FALSE)[,4]
  

genome <- "hg38"
comparisons <- c("Irradiation", "Irradiation+M", "Mutation", "Mutation+I")
names(comparisons) <- c("WT-NIR_vs_WT-IR", "R172K-NIR_vs_R172K-IR", "WT-NIR_vs_R172K-NIR", "WT-IR_vs_R172K-IR")

targetBed <- c(target="C:/PROJECTS/Shane/Harding_240908/h4h/cgmaptools/cgmaptools_dmr/WT-NIR_vs_R172K-NIR_dmr.CHH/HOMER/DMRs_hyper.bed")

overlap_count <- function(targetBed, backgroudBed, ignore.strand = TRUE){
  
  bedimportParams <- setImportParams(
    offset = 0, fix_width = 0, fix_point = "center", norm = FALSE,
    useScore = FALSE, outRle = FALSE, useSizeFactor = FALSE, genome = genome
  )
  
  target_gr <- handle_bed(targetBed, importParams = bedimportParams)$query
  background_gr <- handle_bed(backgroundBed, importParams = bedimportParams)$query
  
  counts <- GenomicRanges::countOverlaps(background_gr, target_gr, ignore.strand = ignore.strand)
  
  return(counts)
}

sink(file.path(outDir, "DMR_in_ERVqPCR_count.tab"))
cat(paste("Comparison", "Context", "Direction", paste(ERVnames, collapse = "\t"), sep = "\t"))
cat("\n")

count_list <- list()
for (comp in names(comparisons)){
  for(context in c("CG", "CHG", "CHH")){
    dir <- file.path(dataDir, paste0(comp,"_dmr.", context))
    for (direction in c("hyper", "hypo")){
      #print(paste(comp, context, direction, sep = " : "))
      stat_file <- file.path(dir, "HOMER", paste0("DMRs_", direction, ".bed"))
      names(stat_file) <- paste(comp, context, direction, sep=":")
      count <- overlap_count(targetBed=stat_file, backgroudBed, ignore.strand = TRUE)
      
      
      count_list[[comp]][[context]][[direction]] <- count
      cat(paste0(paste(c(comp, context, direction, unlist(count)), 
                       collapse="\t"), "\n"))
    }
  }
}

sink()

count_df <- readr::read_tsv(file.path(outDir, "DMR_in_ERVqPCR_count.tab"), 
                            col_types = c(rep("c", 3), rep("n",14)))
print(head(count_df))

annot_df <- as.data.frame(count_df[, c("Comparison", "Context", "Direction")])
data_df <- t(count_df[, -c(1:3)])
colnames(data_df) <- apply(annot_df, 1, function(x){
  paste(unlist(x), collapse="_")
})
vannot <- HeatmapAnnotation(df = annot_df, which = "column")

pdf(file.path(outDir, paste0("DMR_count_in_ERVqPCR_heatmap.pdf")), width = 8, height = 10)
p <- ComplexHeatmap::Heatmap(data_df, 
                             name = "Count",
                             top_annotation = vannot, 
                             show_column_names = TRUE,
                             column_names_rot = 90,
                             column_names_gp = gpar(fontsize = 8),
                             cluster_rows = TRUE,
                             cluster_columns = FALSE,
                             row_names_gp = gpar(fontsize = 8),
                             heatmap_legend_param = list(direction = "horizontal"),
                             col = circlize::colorRamp2(c(0,1,5), 
                                              colors=c("orange", "white", "black"))
                             )
draw(p, heatmap_legend_side = "top", annotation_legend_side = "bottom")
dev.off()


