library(GenomicPlot)

resource_dir <- "C:/PROJECTS/Shane/resource"
gtffile <- file.path(resource_dir, "gencode.v45.primary_assembly.annotation.gtf")
                  
genome <- "GRCh38"
overwrite <- FALSE

if (file.exists(file.path(resource_dir, paste0(genome, "_txdb.sql"))) && !overwrite) {
  txdb <- AnnotationDbi::loadDb(file.path(resource_dir, paste0(genome, "_txdb.sql")))
  print(class(txdb))
} else {
  txdb <- custom_TxDb_from_GTF(gtffile, genome)
  AnnotationDbi::saveDb(txdb, file = file.path(resource_dir, paste0(genome, "_txdb.sql")))
}

if (file.exists(file.path(resource_dir, paste0(genome, "_metaFeatures.rds"))) && !overwrite) {
  metaFeatures <- readRDS(file.path(resource_dir, paste0(genome, "_metaFeatures.rds")))
  print(class(metaFeatures))
} else {
  metaFeatures <- prepare_5parts_genomic_features(txdb = txdb, longest = TRUE, 
                                                  meta = TRUE, nbins = 100, 
                                                  fiveP = -2000, threeP = 2000)
  saveRDS(metaFeatures, file = file.path(resource_dir, paste0(genome, "_metaFeatures.rds")))
}

if (file.exists(file.path(resource_dir, paste0(genome, "_Features.rds"))) && !overwrite) {
  Features <- readRDS(file.path(resource_dir, paste0(genome, "_Features.rds")))
  print(class(Features))
} else {
  Features <- prepare_5parts_genomic_features(txdb = txdb, longest = TRUE, 
                                              meta = FALSE, nbins = 100, 
                                              fiveP = -2000, threeP = 2000)
  saveRDS(Features, file = file.path(resource_dir, paste0(genome, "_Features.rds")))
}

setwd("C:/PROJECTS/Shane/Harding_combined/local")

exp_list <- c(treat="WT-NIR_vs_WT-IR", mutant="WT-NIR_vs_R172K-NIR")
for(exp in names(exp_list)){
  comparison <- exp_list[exp]
  dmr_dirs <- c(Jan = paste0("C:/PROJECTS/Shane/Harding_240124/h4h/cgmaptools/cgmaptools_dmr/", comparison, "_dmr.CG/HOMER"),
                Aug =  paste0("C:/PROJECTS/Shane/Harding_240908/h4h/cgmaptools/cgmaptools_dmr/", comparison, "_dmr.CG/HOMER"),
                Combined =  paste0("C:/PROJECTS/Shane/Harding_combined/h4h/cgmaptools/cgmaptools_dmr/", comparison, "_dmr.CG/HOMER"),
                CombinedDM = paste0("C:/PROJECTS/Shane/Harding_combined/h4h/DMRichR_", exp, "/HOMER"),
                JanDM = paste0("C:/PROJECTS/Shane/Harding_240124/h4h/DMRichR_", exp, "/HOMER"))
  
  
  for (direction in c("hyper", "hypo")){
    queryfiles <- file.path(dmr_dirs, paste0("DMRs_", direction, ".bed"))
    names(queryfiles) <- names(dmr_dirs)
    
    importParams <- setImportParams(
      offset = 0, fix_width = 0, fix_point = "start", norm = FALSE, skip = 0, val = 4,
      useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg38"
    )
    
    plot_locus(
      queryFiles = queryfiles,
      centerFiles = "gene",
      txdb = txdb,
      ext = c(-2000, 2000),
      hl = c(-1000, 300),
      shade = TRUE,
      smooth = TRUE,
      importParams = importParams,
      verbose = FALSE,
      binSize = 20,
      refPoint = "start",
      Xlab = "TSS",
      Ylab = "DMRs/base/gene",
      inputFiles = NULL,
      stranded = TRUE,
      heatmap = TRUE,
      scale = FALSE,
      outPrefix = paste("cgmaptools_dmrCG_at_TSS", direction, exp, sep="_"),
      rmOutlier = 0,
      transform = NA,
      statsMethod = "wilcox.test",
      heatRange = c(0, 0.1, 0.2),
      hw = c(8, 8),
      nc = 5
    )
    
    plot_locus(
      queryFiles = queryfiles,
      centerFiles = queryfiles["Combined"],
      txdb = txdb,
      ext = c(-2000, 2000),
      hl = c(-1000, 1000),
      shade = TRUE,
      smooth = TRUE,
      importParams = importParams,
      verbose = FALSE,
      binSize = 20,
      refPoint = "center",
      Xlab = "DMR",
      Ylab = "DMRs/bin",
      inputFiles = NULL,
      stranded = TRUE,
      heatmap = TRUE,
      scale = FALSE,
      outPrefix = paste("cgmaptools_dmrCG_comparison_", direction, exp, sep="_"),
      rmOutlier = 0,
      transform = NA,
      statsMethod = "wilcox.test",
      heatRange = c(0,1,2),
      hw = c(8, 8),
      nc = 5
    )
  }
}
