args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 3) {
  stop("Usage: Rscript centromere_region_plot.R FALSE GENE CpG")
}
LOCAL = args[1] # TRUE or FALSE
TARGET = args[2] # c("LOCUS", "GENE", "METAGENE", "Centromere", "CpGIsland", "ERV", "ERVqPCR")
CONTEXT = args[3] # one of "CpG", "CHG", "CHH"
COMPARISON = args[4] # "WT-NIR_vs_WT-IR"

local_home_dir <- "C:/PROJECTS/Shane/Harding_240124"
resource_dir <- ifelse(LOCAL, file.path(local_home_dir, "resource"), "~/resource")
data_dir <- ifelse(LOCAL, file.path(local_home_dir, paste0("h4h/cgmaptools/cgmaptools_dmr/", COMPARISON, "_dmr.", CONTEXT, "/HOMER")), 
                   "/cluster/projects/hardinggroup/Shuye/EMseq/Harding_240908/alignment/MethylDackel/report/results/dmr")

wd <- ifelse(LOCAL, file.path(local_home_dir, "local"),
             "/cluster/projects/hardinggroup/Shuye/EMseq/Harding_240908/analysis/GenomicPlot")
setwd(wd)

library(GenomicPlot)

gtffile <- ifelse(LOCAL, file.path(resource_dir, "gencode.v45.primary_assembly.annotation.gtf"),
                  "/cluster/projects/hardinggroup/Shuye/reference_genome/human/GRCh38/gencode.v45.primary_assembly.annotation.gtf")
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

centerfiles <- c(Centromere = file.path(resource_dir, "GRCh38_centromere_UCSC_simplified.bed"),
                 CpGIsland = file.path(resource_dir, "cpgIsland_hg38_ucsc.bed"),
                 ERV = file.path(resource_dir, "ERV_hg38_ucsc.bed"),
                 ERVqPCR = file.path(resource_dir, "ERV_primers_sequence_matched_loci.bed")
                 )
nbins <- c(1000, 100, 100, 100)
fiveP <- c(-50000, -500, -300, -500)
threeP <- c(50000, 500, 300, 500)
names(nbins) <- names(fiveP) <- names(threeP) <- names(centerfiles)

importParams <- setImportParams(
  offset = 0, fix_width = 0, fix_point = "start", norm = FALSE, skip = 1, val = 4,
  useScore = TRUE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg38"
)

queryfiles <- c("DMRs_hyper.bed", "DMRs_hypo.bed")
shortName <- paste(COMPARISON, CONTEXT, c("hyper", "hypo"), sep = "_")
queryfiles <- file.path(data_dir, queryfiles)
names(queryfiles) <- shortName

print(queryfiles)

if(TARGET %in% c("Centromere", "CpGIsland", "ERV", "ERVqPCR")){
  plot_region(
    queryFiles = queryfiles,
    centerFiles = centerfiles[TARGET],
    inputFiles = NULL,
    nbins = nbins[TARGET],
    heatmap = TRUE,
    scale = FALSE,
    regionName = TARGET,
    importParams = importParams,
    verbose = FALSE,
    fiveP = fiveP[TARGET],
    threeP = threeP[TARGET],
    smooth = TRUE,
    transform = NA,
    stranded = FALSE,
    heatRange = c(0,1,2),
    outPrefix = paste0("DMRs_in_", TARGET, "_", CONTEXT, "_", COMPARISON),
    Ylab = paste(CONTEXT, "DMR density"),
    rmOutlier = 0,
    nc = 5
  )
}

if(TARGET == "METAGENE"){
  plot_5parts_metagene(
    queryFiles = queryfiles,
    gFeatures_list = list(metagene = metaFeatures),
    inputFiles = NULL,
    importParams = importParams,
    verbose = FALSE,
    transform = NA,
    smooth = TRUE,
    scale = FALSE,
    stranded = TRUE,
    outPrefix = paste0("Methylation_level_metagene_", CONTEXT),
    heatmap = TRUE,
    heatRange = NULL,
    rmOutlier = 0,
    Ylab = "Methylation/base/gene",
    hw = c(10, 10),
    nc = 5
  )
}

if(TARGET == "GENE"){
  plot_5parts_metagene(
    queryFiles = queryfiles,
    gFeatures_list = list(gene = Features),
    inputFiles = NULL,
    importParams = importParams,
    verbose = FALSE,
    transform = NA,
    smooth = TRUE,
    scale = FALSE,
    stranded = TRUE,
    outPrefix = paste0("Methylation_level_gene_", CONTEXT),
    heatmap = TRUE,
    heatRange = NULL,
    rmOutlier = 0,
    Ylab = "Methylation/base/gene",
    hw = c(10, 10),
    nc = 5
  )
}

if(TARGET == "LOCUS"){
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
    Ylab = "Methylation/base/gene",
    inputFiles = NULL,
    stranded = TRUE,
    heatmap = TRUE,
    scale = FALSE,
    outPrefix = paste0("Methylation_level_TSS_", CONTEXT),
    rmOutlier = 0,
    transform = NA,
    statsMethod = "wilcox.test",
    heatRange = NULL,
    hw = c(8, 8),
    nc = 5
  )
}

if(TARGET == "CHROM"){
  chrs <- circlize::read.chromInfo(species = "hg38")$chromosome
  lapply(chrs, function(chr){
    
    importParams <- setImportParams(
      offset = 0, fix_width = 0, fix_point = "start", norm = FALSE, skip = 1, val = 4,
      useScore = TRUE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg38", chr = chr
    )
    
    plot_region(
      queryFiles = queryfiles,
      centerFiles = centerfiles["Centromere"],
      inputFiles = NULL,
      nbins = 1000,
      heatmap = TRUE,
      scale = FALSE,
      regionName = "Centromere",
      importParams = importParams,
      verbose = TRUE,
      fiveP = -500000,
      threeP = 500000,
      smooth = TRUE,
      transform = NA,
      stranded = FALSE,
      outPrefix = paste0("Methylation_level_in_centromere_", CONTEXT, "_", chr),
      Ylab = paste(CONTEXT, "methylation/chromosome"),
      rmOutlier = 0,
      nc = 5
    )
    
  })
}


