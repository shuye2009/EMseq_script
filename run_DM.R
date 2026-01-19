#
options(error=recover)
library(DMRichR)
setwd("C:/PROJECTS/Shane/Harding_240124/h4h/data")

#Global variables, same as the arguments to DM.R
genome = "hg38"
coverage = 1
perGroup =  0.75
minCpGs =  5
maxPerms =  4
maxBlockPerms = 4
cutoff = 0.05
testCovariate = "Phenotype"
adjustCovariate = "Treatment"
matchCovariate = NULL
cores = 1
GOfuncR = TRUE
sexCheck = FALSE
EnsDb = FALSE
postfix = "*.CG_report.txt.gz"
internet = TRUE

if (Biobase::package.version("dmrseq") %>% stringr::str_remove("1.") %>% 
    as.numeric() < 7.3) {
  warning(paste("Your version of dmrseq is out of date and contains a bug.", 
                "This bug won't affect the DMRichR run but could affect your custom follow up analyses.", 
                "See the install section of the DMRichR README for the code to manually update.", 
                "Read more about the issue: https://github.com/kdkorthauer/dmrseq/issues/37"))
}
options(scipen = 999)
options(readr.num_columns = 0)
stopifnot(genome %in% c("hg38", "hg19", "mm10", "mm9", "rheMac10", 
                        "rheMac8", "rn6", "danRer11", "galGal6", "bosTau9", "panTro6", 
                        "dm6", "susScr11", "canFam3", "TAIR10", "TAIR9"))
stopifnot(!is.null(testCovariate))
stopifnot(coverage >= 1)
nSamples <- openxlsx::read.xlsx("sample_info.xlsx", colNames = TRUE) %>% 
  nrow()
if (nSamples < maxPerms) {
  print(glue::glue("Warning: You have requested {maxPerms} permutations for the DMR analysis, \\\n                   which is more than the {nSamples} samples you have. \\\n                   maxPerms will now be changed to {nSamples}."))
  maxPerms <- nSamples
}
if (nSamples < maxBlockPerms) {
  print(glue::glue("Warning: You have requested {maxBlockPerms} permutations for the block analysis, \\\n                   which is more than the {nSamples} samples you have. \\\n                   maxBlockPerms will now be changed to {nSamples}."))
  maxBlockPerms <- nSamples
}
rm(nSamples)
print(glue::glue("genome = {genome}"))
print(glue::glue("coverage = {coverage}"))
print(glue::glue("perGroup = {perGroup}"))
print(glue::glue("minCpGs = {minCpGs}"))
print(glue::glue("maxPerms = {maxPerms}"))
print(glue::glue("maxBlockPerms = {maxBlockPerms}"))
print(glue::glue("cutoff = {cutoff}"))
print(glue::glue("testCovariate = {testCovariate}"))
print(glue::glue("adjustCovariate = {adjustCovariate}"))
print(glue::glue("matchCovariate = {matchCovariate}"))
print(glue::glue("cores = {cores}"))
print(glue::glue("sexCheck = {sexCheck}"))
print(glue::glue("EnsDb = {EnsDb}"))
print(glue::glue("GOfuncR = {GOfuncR}"))
cat("\n[DMRichR] Selecting annotation databases \t\t", format(Sys.time(), 
                                                              "%d-%m-%Y %X"), "\n")
DMRichR::annotationDatabases(genome = genome, EnsDb = EnsDb)
print(glue::glue("Saving Rdata..."))
dir.create("RData")
settings_env <- ls(all = TRUE)
save(list = settings_env, file = "RData/settings.RData")
bs.filtered <- DMRichR::processBismark(files = list.files(path = getwd(), 
                                       pattern = postfix), 
                                       meta = openxlsx::read.xlsx("sample_info.xlsx", colNames = TRUE) %>% 
                                         dplyr::mutate_if(is.character, as.factor), 
                                       testCovariate = testCovariate, adjustCovariate = adjustCovariate, 
                                       matchCovariate = matchCovariate, coverage = coverage, 
                                       cores = cores, perGroup = perGroup, sexCheck = sexCheck)
print(glue::glue("Saving Rdata..."))
save(bs.filtered, file = "RData/bismark.RData")
print(glue::glue("Building annotations for plotting..."))
if (is(TxDb, "TxDb")) {
  if(internet){
    annoTrack <- dmrseq::getAnnot(genome)
    saveRDS(annoTrack, file.path(getwd(), "RData/annoTrack.rds"))
  }else{
    annoTrack <- readRDS(file.path(getwd(), "RData/annoTrack.rds"))
  }
}else if (is(TxDb, "EnsDb")) {
  annoTrack <- GenomicRanges::GRangesList(CpGs = DMRichR::getCpGs(genome), 
                                          Exons = DMRichR::getExons(TxDb), compress = FALSE)
}
cat("\n[DMRichR] Getting bsseq background regions \t\t", 
    format(Sys.time(), "%d-%m-%Y %X"), "\n")
dir.create("Extra")
DMRichR::getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 1000) %>% 
                      write.table(file = "Extra/bsseq_background.csv", 
                                  sep = ",", quote = FALSE, row.names = FALSE)

cat("\n[DMRichR] Testing for blocks with dmrseq \t\t", format(Sys.time(), 
                                                              "%d-%m-%Y %X"), "\n")
start_time <- Sys.time()
tryCatch({
  blocks <- dmrseq::dmrseq(bs = bs.filtered, cutoff = cutoff, 
                           maxPerms = maxBlockPerms, testCovariate = testCovariate, 
                           adjustCovariate = adjustCovariate, matchCovariate = matchCovariate, 
                           block = TRUE, minInSpan = 500, bpSpan = 50000, maxGapSmooth = 1e+06, 
                           maxGap = 5000, minNumRegion = (minCpGs * 2), BPPARAM = BiocParallel::MulticoreParam(workers = cores))
  print(glue::glue("Selecting significant blocks..."))
  if (length(blocks) != 0) {
    blocks <- blocks %>% plyranges::mutate(direction = dplyr::case_when(stat > 
                                                                          0 ~ "Hypermethylated", stat < 0 ~ "Hypomethylated"), 
                                           difference = round(beta/pi * 100))
  }
  if (sum(blocks$qval < 0.05) == 0 & sum(blocks$pval < 
                                         0.05) != 0) {
    sigBlocks <- blocks %>% plyranges::filter(pval < 
                                                0.05)
  }
  else if (sum(blocks$qval < 0.05) >= 1) {
    sigBlocks <- blocks %>% plyranges::filter(qval < 
                                                0.05)
  }
  else if (sum(blocks$pval < 0.05) == 0 & length(blocks) != 
           0) {
    glue::glue("No significant blocks detected in {length(blocks)} background blocks")
  }
  else if (length(blocks) == 0) {
    glue::glue("No background blocks detected")
  }
  if (length(blocks) != 0) {
    print(glue::glue("Exporting block and background information..."))
    dir.create("Blocks")
    gr2bed(blocks, "Blocks/backgroundBlocks.bed")
    if (sum(blocks$pval < 0.05) > 0) {
      print(glue::glue("{length(sigBlocks)} significant blocks of differential methylation in {length(blocks)} background blocks"))
      gr2bed(sigBlocks, "Blocks/blocks.bed")
      print(glue::glue("Annotating and plotting blocks..."))
      pdf("Blocks/Blocks.pdf", height = 7.5, width = 11.5)
      dmrseq::plotDMRs(bs.filtered, regions = sigBlocks, 
                       testCovariate = testCovariate, annoTrack = annoTrack, 
                       regionCol = "#FF00001A", qval = FALSE, stat = FALSE)
      dev.off()
    }
  }
  print(glue::glue("Blocks timing..."))
  end_time <- Sys.time()
  end_time - start_time
  if (length(blocks) != 0) {
    if (sum(blocks$pval < 0.05) > 0) {
      print(glue::glue("Annotating blocks with gene symbols..."))
      sigBlocks %>% DMRichR::annotateRegions(TxDb = TxDb, 
                                             annoDb = annoDb) %T>% DMRichR::DMReport(regions = blocks, 
                                                                                     bs.filtered = bs.filtered, coverage = coverage, 
                                                                                     name = "blockReport") %>% openxlsx::write.xlsx(file = "Blocks/Blocks_annotated.xlsx")
    }
    print(glue::glue("Annotating background blocks with gene symbols..."))
    blocks %>% DMRichR::annotateRegions(TxDb = TxDb, 
                                        annoDb = annoDb) %>% openxlsx::write.xlsx(file = "Blocks/background_blocks_annotated.xlsx")
  }
  print(glue::glue("Saving RData..."))
  save(blocks, file = "RData/Blocks.RData")
}, error = function(error_condition) {
  print(glue::glue("Warning: Block analysis has produced an error"))
})
cat("\n[DMRichR] Testing for DMRs with dmrseq \t\t\t", format(Sys.time(), 
                                                              "%d-%m-%Y %X"), "\n")
start_time <- Sys.time()
regions <- dmrseq::dmrseq(bs = bs.filtered, cutoff = cutoff, 
                          minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = testCovariate, 
                          adjustCovariate = adjustCovariate, matchCovariate = matchCovariate, 
                          BPPARAM = BiocParallel::MulticoreParam(workers = cores))
print(glue::glue("Selecting significant DMRs..."))
regions <- regions %>% plyranges::mutate(direction = dplyr::case_when(stat > 
                                                                        0 ~ "Hypermethylated", stat < 0 ~ "Hypomethylated"), 
                                         difference = round(beta/pi * 100))
if (sum(regions$qval < 0.05) < 100 & sum(regions$pval < 0.05) != 0) {
  sigRegions <- regions %>% plyranges::filter(pval < 0.05)
}else if (sum(regions$qval < 0.05) >= 100) {
  sigRegions <- regions %>% plyranges::filter(qval < 0.05)
}else if (sum(regions$pval < 0.05) == 0) {
  stop(glue::glue("No significant DMRs detected in {length(regions)} background regions"))
}
print(glue::glue("Exporting DMR and background region information..."))
dir.create("DMRs")
gr2bed(sigRegions, "DMRs/DMRs.bed")
gr2bed(regions, "DMRs/backgroundRegions.bed")
if (sum(sigRegions$stat > 0) > 0 & sum(sigRegions$stat < 
                                       0) > 0) {
  print(glue::glue("Summary: There are {tidySigRegions} DMRs \\\n             ({tidyHyper}% hypermethylated, {tidyHypo}% hypomethylated) \\\n             from {tidyRegions} background regions consisting of {tidyCpGs} CpGs \\\n             assayed at {coverage}x coverage", 
                   tidySigRegions = length(sigRegions), tidyHyper = round(sum(sigRegions$stat > 
                                                                                0)/length(sigRegions), digits = 2) * 100, tidyHypo = round(sum(sigRegions$stat < 
                                                                                                                                                 0)/length(sigRegions), digits = 2) * 100, tidyRegions = length(regions), 
                   tidyCpGs = nrow(bs.filtered)))
}
print(glue::glue("DMR timing..."))
end_time <- Sys.time()
end_time - start_time
print(glue::glue("Saving Rdata..."))
save(regions, sigRegions, file = "RData/DMRs.RData")
print(glue::glue("Annotating DMRs and plotting..."))
pdf("DMRs/DMRs.pdf", height = 4, width = 8)
tryCatch({
  DMRichR::plotDMRs2(bs.filtered, regions = sigRegions, 
                     testCovariate = testCovariate, extend = (end(sigRegions) - 
                                                                start(sigRegions) + 1) * 2, addRegions = sigRegions, 
                     annoTrack = annoTrack, regionCol = "#FF00001A", lwd = 2, 
                     qval = FALSE, stat = FALSE, horizLegend = FALSE)
}, error = function(error_condition) {
  print(glue::glue("Warning: One (or more) of your DMRs can't be plotted, \\\n                      try again later by manually loading R Data and subsetting sigRegions"))
})
dev.off()
cat("\n[DMRichR] Annotating DMRs with gene symbols \t\t", 
    format(Sys.time(), "%d-%m-%Y %X"), "\n")
sigRegions %>% DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %T>% 
  DMRichR::DMReport(regions = regions, bs.filtered = bs.filtered, 
                    coverage = coverage, name = "DMReport") %>% openxlsx::write.xlsx(file = "DMRs/DMRs_annotated.xlsx")
print(glue::glue("Annotating background regions with gene symbols..."))
regions %>% DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>% 
  openxlsx::write.xlsx(file = "DMRs/background_annotated.xlsx")
cat("\n[DMRichR] Smoothing individual methylation values \t\t", 
    format(Sys.time(), "%d-%m-%Y %X"), "\n")
start_time <- Sys.time()
bs.filtered.bsseq <- bsseq::BSmooth(bs.filtered, BPPARAM = BiocParallel::MulticoreParam(workers = cores, 
                                                                                        progressbar = TRUE))
if (genome == "rn6") {
  bs.filtered.bsseq <- GenomeInfoDb::dropSeqlevels(bs.filtered.bsseq, 
                                                   "chrY", pruning.mode = "coarse")
  GenomeInfoDb::seqlevels(bs.filtered.bsseq)
}
bs.filtered.bsseq
print(glue::glue("Extracting individual smoothed methylation values of DMRs..."))
bs.filtered.bsseq %>% DMRichR::smooth2txt(regions = sigRegions, 
                                          txt = "DMRs/DMR_individual_smoothed_methylation.txt")
print(glue::glue("Extracting individual smoothed methylation values of background regions..."))
bs.filtered.bsseq %>% DMRichR::smooth2txt(regions = regions, 
                                          txt = "DMRs/background_region_individual_smoothed_methylation.txt")
print(glue::glue("Individual smoothing timing..."))
end_time <- Sys.time()
end_time - start_time
print(glue::glue("Saving Rdata..."))
save(bs.filtered.bsseq, file = "RData/bsseq.RData")
if (length(grep("genomecenter.ucdavis.edu", .libPaths())) > 
    0 & genome == "hg38") {
  dir.create("LOLA")
  setwd("LOLA")
  dmrList <- sigRegions %>% DMRichR::dmrList()
  LOLA <- function(x) {
    dir.create(names(dmrList)[x])
    setwd(names(dmrList)[x])
    dmrList[x] %>% DMRichR::chromHMM(regions = regions, 
                                     cores = floor(cores/3)) %>% DMRichR::chromHMM_heatmap()
    dmrList[x] %>% DMRichR::roadmap(regions = regions, 
                                    cores = floor(cores/3)) %>% DMRichR::roadmap_heatmap()
    if (file.exists("Rplots.pdf")) {
      file.remove("Rplots.pdf")
    }
  }
  parallel::mclapply(seq_along(dmrList), LOLA, mc.cores = 3, 
                     mc.silent = TRUE)
  setwd("..")
}
sigRegions %>% DMRichR::prepareHOMER(regions = regions)
DMRichR::HOMER(genome = genome, cores = cores)
dir.create("Global")
bs.filtered.bsseq %>% DMRichR::globalStats(genome = genome, 
                                           testCovariate = testCovariate, adjustCovariate = adjustCovariate, 
                                           matchCovariate = matchCovariate) %>% openxlsx::write.xlsx("Global/smoothed_globalStats.xlsx")
windows <- bs.filtered.bsseq %>% DMRichR::windows(goi = goi)
CpGs <- bs.filtered.bsseq %>% DMRichR::CpGs()
plots <- c("windows", "CpGs")
if (genome %in% c("hg38", "hg19", "mm10", "mm9", "rheMac10", 
                  "rheMac8", "rn6", "danRer11", "galGal6", "bosTau9", "panTro6", 
                  "dm6", "susScr11", "canFam3")) {
  CGi <- bs.filtered.bsseq %>% DMRichR::CGi(genome = genome)
  plots <- c("windows", "CpGs", "CGi")
}
purrr::walk(plots, function(plotMatrix, group = bs.filtered.bsseq %>% 
                              pData() %>% dplyr::as_tibble() %>% dplyr::pull(!!testCovariate) %>% 
                              forcats::fct_rev()) {
  title <- dplyr::case_when(plotMatrix == "windows" ~ "20Kb Windows", 
                            plotMatrix == "CpGs" ~ "Single CpG", plotMatrix == 
                              "CGi" ~ "CpG Island")
  plotMatrix %>% get() %>% DMRichR::PCA(testCovariate = testCovariate, 
                                        bs.filtered.bsseq = bs.filtered.bsseq) %>% ggplot2::ggsave(glue::glue("Global/{title} PCA.pdf"), 
                                                                                                   plot = ., device = NULL, width = 11, height = 8.5)
  plotMatrix %>% get() %>% DMRichR::densityPlot(group = group) %>% 
    ggplot2::ggsave(glue::glue("Global/{title} Density Plot.pdf"), 
                    plot = ., device = NULL, width = 11, height = 4)
  Glimma::glMDSPlot(plotMatrix %>% get(), groups = cbind(bsseq::sampleNames(bs.filtered.bsseq), 
                                                         pData(bs.filtered.bsseq)) %>% dplyr::as_tibble() %>% 
                      dplyr::select(-col) %>% dplyr::rename(Name = bsseq..sampleNames.bs.filtered.bsseq.), 
                    path = getwd(), folder = "interactiveMDS", html = glue::glue("{title} MDS plot"), 
                    launch = FALSE)
})
sigRegions %>% DMRichR::smoothPheatmap(bs.filtered.bsseq = bs.filtered.bsseq, 
                                       testCovariate = testCovariate)
cat("\n[DMRichR] Performing DMRichments \t\t\t\t", format(Sys.time(), 
                                                          "%d-%m-%Y %X"), "\n")
DMRich <- function(x) {
  if (genome %in% c("hg38", "hg19", "mm10", "mm9", "rheMac10", 
                    "rheMac8", "rn6", "danRer11", "galGal6", "bosTau9", 
                    "panTro6", "dm6", "susScr11", "canFam3")) {
    print(glue::glue("Running CpG annotation enrichments for {names(dmrList)[x]}"))
    dmrList[x] %>% DMRichR::DMRichCpG(regions = regions, 
                                      genome = genome) %T>% openxlsx::write.xlsx(file = glue::glue("DMRichments/{names(dmrList)[x]}_CpG_enrichments.xlsx")) %>% 
      DMRichR::DMRichPlot(type = "CpG") %>% ggplot2::ggsave(glue::glue("DMRichments/{names(dmrList)[x]}_CpG_enrichments.pdf"), 
                                                            plot = ., width = 4, height = 3)
  }
  print(glue::glue("Running gene region annotation enrichments for {names(dmrList)[x]}"))
  dmrList[x] %>% DMRichR::DMRichGenic(regions = regions, 
                                      TxDb = TxDb, annoDb = annoDb) %T>% openxlsx::write.xlsx(file = glue::glue("DMRichments/{names(dmrList)[x]}_genic_enrichments.xlsx")) %>% 
    DMRichR::DMRichPlot(type = "genic") %>% ggplot2::ggsave(glue::glue("DMRichments/{names(dmrList)[x]}_genic_enrichments.pdf"), 
                                                            plot = ., width = 4, height = 4)
}
dmrList <- sigRegions %>% DMRichR::dmrList()
dir.create("DMRichments")
purrr::walk(seq_along(dmrList), DMRich)
purrr::walk(dplyr::case_when(genome %in% c("hg38", "hg19", 
                                           "mm10", "mm9", "rn6") ~ c("CpG", "genic"), TRUE ~ "genic") %>% 
              unique(), function(type) {
                print(glue::glue("Creating DMRichMultiPlots for {type} annotations"))
                DMRichR::DMparseR(direction = c("All DMRs", "Hypermethylated DMRs", 
                                                "Hypomethylated DMRs"), type = type) %>% DMRichR::DMRichPlot(type = type, 
                                                                                                             multi = TRUE) %>% ggplot2::ggsave(glue::glue("DMRichments/{type}_multi_plot.pdf"), 
                                                                                                                                               plot = ., device = NULL, height = dplyr::case_when(type == 
                                                                                                                                                                                                    "genic" ~ 5, type == "CpG" ~ 3.5), width = 7)
              })
cat("\n[DMRichR] Testing for imprinted gene enrichment \t\t", 
    format(Sys.time(), "%d-%m-%Y %X"), "\n")
dmrList <- sigRegions %>% DMRichR::dmrList()
sink("DMRs/human_imprinted_gene_overlaps.txt")
purrr::walk(seq_along(dmrList), function(x) {
  print(glue::glue("Analyzing {names(dmrList)[x]}"))
  dmrList[x] %>% DMRichR::imprintOverlap(regions = regions, 
                                         TxDb = TxDb, annoDb = annoDb)
})
sink()
tryCatch({
  regions %>% DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>% 
    DMRichR::Manhattan()
}, error = function(error_condition) {
  print(glue::glue("Manhattan plot error"))
})
cat("\n[DMRichR] Performing gene ontology analyses \t\t\t", 
    format(Sys.time(), "%d-%m-%Y %X"), "\n")
dir.create("Ontologies")
if (genome %in% c("hg38", "hg19", "mm10", "mm9")) {
  print(glue::glue("Running GREAT"))
  GREATjob <- sigRegions %>% dplyr::as_tibble() %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>% 
    rGREAT::submitGreatJob(bg = regions, species = genome, 
                           rule = "oneClosest", request_interval = 1, version = "4.0.4")
  print(glue::glue("Saving and plotting GREAT results"))
  GREATjob %>% rGREAT::getEnrichmentTables(category = "GO") %T>% 
    openxlsx::write.xlsx(file = glue::glue("Ontologies/GREAT_results.xlsx")) %>% 
    DMRichR::slimGO(tool = "rGREAT", annoDb = annoDb, 
                    plots = FALSE) %T>% openxlsx::write.xlsx(file = glue::glue("Ontologies/GREAT_slimmed_results.xlsx")) %>% 
    DMRichR::GOplot() %>% ggplot2::ggsave(glue::glue("Ontologies/GREAT_plot.pdf"), 
                                          plot = ., device = NULL, height = 8.5, width = 10)
}
if (GOfuncR == TRUE) {
  print(glue::glue("Running GOfuncR"))
  sigRegions %>% DMRichR::GOfuncR(regions = regions, n_randsets = 1000, 
                                  upstream = 5000, downstream = 1000, annoDb = annoDb, 
                                  TxDb = TxDb) %T>% openxlsx::write.xlsx(glue::glue("Ontologies/GOfuncR.xlsx")) %>% 
    DMRichR::slimGO(tool = "GOfuncR", annoDb = annoDb, 
                    plots = FALSE) %T>% openxlsx::write.xlsx(file = glue::glue("Ontologies/GOfuncR_slimmed_results.xlsx")) %>% 
    DMRichR::GOplot() %>% ggplot2::ggsave(glue::glue("Ontologies/GOfuncR_plot.pdf"), 
                                          plot = ., device = NULL, height = 8.5, width = 10)
}
if (genome != "TAIR10" & genome != "TAIR9") {
  tryCatch({
    print(glue::glue("Running enrichR"))
    enrichR:::.onAttach()
    dbs <- c("GO_Biological_Process_2018", "GO_Cellular_Component_2018", 
             "GO_Molecular_Function_2018", "KEGG_2019_Human", 
             "Panther_2016", "Reactome_2016", "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")
    if (genome %in% c("mm10", "mm9", "rn6")) {
      dbs %>% gsub(pattern = "Human", replacement = "Mouse")
    }
    else if (genome %in% c("danRer11", "dm6")) {
      if (genome == "danRer11") {
        enrichR::setEnrichrSite("FishEnrichr")
      }
      else if (genome == "dm6") {
        enrichR::setEnrichrSite("FlyEnrichr")
      }
      dbs <- c("GO_Biological_Process_2018", "GO_Cellular_Component_2018", 
               "GO_Molecular_Function_2018", "KEGG_2019")
    }
    sigRegions %>% DMRichR::annotateRegions(TxDb = TxDb, 
                                            annoDb = annoDb) %>% dplyr::select(geneSymbol) %>% 
      purrr::flatten() %>% enrichR::enrichr(dbs) %>% 
      purrr::set_names(names(.) %>% stringr::str_trunc(31, 
                                                       ellipsis = "")) %T>% openxlsx::write.xlsx(file = glue::glue("Ontologies/enrichr.xlsx")) %>% 
      DMRichR::slimGO(tool = "enrichR", annoDb = annoDb, 
                      plots = FALSE) %T>% openxlsx::write.xlsx(file = glue::glue("Ontologies/enrichr_slimmed_results.xlsx")) %>% 
      DMRichR::GOplot() %>% ggplot2::ggsave(glue::glue("Ontologies/enrichr_plot.pdf"), 
                                            plot = ., device = NULL, height = 8.5, width = 10)
  }, error = function(error_condition) {
    print(glue::glue("Warning: enrichR did not finish. \\\n                      The website may be down or there are internet connection issues."))
  })
}
tryCatch({
  methylLearnOutput <- DMRichR::methylLearn(bs.filtered.bsseq = bs.filtered.bsseq, 
                                            sigRegions = sigRegions, testCovariate = testCovariate, 
                                            TxDb = TxDb, annoDb = annoDb, topPercent = 1, output = "all", 
                                            saveHtmlReport = TRUE)
  if (!dir.exists("./Machine_learning")) {
    dir.create("./Machine_learning")
  }
  if (length(methylLearnOutput) == 1) {
    openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput), 
                         file = "./Machine_learning/Machine_learning_output_one.xlsx")
  }
  else {
    openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput$`Annotated common DMRs`, 
                              RF_Ranking_All_DMRs = methylLearnOutput$`RF ranking`, 
                              SVM_Ranking_All_DMRs = methylLearnOutput$`SVM ranking`), 
                         file = "./Machine_learning/Machine_learning_output_all.xlsx")
  }
  print(glue::glue("Saving RData..."))
  save(methylLearnOutput, file = "RData/machineLearning.RData")
}, error = function(error_condition) {
  print(glue::glue("Warning: methylLearn did not finish. \\\n                      You may have not had enough top DMRs across algrothims."))
})
cat("\n[DMRichR] Summary \t\t\t\t\t", format(Sys.time(), 
                                             "%d-%m-%Y %X"), "\n")
print(glue::glue("Summary: There were {dmrLength} DMRs that covered {sigRegionPercent} of the genome. \\\n                   The DMRs were identified from {backgroundLength } background regions that covered {regionPercent} of the genome.\n                   {tidyHyper} of the DMRs were hypermethylated, and {tidyHypo} were hypomethylated. \\\n                   The methylomes consisted of {tidyCpGs} CpGs.", 
                 dmrLength = sigRegions %>% length() %>% formatC(format = "d", 
                                                                 big.mark = ","), backgroundLength = regions %>% length() %>% 
                   formatC(format = "d", big.mark = ","), tidyHyper = (sum(sigRegions$stat > 
                                                                             0)/length(sigRegions)) %>% scales::percent(), tidyHypo = (sum(sigRegions$stat < 
                                                                                                                                             0)/length(sigRegions)) %>% scales::percent(), tidyCpGs = nrow(bs.filtered) %>% 
                   formatC(format = "d", big.mark = ","), genomeSize = goi %>% 
                   seqinfo() %>% GenomeInfoDb::keepStandardChromosomes() %>% 
                   as.data.frame() %>% purrr::pluck("seqlengths") %>% 
                   sum(), dmrSize = sigRegions %>% dplyr::as_tibble() %>% 
                   purrr::pluck("width") %>% sum(), backgroundSize = regions %>% 
                   dplyr::as_tibble() %>% purrr::pluck("width") %>% 
                   sum(), sigRegionPercent = (dmrSize/genomeSize) %>% 
                   scales::percent(accuracy = 0.01), regionPercent = (backgroundSize/genomeSize) %>% 
                   scales::percent(accuracy = 0.01)))
try(if (sum(blocks$pval < 0.05) > 0 & length(blocks) != 0) {
  print(glue::glue("{length(sigBlocks)} significant blocks of differential methylation \\\n           in {length(blocks)} background blocks"))
}, silent = TRUE)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}
print(glue::glue("Done..."))