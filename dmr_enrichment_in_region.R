library(GenomicPlot)
library(bit64)
library(dplyr)
library(ggplot2)



for (dir in c("Harding_240124", "Harding_240908", "Harding_combined")){
  print(dir)
  hd <- file.path("C:/PROJECTS/Shane", dir, "h4h/cgmaptools")
  setwd(hd)
  
  dataDir <- file.path(hd, "cgmaptools_dmr")
  outDir <- file.path(hd, "results")
  if(!dir.exists(outDir)) dir.create(outDir)
  
  dmrDirs <- list.dirs(dataDir, recursive = FALSE)
  
  backgroundBed1 <- c(CpGi="C:/PROJECTS/Shane/resource/cpgIsland_hg38_ucsc.bed")
  backgroundBed2 <- c(ERV="C:/PROJECTS/Shane/resource/ERV_hg38_ucsc.bed")
  backgroundBed3 <- c(Centromere="C:/PROJECTS/Shane/resource/GRCh38_centromere_UCSC.bed")
  backgroundBed4 <- c(ERVqPCR="C:/PROJECTS/Shane/resource/ERV_primers_sequence_matched_loci.bed")
  
  genome <- "hg38"
  comparisons <- c("Irradiation", "Irradiation+M", "Mutation", "Mutation+I")
  names(comparisons) <- c("WT-NIR_vs_WT-IR", "R172K-NIR_vs_R172K-IR", "WT-NIR_vs_R172K-NIR", "WT-IR_vs_R172K-IR")
  
  #targetBed <- c(target="C:/PROJECTS/Shane/Harding_240908/h4h/cgmaptools/cgmaptools_dmr/WT-NIR_vs_R172K-NIR_dmr.CHH/HOMER/DMRs_hyper.bed")
  
  overlap_sig <- function(targetBed, backgroudBed, genome, ignore.strand = TRUE){
    
    bedimportParams <- setImportParams(
      offset = 0, fix_width = 0, fix_point = "center", norm = FALSE,
      useScore = FALSE, outRle = FALSE, useSizeFactor = FALSE, genome = genome
    )
    
    target_gr <- GenomicRanges::reduce(handle_bed(targetBed, 
                                                  importParams = bedimportParams)$query)
    background_gr <- GenomicRanges::reduce(handle_bed(backgroundBed, 
                                                      importParams = bedimportParams)$query)
    
    chromInfo <- circlize::read.chromInfo(species = genome)$df
    genomeSize <- sum(chromInfo$end)
    
    targetSize <- sum(width(target_gr))
    medianSize <- median(width(target_gr))
    backgroundSize <- sum(width(background_gr))
    
    intersect_gr <- GenomicRanges::intersect(target_gr, background_gr, 
                                             ignore.strand = ignore.strand)
    intersectSize <- sum(width(intersect_gr))
    
    #contingency <- matrix(c(intersectSize, targetSize-intersectSize,  backgroundSize-intersectSize, 
    #                      genomeSize+intersectSize-targetSize-backgroundSize), nrow=2, byrow = TRUE,
    #                      dimnames=list(c("target", "nonTarget"), c("background", "nonBackground")))
    #cht <- chisq.test(contingency)
    
    numberPeaks <- length(intersect_gr)
   
    background_p <- backgroundSize/genomeSize
    expected <- targetSize * background_p
    
    log2ratio <- log2((intersectSize+1)/(expected+1))
    
    #logp <- dhyper(intersectSize, backgroundSize, genomeSize-backgroundSize, 
    #               targetSize, log = TRUE) # hypergeometric distribution
    logp <- pbinom(intersectSize, targetSize, background_p, lower.tail=FALSE, log.p = TRUE)
    
    if(intersectSize <= expected){
      logp <- -1*pbinom(intersectSize, targetSize, background_p, lower.tail=TRUE, log.p = TRUE)
    }
    #print(contingency)
    #print(logp)
    return(list(nPeaks=numberPeaks, totalSize=backgroundSize, l2r=log2ratio, logP=logp, 
                targetSize=targetSize, targetWidth = medianSize, 
                intersectSize=intersectSize, expectedSize=expected))
  }
  
  LOGP <- FALSE
  WIDTH <- TRUE
  
  if(LOGP){
    print("plot log-pvalue")
    sink(file.path(outDir, "ERV_enrichment_test.tab"))
    cat(paste0(paste("Comparison", "Context", "Direction", "Annotation",
                     "Number of peaks",	"Total size (bp)", "Log2 Ratio (obs/exp)",
                     "LogP enrichment (+values depleted)", "targetSize", "targetWidth", "intersectSize", "expectedSize", sep="\t"), "\n"))
    backgroundBeds <- c(backgroundBed1, backgroundBed2, backgroundBed3, backgroundBed4)
    for (bname in names(backgroundBeds)){
      backgroundBed <- backgroundBeds[bname]
      enrichment_list <- list()
      for (comp in names(comparisons)){
        for(context in c("CG", "CHG", "CHH")){
          dir <- file.path(dataDir, paste0(comp,"_dmr.", context))
          for (direction in c("hyper", "hypo")){
            #print(paste(comp, context, direction, sep = " : "))
            stat_file <- file.path(dir, "HOMER", paste0("DMRs_", direction, ".bed"))
            names(stat_file) <- paste(comp, context, direction, sep=":")
            enrichment <- overlap_sig(targetBed=stat_file, backgroudBed, genome, 
                                      ignore.strand = TRUE)
            #names(enrichment) <- paste(comp, context, direction, sep=":")
            
            enrichment_list[[comp]][[context]][[direction]] <- enrichment
            cat(paste0(paste(c(comp, context, direction, bname, unlist(enrichment)), 
                             collapse="\t"), "\n"))
          }
        }
      }
      
      enrichment_vec <- unlist(enrichment_list)
      enrichment_vec
    }
    sink()
  }
  
  
  if(WIDTH){
    print("plot dmr width and count")
    sink(file.path(outDir, "DMR_median_width.tab"))
    pdf(file.path(outDir, "DMR_width_boxplot.pdf"), height = 8, width = 8)
    width_df <- NULL
    for (comp in names(comparisons)){
      for(context in c("CG", "CHG", "CHH")){
        dir <- file.path(dataDir, paste0(comp,"_dmr.", context))
        for (direction in c("hyper", "hypo")){
          #print(paste(comp, context, direction, sep = " : "))
          stat_file <- file.path(dir, "HOMER", paste0("DMRs_", direction, ".bed"))
          
          dmr_df <- read.delim(stat_file, header = FALSE) %>%
            dplyr::mutate(width = V3 - V2)
          Width <- dmr_df$width
          medianWidth <- median(Width)
          df <- data.frame(Comparison = rep(comp, length(Width)),
                           Context = rep(context, length(Width)),
                           Direction = rep(direction, length(Width)),
                           Width)
          width_df <- rbind(width_df, df)
          cat(paste0(paste(c(comp, context, direction, medianWidth), 
                           collapse="\t"), "\n"))
        }
      }
    }
    sink()
    
    
    p1 <- ggplot(width_df, aes(x=Context, y=Width, color=Direction)) +
      geom_boxplot() +
      facet_wrap(vars(Comparison)) +
      theme_classic()
    print(p1)
    
    
    p2 <- ggplot(width_df, aes(x=Context, fill=Direction)) +
      geom_bar(position="dodge") +
      scale_y_continuous(name="Count", limits=c(1, 6e+6), transform = "log10") +
      facet_wrap(vars(Comparison)) + 
      theme_classic()
    print(p2)
    
    p3 <- ggplot(width_df, aes(x=Context, fill=Direction)) +
      geom_bar(position="dodge") +
      scale_y_continuous(name="Count", limits=c(0, 4e+6), transform = "identity") +
      facet_wrap(vars(Comparison)) + 
      theme_classic()
    print(p3)
    
    p4 <- ggplot(width_df, aes(x=Context, fill=Direction)) +
      geom_bar(position="dodge") +
      facet_wrap(vars(Comparison)) + 
      theme_classic()
    print(p4)
    
    
    dev.off()
  }
  
}