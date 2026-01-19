rm(list=ls())

library(clusterProfiler)
library(bit64)
library(dplyr)
library(ggplot2)
# to plot and transform data
library(readxl)
library(purrr)
library(tibble)
library(tidyr)
library(readr)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(enrichplot)
library(GOSemSim)
library(GenomicRanges)
library(plyranges)
library(rGREAT)

source("C:/PROJECTS/Repositories/Functional_enrichment/GO_KEGG_enrichment_lib.R")

ap_cutoff <- 0.05
sim_cutoff <- 0.7
SemList <- list()
SemList[["BP"]] <- godata('org.Hs.eg.db', ont="BP", computeIC=TRUE)
SemList[["CC"]] <- godata('org.Hs.eg.db', ont="CC", computeIC=TRUE)
SemList[["MF"]] <- godata('org.Hs.eg.db', ont="MF", computeIC=TRUE)

## GSEA and GO analysis ####
gene_list_for_compare_GO <- list()
for (dir in c("Harding_240908", "Harding_240124", "Harding_combined")){
  print(dir)
  hd <- file.path("C:/PROJECTS/Shane", dir, "h4h/cgmaptools")
  setwd(hd)
  
  dataDir <- file.path(hd, "cgmaptools_dmr")
  outDir <- file.path(hd, "GSEA_results")
  if(!dir.exists(outDir)) dir.create(outDir)
  
  dmrDirs <- list.dirs(dataDir, recursive = FALSE)
  
  genome <- "hg38"
  comparisons <- c("Irradiation", "Irradiation+M", "Mutation", "Mutation+I")
  sample_pair<- c("WT-NIR_vs_WT-IR", "R172K-NIR_vs_R172K-IR", "WT-NIR_vs_R172K-NIR", "WT-IR_vs_R172K-IR")
  names(comparisons) <- sample_pair
  
  context <- c("CG", "CHG", "CHH")
  
  for (apair in sample_pair){
    for (cx in context){
      dmrdir <- file.path(dataDir, paste0(apair, "_dmr.", cx), "DMRs")
      gseadir <- file.path(outDir, paste0(apair, "_dmr.", cx))
      if(!dir.exists(gseadir)) dir.create(gseadir)
      
      dmrs <- file.path(dmrdir, "DMRs_annotated.xlsx")
      bgns <- file.path(dmrdir, "background_annotated.xlsx")
      
      annotated_bgn <- readxl::read_xlsx(bgns)
      annotated_dmr <- annotated_bgn %>%
        filter(q.value < 0.05)
      
      openxlsx::write.xlsx(annotated_dmr, file = dmrs)
      
      summarized_dmr <- process_DMR(dmrFile=dmrs) %>%
        arrange(desc(mscore))
      
      hyper <- summarized_dmr %>%
        filter(mscore > 0)
      hypo <- summarized_dmr %>%
        filter(mscore < 0)
      
      summarized_bgn <- process_DMR(dmrFile=bgns)
      
      
      summarized_bgn <- summarized_bgn[!summarized_bgn$geneSymbol %in% summarized_dmr$geneSymbol,] %>%
        arrange(desc(mscore))
      
      gene_list_for_compare_GO[[dir]][[apair]][[cx]][["hyper"]] <- hyper$geneSymbol
      gene_list_for_compare_GO[[dir]][[apair]][[cx]][["hypo"]] <- hypo$geneSymbol
      gene_list_for_compare_GO[[dir]][[apair]][[cx]][["dmr"]] <- union(hyper$geneSymbol, hypo$geneSymbol)
      
      
      rank_table <- rbind(summarized_dmr, summarized_bgn)
      rank_list <- rank_table$mscore
      names(rank_list) <- rank_table$geneSymbol
      
      if(0){
      setwd(gseadir)
        
      ## clusterProfiler analysis  
      run_gseGO_simpleList(gene_list=rank_list, dataName="DMR", adjp_cutoff=0.25, simplify=TRUE, ont="BP")
      run_gseKEGG_simpleList(gene_list=rank_list, dataName="DMR", adjp_cutoff=0.25)
        
      run_enrichGO_simpleList(geneList=hyper$geneSymbol, ont="BP", dataName="hypermethylated", adjp_cutoff=0.05)
      run_enrichKEGG_simpleList(geneList=hyper$geneSymbol, dataName="hypermethylated", adjp_cutoff=0.05)
      
      run_enrichGO_simpleList(geneList=hypo$geneSymbol, ont="BP", dataName="hypomethylated", adjp_cutoff=0.05)
      run_enrichKEGG_simpleList(geneList=hypo$geneSymbol, dataName="hypomethylated", adjp_cutoff=0.05)
      
      run_enrichGO_simpleList(geneList=summarized_dmr$geneSymbol, ont="BP", dataName="DMR", adjp_cutoff=0.05)
      run_enrichKEGG_simpleList(geneList=summarized_dmr$geneSymbol, dataName="DMR", adjp_cutoff=0.05)
      }
      
      if(0){
      ## GOfuncR analysis
      gene_df <- data.frame(gene = c(summarized_dmr$geneSymbol, summarized_bgn$geneSymbol),
                            dmr = c(rep(1, length(summarized_dmr$geneSymbol)), rep(0, length(summarized_bgn$geneSymbol)))) %>%
        drop_na()
      hyper_df <- data.frame(gene = c(hyper$geneSymbol, summarized_bgn$geneSymbol),
                            dmr = c(rep(1, length(hyper$geneSymbol)), rep(0, length(summarized_bgn$geneSymbol)))) %>%
        drop_na()
      hypo_df <- data.frame(gene = c(hypo$geneSymbol, summarized_bgn$geneSymbol),
                             dmr = c(rep(1, length(hypo$geneSymbol)), rep(0, length(summarized_bgn$geneSymbol)))) %>%
        drop_na()
      
      df_list <- list(all=gene_df, hypermethylated=hyper_df, hypomethylated=hypo_df)
      
      lapply(names(df_list), function(aset){
        a_df <- df_list[[aset]]
        go_res <- GOfuncR::go_enrich(a_df, test="hyper")
        
        by(go_res[[1]], go_res[[1]][,'ontology'], head)
        
        res_list <- lapply(c("cellular_component", "molecular_function", "biological_process"), function(ont){
          go_res_ordered <- go_res[[1]] %>%
            filter(ontology == ont) %>%
            arrange(raw_p_overrep) %>%
            filter(FWER_underrep < 0.01 | FWER_overrep < 0.01)
        })
        
        res_df <- bind_rows(res_list)
        go_file <- file.path(gseadir, paste0("GO_enrichment_", aset, ".xlsx"))
        openxlsx::write.xlsx(res_df, file = go_file)
      })
      }
    }
  }
}


## ClusterProfiler analysis ####
for (dir in c("Harding_240908", "Harding_240124", "Harding_combined")){
  print(dir)
  hd <- file.path("C:/PROJECTS/Shane", dir, "h4h/cgmaptools")
  setwd(hd)
  
  dataDir <- file.path(hd, "cgmaptools_dmr")
  outDir <- file.path(hd, "GSEA_results")
  if(!dir.exists(outDir)) dir.create(outDir)
  
  dmrDirs <- list.dirs(dataDir, recursive = FALSE)
  
  genome <- "hg38"
  comparisons <- c("Irradiation", "Irradiation+M", "Mutation", "Mutation+I")
  sample_pair<- c("WT-NIR_vs_WT-IR", "R172K-NIR_vs_R172K-IR", "WT-NIR_vs_R172K-NIR", "WT-IR_vs_R172K-IR")
  names(comparisons) <- sample_pair
  
  context <- c("CG", "CHG", "CHH")
  
  gene_list <- gene_list_for_compare_GO[[dir]]
  
  for (c in context){
    df_list <- list()
    entrez_list <- list()
    for (d in c("hyper", "hypo", "dmr")){  # leave out "dmr"
      d_list <- lapply(sample_pair, function(p){
        gl <- gene_list[[p]][[c]][[d]]
        gl <- bitr(gl, "SYMBOL", "ENTREZID", OrgDb="org.Hs.eg.db")
        gl$ENTREZID
      })
      names(d_list) <- comparisons
      entrez_list[[d]] <- d_list
      dflist_d <- lapply(names(d_list), function(x){
        genes <- d_list[[x]]
        comparison <- rep(x, length(genes))
        direction <- rep(d, length(genes))
        df <- data.frame(gene=genes, comparison=comparison, direction=direction)
      })
      df_list[[d]] <- bind_rows(dflist_d)
    }
    
    
    dmr_df <- bind_rows(df_list)
    
    kegg_res <- compareCluster(gene~comparison+direction, data=dmr_df, fun="enrichKEGG",
                               organism="hsa", pvalueCutoff=0.05)
    
    go_res <- compareCluster(gene~comparison+direction, data=dmr_df, fun="enrichGO",
                             OrgDb='org.Hs.eg.db', ont = "BP")
    kegg_table <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    go_table <- setReadable(go_res, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    
    write.table(head(kegg_table, n=nrow(kegg_table)), file.path(outDir, paste0(c, "_KEGG_table_withDMR.tsv")), 
                row.names = FALSE, col.names=TRUE, sep="\t", quote = FALSE)
    write.table(head(go_table, n=nrow(go_table)), file.path(outDir, paste0(c, "_GOBP_table_withDMR.tsv")), 
                row.names = FALSE, col.names=TRUE, sep="\t", quote = FALSE)
    if(1){
    kegg_hyper <- compareCluster(entrez_list[["hyper"]], fun="enrichKEGG",
                         organism="hsa", pvalueCutoff=0.05)
    kegg_hypo <- compareCluster(entrez_list[["hypo"]], fun="enrichKEGG",
                         organism="hsa", pvalueCutoff=0.05)
    
    go_hyper <- compareCluster(entrez_list[["hyper"]], fun="enrichGO",
                          OrgDb='org.Hs.eg.db', ont = "BP")
    go_hypo <- compareCluster(entrez_list[["hypo"]], fun="enrichGO",
                               OrgDb='org.Hs.eg.db', ont = "BP")
    hypok <- setReadable(kegg_hypo, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    hypog <- setReadable(go_hypo, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    
    }
    
    pdf(file.path(outDir, paste0(c,"_compareCluster_withDMR.pdf")), height = 10, width = 12)
    print(dotplot(kegg_res, x="comparison", font.size=10, title="KEGG", label_format=40) + facet_grid(~direction))
    print(dotplot(kegg_res, x="direction", font.size=10, title="KEGG", label_format=40) + facet_grid(~comparison))
    print(dotplot(go_res, x="comparison", font.size=10, title="GO:BP", label_format=40) + facet_grid(~direction))
    print(dotplot(go_res, x="direction", font.size=10, title="GO:BP", label_format=40) + facet_grid(~comparison))
    
    print(dotplot(go_hyper, title="GO:BP:Hypermethylated"))
    print(dotplot(go_hypo, title="GO:BP:Hypomethylated"))
    print(dotplot(kegg_hyper, title="KEGG:Hypermethylated"))
    print(dotplot(kegg_hypo, title="KEGG:Hypomethylated"))
    dev.off()
  }
}

## rGREAT analysis ####

for (dir in c("Harding_240908", "Harding_240124", "Harding_combined")){
  print(dir)
  hd <- file.path("C:/PROJECTS/Shane", dir, "h4h/cgmaptools")
  setwd(hd)
  
  dataDir <- file.path(hd, "cgmaptools_dmr")
  outDir <- file.path(hd, "GREAT_results")
  if(!dir.exists(outDir)) dir.create(outDir)
  
  dmrDirs <- list.dirs(dataDir, recursive = FALSE)
  
  genome <- "hg38"
  comparisons <- c("Irradiation", "Irradiation+M", "Mutation", "Mutation+I")
  sample_pair<- c("WT-NIR_vs_WT-IR", "R172K-NIR_vs_R172K-IR", "WT-NIR_vs_R172K-NIR", "WT-IR_vs_R172K-IR")
  names(comparisons) <- sample_pair
  
  context <- c("CG", "CHG", "CHH")
  
  for (apair in sample_pair){
    for (cx in context){
      dmrdir <- file.path(dataDir, paste0(apair, "_dmr.", cx), "HOMER")
      greatdir <- file.path(outDir, paste0(apair, "_dmr.", cx))
      if(!dir.exists(greatdir)) dir.create(greatdir)
      
      hyper <- read.delim(file.path(dmrdir, "DMRs_hyper.bed"), header = FALSE)
      hypo <- read.delim(file.path(dmrdir, "DMRs_hypo.bed"), header = FALSE)
      bgns <- read.delim(file.path(dmrdir, "background.bed"), header = FALSE)
      
      colnames(hyper) <- colnames(hypo) <- colnames(bgns) <- c("chr", "start", "end")
      hyper_gr <- makeGRangesFromDataFrame(hyper)
      hypo_gr <- makeGRangesFromDataFrame(hypo)
      bgns_gr <- makeGRangesFromDataFrame(bgns)
      
      bgns_gr <- plyranges::setdiff_ranges(bgns_gr, hyper_gr)
      bgns_gr <- plyranges::setdiff_ranges(bgns_gr, hypo_gr)
      
      dmr_gr <- c(hyper_gr, hypo_gr)
      
      setwd(greatdir)
      
      #gr_list <- list(hyper=hyper_gr, hypo=hypo_gr, background=bgns_gr)
      gr_list <- list(dmr=dmr_gr) 
      
      lapply(names(gr_list), function(n){
        gr <- gr_list[[n]]
        
        for(geneset in c("GO:BP", "GO:MF", "GO:CC", "msigdb:C2:CP:KEGG")){
        
          res <- rGREAT::great(gr, geneset, "GREAT:hg38")
          
          genesetName <- gsub(":", "_", geneset, fixed = TRUE)
          pdf(paste0(cx, "_GREAT_analysis_", genesetName, "_", n, ".pdf"))
          plotRegionGeneAssociations(res)
          plotVolcano(res)
          dev.off()
          
          res_table <- getEnrichmentTable(res) %>%
            dplyr::filter(p_adjust < 0.05)
          great_file <- file.path(greatdir, paste0(cx, "_GREAT_analysis_", genesetName, "_", n, ".xlsx"))
          openxlsx::write.xlsx(res_table, file = great_file)
        }
      })
      
    }
  }
}


process_DMR <- function(dmrFile){
  annotated <- read_xlsx(dmrFile)
  pseudoP <- 1e-100
  selected_df <- annotated %>%
    mutate(multiplier = ifelse(direction == "Hypermethylated", 1, -1)) %>%
    mutate(score = -log(p.value+pseudoP) * multiplier) %>%
    dplyr::select(score, p.value, direction, annotation, geneSymbol)
  
  mixed_gene_df <- selected_df %>%
    group_by(geneSymbol)%>%
    summarise(directionCount = length(unique(direction))) %>%
    filter(directionCount > 1)
  
  
  summarized_df <- selected_df %>%
    group_by(geneSymbol)%>%
    summarise(mscore = mean(score))
  
  nmixed_genes <- nrow(mixed_gene_df)
  all_genes <- nrow(summarized_df)
  
  percent <- round(nmixed_genes*100/all_genes, digits = 2)
  
  message("Portion of mixed genes in ", dmrFile, " is: ", percent, "%")
  
  return(summarized_df)
    
}

