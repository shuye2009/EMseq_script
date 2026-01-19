# two factor model analysis of dmr #

library(DSS)
library(bsseq)
library(comethyl)
library(dplyr)
library(BiocFileCache)
library(GenomicPlot)
library(annotatr)
library(methylKit)
library(ComBatMet)

# fetch gtf file ####
genome <- "hg38"
bfc <- BiocFileCache()
if(genome == "hg38"){
  url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.primary_assembly.annotation.gtf.gz"
}else if(genome == "hg38"){
  url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
}else{
  stop("Genome is not supported!")
}
gz <- bfcrpath(bfc, url)
gtfFile <- gsub(".gz", "", gz)
if(!file.exists(gtfFile)) gtfFile <- R.utils::gunzip(gz, remove=FALSE)

# load bs ####
bsdir <- "C:\\PROJECTS\\Shane\\Harding_combined\\local\\DMRcate_res"

bs <- readRDS(file.path(bsdir, "bss.RDS"))
colnames(bs) <- sapply((strsplit(colnames(bs), split=".", fixed = T)), function(x)x[1])


standardChr <- c(paste0("chr", seq.int(22)), "chrX", "chrY")
bsf <- bs[seqnames(bs) %in% standardChr, ] # remove non-standard chromosomes
seqlevels(bsf) <- standardChr

rm(bs)
gc()


## with batch correction
if(0){
design <- as.data.frame(Reduce(rbind, strsplit(colnames(bsf), split="-")))
colnames(design) <- c("pheno", "treat", "rep")
rownames(design) <- colnames(bsf)
design$batch <- rep(c(1,2,2), by=4)

design <- design |>
  mutate(treat=factor(treat, levels=c("NIR", "IR"))) |>
  mutate(pheno=factor(pheno, levels= c("WT", "R172K"))) |>
  mutate(batch=factor(batch, levels=c(1,2))) |>
  mutate(group=factor(paste(pheno, treat, sep="-"))) |>
  mutate(group=relevel(group, "WT-NIR"))

group <- design$group
batch <- design$batch


beta_mat <-  bsseq::getMeth(bsf, type="smooth")


# batch correction ####
# Adjust for batch effects including biological conditions
adj_beta_mat <- ComBat_met(vmat = beta_mat, 
                           dtype = "b-value",
                              batch = batch, 
                              group = group, 
                              covar_mod = NULL,
                              full_mod = TRUE,
                              shrink = FALSE, 
                              mean.only = FALSE,
                              feature.subset.n = NULL
)
saveRDS(adj_beta_mat, "adj_beta_mat.RDS")
head(adj_beta_mat)
# adjust counts of methylated C
bsf@assays@data$M <- matrix(as.integer(round(bsf@assays@data$Cov * adj_beta_mat)),
                          nrow=nrow(adj_beta_mat), 
                          ncol=ncol(adj_beta_mat), 
                          byrow=FALSE,
                          dimnames = dimnames(adj_beta_mat)
                     )
}

## without batch correction, only use Aug2024 batch 
bsf <- bsf[!grepl("_0", colnames(bsf)), ]
design <- as.data.frame(Reduce(rbind, strsplit(colnames(bsf), split="-")))
colnames(design) <- c("pheno", "treat", "rep")
rownames(design) <- colnames(bsf)
design$batch <- rep(c(1,2,2), by=4)

design <- design |>
  mutate(treat=factor(treat, levels=c("NIR", "IR"))) |>
  mutate(pheno=factor(pheno, levels= c("WT", "R172K"))) |>
  mutate(batch=factor(batch, levels=c(1,2))) |>
  mutate(group=factor(paste(pheno, treat, sep="-"))) |>
  mutate(group=relevel(group, "WT-NIR"))

# call dmr ####
comparisons <- list(treatWT=c("WT-IR-1", "WT-IR-2", "WT-NIR-1", "WT-NIR-2"),
                 phenoNIR=c("R172K-NIR-1", "R172K-NIR-2", "WT-NIR-1", "WT-NIR-2"))

for(aname in names(comparisons)){
  wd <- paste0("C:\\PROJECTS\\Shane\\Harding_combined\\local\\DSS_res_", aname)
  if(!dir.exists(wd)) dir.create(wd)
  setwd(wd)
  

  if(!file.exists(paste0(aname,"_DML.RDS"))){
    DML= DMLtest(bs[, comparisons[[aname]]], 
                 group1=comparisons[[aname]][3:4],
                 group2=comparisons[[aname]][1:2], 
                 smoothing=TRUE)
    saveRDS(DML, paste0(aname,"_DML.RDS"))
  }else{
    DML <- readRDS(paste0(aname,"_DML.RDS"))
  }
  
  if(!file.exists(paste0(aname,"_DMLs.RDS"))){
    DMLs = callDML(DML, p.threshold = 0.001, delta = 0.1)
    DMLs <- na.omit(DMLs)
    saveRDS(DMLs, paste0(aname,"_DMLs.RDS"))
  }else{
    DMLs <- readRDS(paste0(aname,"_DMLs.RDS"))
  }
  
  if(!file.exists(paste0(aname,"_DMRs.RDS"))){
    dmrs <- callDMR(DML, p.threshold = 0.001, delta = 0.1, pct.sig = 0.5) |>
      filter(abs(diff.Methy) > 0.1)
    saveRDS(dmrs, paste0(aname,"_DMRs.RDS"))
  }else{
    dmrs <- readRDS(paste0(aname,"_DMRs.RDS"))
  }
  
  dmrBed <- dmrs
  dmrBed$name <- rownames(dmrBed)
  dmrBed$strand <- rep("*", nrow(dmrBed))
  dmrBed <- dmrBed[, c("chr", "start", "end", "name", "areaStat", "strand")]
  
  write.table(dmrBed, paste0("DMR_", aname, ".bed"), row.names = FALSE, 
              col.names = FALSE, sep="\t", quote=FALSE)
  
  
  write.table(dmrs, paste0("DMR_", aname, ".tsv"), row.names = TRUE, 
              col.names = NA, sep="\t", quote=FALSE)
  
  pdf("dmr_examples.pdf")
  for(i in 1:10)(
    showOneDMR(dmrs[i,], bs[, comparisons[[aname]]])
  )
  dev.off()
}

# annotate dmr ####
annot <- annotatr::builtin_annotations()[grep("hg38", annotatr::builtin_annotations())]
annotations <- build_annotations(genome = genome, annotations = annot)

source("C:/PROJECTS/Shane/script/GO_KEGG_enrichment.R")
for(aname in names(comparisons)){
  wd <- paste0("C:\\PROJECTS\\Shane\\Harding_combined\\local\\DSS_res_", aname)
  if(!dir.exists(wd)) dir.create(wd)
  setwd(wd)
  
  pdf(paste0("DMR_annotations_", aname, ".pdf"))
  
  dm_file = paste0("DMR_", aname, ".bed")
 
  dm_regions = read_regions(con = dm_file, genome = genome, rename_name = 'DM_status', rename_score = 'areaStat')
  dm_regions$DM_status <- ifelse(dm_regions$areaStat>0, "hyper", "hypo")

  table(dm_regions$DM_status)

  dm_random_regions = randomize_regions(
    regions = dm_regions,
    allow.overlaps = TRUE,
    per.chromosome = TRUE)
  
  # Intersect the regions we read in with the annotations
  dm_annotated = annotate_regions(
    regions = dm_regions,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
  # A GRanges object is returned
  print(dm_annotated)
  
  # Coerce to a data.frame
  df_dm_annotated = data.frame(dm_annotated)
  
  # Perform enrichment analysis
  directions <- c("hyper", "hypo")
  
  dm_genes <- lapply(directions, function(x){
    dm_gene <- df_dm_annotated |>
      filter(DM_status == x) |>
      filter(!is.na(annot.symbol)) |>    
      pull(annot.symbol) |>
      unique()
  })
  names(dm_genes) <- directions
  
  gene_list <- list()
  gene_list[[aname]] <- dm_genes
  compareCluster_enrichment(gene_list, names(gene_list), c="DSS", height = 10, width = 6)
  
  # See the GRanges column of dm_annotaed expanded
  print(head(df_dm_annotated, n=50))
  
  dm_random_annotated = annotate_regions(
    regions = dm_random_regions,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = TRUE)
  
  # Find the number of regions per annotation type
  
  dm_annsum_rnd = summarize_annotations(
    annotated_regions = dm_annotated,
    annotated_random = dm_random_annotated,
    quiet = TRUE)
  print(dm_annsum_rnd)
  
  # View the number of regions per annotation and include the annotation
  # of randomized regions
  annots_order = c(
    'hg38_cpg_islands',
    'hg38_cpg_shores',
    'hg38_cpg_shelves',
    'hg38_cpg_inter',
    'hg38_genes_1to5kb',
    'hg38_genes_promoters',
    'hg38_genes_5UTRs',
    'hg38_genes_exons',
    'hg38_genes_introns',
    'hg38_genes_3UTRs',
    'hg38_genes_intergenic',
    'hg38_lncrna_gencode')
  
  dm_vs_kg_annotations_wrandom = plot_annotation(
    annotated_regions = dm_annotated,
    annotated_random = dm_random_annotated,
    annotation_order = annots_order,
    plot_title = 'Dist. of Sites Tested for DM (with rndm.)',
    x_label = 'Annotations',
    y_label = 'Count')
  print(dm_vs_kg_annotations_wrandom)
  
  # View the counts of CpG annotations in data classes
  
  # The orders for the x-axis labels. This is also a subset
  # of the labels (hyper, hypo, none).
  x_order = c(
    'hyper',
    'hypo')
  # The orders for the fill labels. Can also use this
  # parameter to subset annotation types to fill.
  fill_order = c(
    'hg38_cpg_islands',
    'hg38_cpg_shores',
    'hg38_cpg_shelves',
    'hg38_cpg_inter')
  # Make a barplot of the data class where each bar
  # is composed of the counts of CpG annotations.
  dm_vs_cpg_cat = plot_categorical(
    annotated_regions = dm_annotated, 
    annotated_random = dm_random_annotated,
    x='DM_status', fill='annot.type',
    x_order = x_order, fill_order = fill_order, position='stack',
    plot_title = 'DM Status by CpG Annotation Counts',
    legend_title = 'Annotations',
    x_label = 'DM status',
    y_label = 'Count')
  print(dm_vs_cpg_cat)
  
  # View the proportions of data classes in knownGene annotations
  
  # The orders for the x-axis labels.
  x_order = c(
    'hg38_genes_1to5kb',
    'hg38_genes_promoters',
    'hg38_genes_5UTRs',
    'hg38_genes_exons',
    'hg38_genes_introns',
    'hg38_genes_3UTRs',
    'hg38_genes_intergenic')
  # The orders for the fill labels.
  fill_order = c(
    'hyper',
    'hypo')
  dm_vs_kg_cat = plot_categorical(
    annotated_regions = dm_annotated, x='annot.type', fill='DM_status',
    x_order = x_order, fill_order = fill_order, position='stack',
    legend_title = 'DM Status',
    x_label = 'knownGene Annotations',
    y_label = 'Count')
  print(dm_vs_kg_cat)
  
  dev.off()
}



