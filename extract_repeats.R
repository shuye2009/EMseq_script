wd <- "C:\\PROJECTS\\Shane\\resource"
setwd(wd)

library(dplyr)
library(ggplot2)

# repeats are downloaded from:
# https://genome.ucsc.edu/cgi-bin/hgTables?hgta_doMainPage=1&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgsid=2162373060_5CI3uk5DfQpNSO8nW0Y8nMw3Ri4f

repeats_df <- read.delim("repeats_hg38_ucsc.tab", header = TRUE)
dim(repeats_df)

bed_df <- repeats_df %>%
  mutate(length = genoEnd - genoStart) %>%
  filter(repClass == "LTR") %>%  #ERV belongs to LTR class
  select(genoName, genoStart, genoEnd, repFamily, length, strand)


dim(bed_df)
head(bed_df)
summary(bed_df)
unique(bed_df$repFamily)

erv_df <- bed_df %>%
  filter(repFamily %in% c("ERVL-MaLR", "ERV1", "ERVK", "ERVL")) 
colnames(erv_df) <- c("chr", "start", "end", "name", "score", "strand")

dim(erv_df)
head(erv_df)
summary(erv_df)
unique(erv_df$name)

write.table(erv_df, "ERV_hg38_ucsc.bed", col.names = FALSE, row.names = FALSE, sep = "\t")

erv_df <- bed_df %>%
  filter(repFamily %in% c("ERVL")) 
colnames(erv_df) <- c("chr", "start", "end", "name", "score", "strand")

dim(erv_df)
head(erv_df)
summary(erv_df)
unique(erv_df$name)

write.table(erv_df, "ERVL_hg38_ucsc.bed", col.names = FALSE, row.names = FALSE, sep = "\t")


## bed for repName 
all_df <- repeats_df %>%
  mutate(length = genoEnd - genoStart) %>%
  select(genoName, genoStart, genoEnd, repName, length, strand)

write.table(all_df, "All_repeats_hg38_ucsc.bed", col.names = FALSE, row.names = FALSE, sep = "\t")

bed_df <- repeats_df %>%
  mutate(length = genoEnd - genoStart) %>%
  filter(repClass == "LTR") %>%  #ERV belongs to LTR class
  select(genoName, genoStart, genoEnd, repName, length, strand)

write.table(bed_df, "LTR_hg38_ucsc.bed", col.names = FALSE, row.names = FALSE, sep = "\t")
dim(bed_df)


dim(bed_df)
head(bed_df)
summary(bed_df)
LTR_name <- unique(bed_df$repName)


for (aName in c("MER4D", "MER21C", "MER57", "MLTA10", "MLT1")){
  NAME <- LTR_name[grepl(aName, LTR_name)]
  print(NAME)
  
  erv_df <- bed_df %>%
    filter(grepl(aName, repName)) 
  colnames(erv_df) <- c("chr", "start", "end", "name", "score", "strand")
  
  dim(erv_df)
  head(erv_df)
  summary(erv_df)
  unique(erv_df$name)
  
  write.table(erv_df, paste0(aName, "_ERV_hg38_ucsc.bed"), col.names = FALSE, row.names = FALSE, sep = "\t")
}


# get qPCR targeted ERVs ####
library(readxl)

qPCR <- readxl::read_excel("ERV primers sequence.xlsx")
ltr_fasta <- read.delim("LTR_hg38.fasta.tab", header = FALSE)
colnames(ltr_fasta) <- c("Name", "Sequence")

matched_df <- NULL

for(i in seq_len(nrow(qPCR))){
  erv <- qPCR[i, 1]
  forward <- qPCR[i, 2]
  reverse <- qPCR[i, 3]
  message(paste(erv, forward, reverse))
  
  for(j in seq_len(nrow(ltr_fasta))){
    name <- ltr_fasta[j, 1]
    sequence <- ltr_fasta[j, 2]
    
    if(grepl(forward, sequence)){
      print(paste(erv, name, "forward"))
      matched_df <- rbind(matched_df, c(erv, name, "forward"))
    }
    if(grepl(reverse, sequence)) {
      print(paste(erv, name, "reverse"))
      matched_df <- rbind(matched_df, c(erv, name, "reverse"))
    }
  }
}

colnames(matched_df) <- c("geneName", "loci", "direction")
write.table(matched_df, "ERV_primers_sequence_matched_loci.tab",
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

matched_bed <- lapply(seq_len(nrow(matched_df)), function(i){
  x <- matched_df[i,]
  loci <- unlist(x[2])
  direction <- unlist(x[3])
  print(loci)
  if(1){
  chr <- unlist(strsplit(loci, split = ":"))[3]
  coor <- unlist(strsplit(loci, split = ":"))[4]
  coors <- as.numeric(unlist(strsplit(coor, split = "-")))
  length <- coors[2] - coors[1]
  strand <- ifelse(direction == "forward", "+", "-")
  return(c(chr, coors, loci, length, strand))
  }
})

matched_bed <- do.call(rbind, matched_bed)

write.table(matched_bed, "ERV_primers_sequence_matched_loci.bed",
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
