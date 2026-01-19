library(bsseq)
library(dplyr)
library(tidyr)
library(ggplot2)

# load data
bs <- readRDS("C:/PROJECTS/Shane/Harding_250611/wo_chrY/DSS/bs.RDS")

# get methylation data
meth <- getMeth(bs, type="raw", what="perBase", withDimnames=TRUE)
head(meth)

# get coverage data
coverage <- getCoverage(bs, type="Cov", what="perBase", withDimnames=TRUE)
head(coverage)

pos <- GenomicRanges::granges(bs)

pos_df <- data.frame(
          chr = as.character(seqnames(pos)),
          start = start(pos),
          end = end(pos),
          strand = strand(pos),
          stringsAsFactors = FALSE) %>%
          mutate(CpG = paste0(chr, ":", start, "-", end))

head(pos_df)

rownames(meth) <- pos_df$CpG
rownames(coverage) <- pos_df$CpG

head(meth)
head(coverage)

# filter meth_df by coverage
coverage_threshold <- 25
sample_frac <- 0.75
keep <- rowSums(coverage > coverage_threshold) > as.integer(sample_frac * ncol(coverage))
filtered_meth <- as.data.frame(meth[keep, ])
filtered_cov <- as.data.frame(coverage[keep, ])
# generate a violin plot of filtered meth, with a boxplot on top colored by group

filtered_meth_long <- pivot_longer(filtered_meth, cols = everything(), names_to = "Sample", values_to = "Methylation") %>%
  mutate(group = sapply(strsplit(Sample, "-"), function(x) x[1]))
  
ggplot(filtered_meth_long, aes(x = Sample, y = Methylation)) +
  geom_violin() +
  geom_boxplot(aes(fill = group), width = 0.1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(title = "Filtered Methylation, coverage threshold = 25, sample fraction = 0.75")
ggsave("C:/PROJECTS/Shane/Harding_250611/wo_chrY/DSS/filtered_methylation.pdf", width = 10, height = 10)

# generate a violin plot of filtered coverage, with a boxplot on top colored by group

filtered_cov_long <- pivot_longer(filtered_cov, cols = everything(), names_to = "Sample", values_to = "Coverage") %>%
  mutate(group = sapply(strsplit(Sample, "-"), function(x) x[1]))
  
ggplot(filtered_cov_long, aes(x = Sample, y = Coverage)) +
  geom_violin() +
  geom_boxplot(aes(fill = group), width = 0.1) +
  scale_y_log10() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(title = "Filtered Coverage, coverage threshold = 25, sample fraction = 0.75")
ggsave("C:/PROJECTS/Shane/Harding_250611/wo_chrY/DSS/filtered_coverage.pdf", width = 10, height = 10)

# add rownames as the first column of meth
meth_df <- data.frame(CpG = rownames(filtered_meth), filtered_meth)
head(meth_df)
dim(meth_df)
write.table(meth_df, "C:/PROJECTS/Shane/Harding_250611/wo_chrY/DSS/methylation_data.tsv", row.names = FALSE, sep = "\t", quote = FALSE)

# Check the group assignment
print("Sample names:")
print(colnames(filtered_meth))
print("Group assignment:")
print(table(filtered_meth_long$group, filtered_meth_long$Sample))
