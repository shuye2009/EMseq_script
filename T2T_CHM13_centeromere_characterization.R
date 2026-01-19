rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(readr)
library(tibble)
library(GenomicRanges)
# load centromere file
centromere_file <- "C:\\PROJECTS\\resource\\T2T_CHM13\\chm13v2.0_censat_v2.1.bed"
centromere_data <- read_tsv(centromere_file, col_names = FALSE, col_types = cols(), skip = 1)
head(centromere_data)

# only keep the first 6 columns
centromere_data <- centromere_data[, 1:6]
colnames(centromere_data) <- c("chr", "start", "end", "name", "score", "strand")

# replace strand with * for all rows
centromere_data$strand <- "*"   

# extract satellite names and chromosomes from column 4 like "hor_5_6(S1C1/5/19H1L)", where "hor" is the satellite name, "5" is the chromosome, and "6" is the bin number
sat_names <- sapply(centromere_data[[4]], function(x) strsplit(x, "_")[[1]][1])
sat_chrs <- sapply(centromere_data[[4]], function(x) strsplit(x, "_")[[1]][2])

centromere_data$sat_name <- sat_names

# create a bed file such that each chromosome has a single row by summarizing the centromere data by chromosome using the min and max of the start and end columns

centromere_bed <- centromere_data %>% group_by(chr) %>% summarise(start = min(start), end = max(end))
write_tsv(centromere_bed, "C:\\PROJECTS\\resource\\T2T_CHM13\\chm13v2.0_centromere.bed", col_names = FALSE)

# get average length by sat_names

centromere_data$length <- centromere_data$end - centromere_data$start
centromere_data$chr <- factor(centromere_data$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY"))

# order by chr names using numeric values followed by X and Y and sort by start
centromere_data <- centromere_data[order(centromere_data$chr, centromere_data$start), ]

sat_sum_length <- centromere_data %>% group_by(chr, sat_name) %>% summarise(sum_length = as.integer(sum(length)))
print(sat_sum_length)
# order by chr names using numeric values followed by X and Y
sat_sum_length$chr <- factor(sat_sum_length$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY"))
sat_sum_length <- sat_sum_length[order(sat_sum_length$chr), ]


# plot as stacked bar plot
p <- sat_sum_length %>% ggplot(aes(x = chr, y = sum_length, fill = sat_name)) +
  geom_bar(stat = "identity") +
  labs(title = "Satellite Lengths by Chromosome", x = "Chromosome", y = "Length")
print(p)

# save plot
ggsave("C:\\PROJECTS\\resource\\T2T_CHM13\\chm13v2.0_satellite_lengths.png", width = 10, height = 6)
ggsave("C:\\PROJECTS\\resource\\T2T_CHM13\\chm13v2.0_satellite_lengths.pdf", width = 10, height = 6)

# create a GRanges object from the centromere data
centromere_gr <- GRanges(seqnames = centromere_data$chr, ranges = IRanges(centromere_data$start, centromere_data$end), strand = centromere_data$strand)

mcols(centromere_gr) <- DataFrame(name = centromere_data$name, sat_name = centromere_data$sat_name)
centromere_gr

# get size distribution for each satellite name
# Width is already calculated as 'length' column (end - start)
# Plot histogram/density for each sat_name

# Create histogram with facets for each satellite name
p_size_dist <- centromere_data %>% 
  ggplot(aes(x = length, fill = sat_name)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  scale_x_log10() +
  facet_wrap(~sat_name, scales = "free", ncol = 4) +
  labs(title = "Size Distribution of Satellite Elements by Type",
       x = "Length (bp, log10 scale)",
       y = "Count") +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = 9, face = "bold"))

print(p_size_dist)

ggsave("C:\\PROJECTS\\resource\\T2T_CHM13\\chm13v2.0_satellite_size_distribution.png", 
       p_size_dist, width = 16, height = 12)
ggsave("C:\\PROJECTS\\resource\\T2T_CHM13\\chm13v2.0_satellite_size_distribution.pdf", 
       p_size_dist, width = 16, height = 12)

# Get summary statistics for each satellite name
size_stats <- centromere_data %>% 
  group_by(sat_name) %>% 
  summarise(
    count = n(),
    min_length = min(length),
    q25 = quantile(length, 0.25),
    median_length = median(length),
    mean_length = mean(length),
    q75 = quantile(length, 0.75),
    max_length = max(length),
    total_length = sum(length)
  ) %>%
  arrange(desc(count))

print(size_stats)

# Save summary statistics
write_tsv(size_stats, "C:\\PROJECTS\\resource\\T2T_CHM13\\chm13v2.0_satellite_size_stats.tsv")

# Plot individual elements: name vs length, colored by sat_name, faceted by chr
# Create custom color palette for sat_names with 'hor' as red
all_sat_names <- unique(centromere_data$sat_name)

# Define a palette of 12 highly distinguishable colors
distinct_colors <- c(
  "#E41A1C",  # Red
  "#377EB8",  # Blue
  "#4DAF4A",  # Green
  "#984EA3",  # Purple
  "#FF7F00",  # Orange
  "#FFFF33",  # Yellow
  "#A65628",  # Brown
  "#F781BF",  # Pink
  "#00CED1",  # Turquoise
  "#8B4513",  # Saddle Brown
  "#32CD32",  # Lime Green
  "#700b59"   # Crimson
)

# Assign colors: 'hor' gets red, others get remaining colors
sat_colors <- rep(distinct_colors, length.out = length(all_sat_names))
names(sat_colors) <- all_sat_names

# Ensure 'hor' gets red color if it exists
if("hor" %in% all_sat_names) {
  sat_colors["hor"] <- "#E41A1C"  # Red
  # Redistribute other colors to remaining satellites
  other_sats <- setdiff(all_sat_names, "hor")
  sat_colors[other_sats] <- distinct_colors[2:(length(other_sats) + 1)]
}

# Set name as factor to preserve order as they appear in sorted data
centromere_data$name <- factor(centromere_data$name, levels = unique(centromere_data$name))

p_elements <- centromere_data %>%
  ggplot(aes(xmin=start, xmax=end, ymin=0, ymax=score, fill=sat_name)) +
  geom_rect(alpha = 0.8) +
  scale_fill_manual(values = sat_colors) +
  facet_wrap(~chr, scales = "free_x", ncol = 4) +
  labs(title = "Centromere Satellite by Chromosome",
       x = "Element Name",
       y = "Score",
       fill = "Satellite Name") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 5))

print(p_elements)

ggsave("C:\\PROJECTS\\resource\\T2T_CHM13\\chm13v2.0_elements_by_chr.png", 
       p_elements, width = 16, height = 14)
ggsave("C:\\PROJECTS\\resource\\T2T_CHM13\\chm13v2.0_elements_by_chr.pdf", 
       p_elements, width = 16, height = 14)
