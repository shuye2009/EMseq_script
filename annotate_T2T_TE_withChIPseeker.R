rm(list = ls())

library(ChIPseeker)
library(org.Hs.eg.db)

wd <- "C:/PROJECTS/resource/T2T_CHM13/"
setwd(wd)

# Create TxDb from GFF3
txdb <- makeTxDbFromGFF("Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf",
                        format = "gtf",
                        organism = "Homo sapiens")

TE_bed <- "chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14_dup.bed"
names(TE_bed) <- "T2T_CHM13_rmsk_TE"

# ChIPseeker annotation

# Read and annotate peaks
peaks <- readPeakFile(TE_bed)
peakAnno <- annotatePeak(peaks, 
                         tssRegion=c(-1000, 0),
                         TxDb=txdb, 
                         annoDb="org.Hs.eg.db")

# View annotation summary
print(peakAnno)

# Plot annotation distribution
pdf("T2T_CHM13_rmsk_TE_chipseeker_annotation.pdf", width=10, height=8)
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
plotDistToTSS(peakAnno, title="Distribution of peaks relative to TSS")
dev.off()

# Save annotated peaks to file
annotated_peaks <- as.data.frame(peakAnno)
write.table(annotated_peaks,    
            file="T2T_CHM13_rmsk_TE_annotated_ChIPseeker.tsv",
            sep="\t", 
            row.names=FALSE, 
            quote=FALSE)


# Diagnostic: Check which chromosomes are in peaks vs annotated

standard_chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")

# Check for peaks that might have been filtered
cat("\nChromosomes in original peaks:\n")
table(seqnames(peaks))[standard_chromosomes]

cat("\nChromosomes in annotated peaks:\n")
table(annotated_peaks$seqnames)[standard_chromosomes]

# Check for peaks that might have been filtered
cat("\nNumber of peaks:", length(peaks), "\n")
cat("Number of annotated peaks:", nrow(annotated_peaks), "\n")
cat("Difference:", length(peaks) - nrow(annotated_peaks), "\n")
