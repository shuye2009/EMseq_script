library(CMplot)

setwd("C:/PROJECTS/Shane/Harding_240124/h4h/result")
pmap <- read.delim("mmbin_CMplot.tab")
CMplot(pmap, LOG10=FALSE, plot.type = "m", file = "pdf", 
       ylab = "Methylation level", multraits = TRUE, multracks = TRUE,
       file.name = "Manhattan_plot", file.output = TRUE,
       bin.size = 5e6, c("#4197d8","#413496", "#f8c120", "#d60b6f"))

