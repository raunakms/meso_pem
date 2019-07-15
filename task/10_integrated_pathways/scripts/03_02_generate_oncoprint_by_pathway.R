### Load LIBRARIES -------
library("stringr")
library("ggplot2")
library("gridExtra")

### SetDirectories -----
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.task <- file.path(dir.wrk, "task/10_integrated_pathways")
dir.data <- file.path(dir.task, "data")
dir.analysis <- file.path(dir.task, "output")
dir.plot <- file.path(dir.task, "plots")

### DEFINE FILES ---
file.dat <- file.path(dir.data, "aggregate_alteration_oncoprint_by_pathway_input.tsv")
file.samples <- file.path(dir.data, "order_samples.txt")
file.genes <- file.path(dir.data, "order_genes.txt")
file.script <- file.path("/Data/Raunak/softwares/bdvtools/oncoprint/get_oncoprint_plot.R")

### LOAD DATA ---
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)

### GENERATE ONCOPRINT ---
source(file.script)
file.plot <- file.path(dir.plot, "oncoprint_pathway_final.pdf")  
pdf(file.plot, width=4, height=8)
	oncoprint(dat, orderGenes="costom.sort", file.genes=file.genes, file.samples=file.samples, circularize=FALSE, plot.title.name="", cex.text=5, cex.axis=5, cex.title=5)
dev.off()


