### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/04_driver_stats")
dir.output <- file.path(dir.wrk, "output")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")
dir.enrichment <- file.path(dir.wrk, "enrichment")

### DEFINE LIBRARIES ----
library("stringr")
library("ggplot2")

### DEFINE FILES ---
file.dat <- file.path(dir.enrichment, "driver_genes/enrichment_h.all.v5.0.symbols.txt")

### LOAD DATA ---
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
dat$Category <- lapply(str_split(dat$Category, "_"), function(x) paste(x[2:length(x)], collapse=" "))
dat$nlogfdr <- -log(dat$fdr)

dat <- dat[order(dat$nlogfdr, decreasing=F),]

### Enrichment Plot : GGPLOT2
dat$Category <- factor(dat$Category, levels=dat$Category)

file.plot <- file.path(dir.plot, "enrichment_drivers_hallmarks.pdf")
pdf(file.plot, height=2, width=3)
ggplot(dat, aes(x=Category, y=nlogfdr)) + 
	geom_bar(fill="grey30", stat="identity") +
	coord_flip() +
	theme(
		axis.text = element_text(size = 6, color="black"),
		axis.title = element_text(size = 6, color="black"),
		plot.title = element_text(size = 5, color="black"),
		axis.ticks.y =element_blank(),
		panel.background = element_rect(fill="grey90"),
		legend.position="none") +
		xlab("") +
		ylab("Enrichment Score [-log10(pvalue)]") + ggtitle("") 
dev.off()	
