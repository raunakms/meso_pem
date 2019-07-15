### SetDirectories -----
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.task <- file.path(dir.wrk, "task/10_integrated_pathways")
dir.data <- file.path(dir.task, "data")
dir.analysis <- file.path(dir.task, "output")
dir.plot <- file.path(dir.task, "plots")

### Load Libraries ----------------------------------------------------------
library("stringr")
library("ggplot2")
library("gridExtra")

### DEFINe FILES ---
file.dat <- file.path(dir.data, "aggregate_alteration_oncoprint.tsv")
file.geneset <- file.path(dir.data, "special_geneset.tsv")
file.script <- file.path("/Data/Raunak/softwares/bdvtools/oncoprint/get_oncoprint_plot.R")

### LOAD DATA ---
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
geneset <- read.delim(file.geneset, header=T, stringsAsFactors=F)


### GET ONCOPRINT ---
source(file.script)
list.plot <- list()
ctr <- 1
for(i in 1:nrow(geneset)){
	if(i == 9) next
	genes <- str_split(geneset$Genesets[i], ",")[[1]]
	dat.geneset <- subset(dat, dat$Gene %in% genes)

	if(nrow(dat.geneset) == 0) next

	list.plot[[ctr]] <- oncoprint(dat.geneset, orderGenes="memo.sort", file.genes=NA, file.samples=NA, circularize=FALSE, plot.title.name=geneset$Category[i], cex.text=5, cex.axis=5, cex.title=5)
	ctr <- ctr + 1
}


### PLOT ONCOPRINT -----------------------------------------------------------------------
file.plot <- file.path(dir.plot, "oncoprint_special_pathways.pdf")  
pdf(file.plot, width=3, height=6)
for(ctr in 1:length(list.plot)){
	grid.arrange(list.plot[[ctr]], ncol=1)
}	
dev.off()


### FOR LARGE SET --
i <- 9
genes <- str_split(geneset$Genesets[i], ",")[[1]]
dat.geneset <- subset(dat, dat$Gene %in% genes)

file.plot <- file.path(dir.plot, "oncoprint_special_pathways_immune_infiltration.pdf")  
pdf(file.plot, width=3, height=35)
oncoprint(dat.geneset, orderGenes="memo.sort", file.genes=NA, file.samples=NA, circularize=FALSE, plot.title.name=geneset$Category[i], cex.text=5, cex.axis=5, cex.title=5)
dev.off()
