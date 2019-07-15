### SetDirectories ----------------------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/cnv")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")

### Define Files ------------------------------------------------------------
#file.dat <- file.path(dir.data, "meso_cnv_nexus_table_2016_07_06.txt")
#file.des <- file.path(dir.data, "samplelist_table.tsv")
file.dat <- file.path(dir.data, "meso_cnv_calls_nexus.tsv")
#file.gistic <- file.path(dir.data, "meso_gistic_calls.txt")

### Load Libraries ----------------------------------------------------------
library("ggplot2")
library("gridExtra")
library("stringr")
library("foreach")
library("doParallel")
#library("GenomicRanges")

### Load CNA Data --------------------------------------------------------
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
#des <- read.delim(file.des, header=T, stringsAsFactors=F)


### Format Data --------------------------------------------------------
#dat$SampleID <- ""
#for(i in 1:nrow(des)){
#	dat$SampleID[which(dat$Sample == des$SampleName[i])] <- des$SampleID[i]
#}
#dat <- dat[order(dat$SampleID, decreasing=F),]
#dat$Chromosome.Region <- str_replace_all(dat$Chromosome.Region, ",", "")
#dat$Chr <- unlist(lapply(str_split(dat$Chromosome.Region, ":"), function(x) x[1]))
#dat$Start <- unlist(lapply(str_split(unlist(lapply(str_split(dat$Chromosome.Region, ":"), function(x) x[2])), "-"), function(x) as.numeric(x[1])))
#dat$End <- unlist(lapply(str_split(unlist(lapply(str_split(dat$Chromosome.Region, ":"), function(x) x[2])), "-"), function(x) as.numeric(x[2])))
#dat$Cytoband <- unlist(lapply(str_split(dat$Cytoband, " - "), function(x) paste(x, collapse="-")))
#dat$Gene.Symbols <- unlist(lapply(str_split(dat$Gene.Symbols, ", "), function(x) paste(x, collapse=",")))

#dat <- dat[,c(13,14,15,16,4,6,5,7,8,11,12)]
#file.dat <- file.path(dir.data, "meso_cnv_calls_nexus.tsv")
#write.table(dat, file.dat, sep="\t", row.names=F, col.names=T, quote=F)

### Declate Cluster --------------------------------------------------------------
#no_cores <- 25
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

### Get CNA Genes --------------------------------------------------------------
lpar <- foreach(i = 1:nrow(dat), .combine=rbind, .packages="stringr") %dopar%{
					genes <- unique(str_split(dat$Gene.Symbol[i], ",")[[1]])
					list.dat <- data.frame(SampleID=dat$SampleID[i],
									Chr=dat$Chr[i],
									Start=dat$Start[i],
									End=dat$End[i],
									Status=dat$Event[i],
									Length=dat$Length[i],
									Cytoband=dat$Cytoband[i],
									CNV.Overlap.Percentage=dat$CNV.Overlap.Percentage[i],
									Probe.Median=dat$Probe.Median[i],
									Gene=genes)
}

stopCluster(cl)

lpar$SampleID <- as.character(lpar$SampleID)
lpar$Chr <- as.character(lpar$Chr)
lpar$Status <- as.character(lpar$Status)
lpar$Cytoband <- as.character(lpar$Cytoband)
lpar$Gene <- as.character(lpar$Gene)

del.index1 <- which(str_detect(lpar$Gene, "LINC") == TRUE)
del.index2 <- which(str_detect(lpar$Gene, "MIR") == TRUE)
del.index3 <- which(str_detect(lpar$Gene, "orf") == TRUE)
del.index4 <- which(str_detect(lpar$Gene, "LOC") == TRUE)
del.index5 <- which(str_detect(lpar$Gene, "RNU") == TRUE)
del.index6 <- which(str_detect(lpar$Gene, "FAM") == TRUE)
del.index7 <- which(str_detect(lpar$Gene, "^OR") == TRUE)
del.index <- c(del.index1, del.index2, del.index3, del.index4, del.index5, del.index6, del.index7)

lpar <- lpar[-del.index,]


lpar <- subset(lpar, as.character(lpar$SampleID) != "MESO-20")

file.output <- file.path(dir.analysis, "meso_cnv_genes.tsv")
write.table(lpar, file.output, sep="\t", row.names=F, col.names=T, quote=F)

genes <- unique(lpar$Gene)
file.output <- file.path(dir.analysis, "meso_cnv_genelist.tsv")
write.table(genes, file.output, sep="\t", row.names=F, col.names=F, quote=F)


### (1) Stat: CNV --------------------------------------------------------------
dstat1 <- data.frame(SampleID=names(table(dat$SampleID)), 
					FreqCNV=as.numeric(table(dat$SampleID)),
					FreqGene=as.numeric(table(lpar$SampleID)))
dstat1$SampleID <- as.character(dstat1$SampleID)

file.output <- file.path(dir.analysis, "cnv_samplewise_stats.tsv")
write.table(dstat1, file.output, sep="\t", row.names=F, col.names=T, quote=F)


### Plot Graph --------------------------------------------------------------
dstat1 <- dstat1[-which(dstat1$SampleID %in% "MESO-20"),]
dstat1 <- dstat1[order(dstat1$FreqCNV, decreasing=T),]
dstat1$SampleID <- factor(dstat1$SampleID, levels=dstat1$SampleID)

p1 <- ggplot(dstat1, aes(y=FreqCNV, x=SampleID)) +
		geom_bar(stat="Identity") +
		theme(
			axis.text.x = element_text(size = 7, color="black", angle=90, hjust=1, vjust=0.5),
			axis.text.y = element_text(size = 7, color="black"),
			axis.title = element_text(size = 7, color="black"),
			plot.title = element_text(size = 10, color="black"),
			axis.ticks = element_line(size=0.4),
			legend.position = "none") +
		ylab("No. of CNA loci") +
		xlab("") + 
		ggtitle("Distribution of CNA Events")
		
p2 <- ggplot(dstat1, aes(y=FreqGene, x=SampleID)) +
		geom_bar(stat="Identity") +
		theme(
			axis.text.x = element_text(size = 7, color="black", angle=90, hjust=1, vjust=0.5),
			axis.text.y = element_text(size = 7, color="black"),
			axis.title = element_text(size = 7, color="black"),
			plot.title = element_text(size = 10, color="black"),
			axis.ticks = element_line(size=0.4),
			legend.position = "none") +
		ylab("No. of CNA Genes") +
		xlab("") + 
		ggtitle("Distribution of CNA Genes")

file.plot <- file.path(dir.plot, "cnv_samplewise_stat.pdf")
pdf(file.plot, width=7, height=4)		
	grid.arrange(p1, p2, ncol=1, nrow=2)
dev.off()


### Filter by GISTIC ------------------------------------------------------------
#gr.lpar <- makeGRangesFromDataFrame(lpar, keep.extra.columns=TRUE, ignore.strand=TRUE)
#gr.gistic <- makeGRangesFromDataFrame(dat.gistic, keep.extra.columns=TRUE, ignore.strand=TRUE)

#gr.sub <- subsetByOverlaps(gr.lpar, gr.gistic, type="any")

#dat.sub <- as.data.frame(gr.sub)

### Stat: CNA Genes ------------------------------------------------------------
lpar <- subset(lpar, !(lpar$SampleID %in% "MESO-20"))
dstat2 <- data.frame(Gene=names(table(lpar$Gene)), Freq=as.numeric(table(lpar$Gene)))
dstat2 <- dstat2[order(dstat2$Freq, decreasing=T),]

dstat2 <- dstat2[1:100,]
dstat2$Gene <- factor(dstat2$Gene , levels=dstat2$Gene)

p3 <- ggplot(dstat2, aes(y=Freq, x=Gene)) +
		geom_bar(stat="Identity") +
		theme(
			axis.text.x = element_text(size = 7, color="black", angle=90, hjust=1, vjust=0.5),
			axis.text.y = element_text(size = 7, color="black"),
			axis.title = element_text(size = 7, color="black"),
			plot.title = element_text(size = 10, color="black"),
			axis.ticks = element_line(size=0.4),
			legend.position = "none") +
		ylab("No. of Samples") +
		xlab("") + 
		ggtitle("No. of Samples with CNA Genes")

file.plot <- file.path(dir.plot, "cna_freq_stat.pdf")
pdf(file.plot, width=15, height=2.5)		
	grid.arrange(p3, ncol=1, nrow=1)
dev.off()

### Prepare Data for HIT'nDRIVE -------------------------------------------------
#dat <- lpar
#dat <- subset(dat, dat$SampleID != "MESO-20")
#dat$CNA.Status <- dat$Status
#dat$Status <- ifelse((dat$CNA.Status == "CN Gain") | (dat$CNA.Status == "High Copy Gain"), "AMP;", "DEL;")

#df.alt <- dat[,c(1,10,5,11)]

#file.output <- file.path(dir.analysis, "alterations_samplewise_cna.tsv")
#write.table(df.alt, file.output, sep="\t", row.names=F, col.names=T, quote=F)

