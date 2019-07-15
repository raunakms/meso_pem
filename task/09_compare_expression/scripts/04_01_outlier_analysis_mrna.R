### DEFINE LIBRARIES ----------------------------------------------------------
library("stringr")
library("ggplot2")
library("reshape2")
library("gridExtra")

### DEFINE PATH ----------------------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.analysis <- file.path(dir.wrk, "data/expression/analysis")
dir.task <- file.path(dir.wrk, "task/09_compare_expression")
dir.data <- file.path(dir.task, "data")
dir.output <- file.path(dir.task, "output")
dir.plot <- file.path(dir.task, "plot")

### DEFINE FILES ------------------------------------------------------------
file.expr <- file.path(dir.analysis, "meso_peritoneal_gene_expression_proteincoding.tsv.gz")
file.outlier <- file.path("/Data/Raunak/softwares/bdvtools/git/gesd/runGESD.R")

## Remove Genes with Expr = 0 in > n% of the samples ------
remove.0.matrix <- function(dat, cutoff){
	y <- apply(dat, 1, function(x) length(which(x == 0)))
	del.index <- which(y > (ncol(dat) * cutoff))
	if(length(del.index) != 0){
		dat <- dat[-del.index,]
	} else{
		return(dat)
	}
	return(dat)
}

### Load Expression Data -----------------------------------------------------
expr <- read.delim(file.expr, header=T, row.names=1,  stringsAsFactors=F)
colnames(expr) <- str_replace(colnames(expr), "[.]", "-")
expr <- remove.0.matrix(dat=expr, cutoff=0.25)

### GET OUTLIERS ----
source(file.outlier)
dat.output <- t(apply(expr, 1, function(x) gesd(x, alpha=0.05, value.zscore="NO", r=NA)))
#z <- as.numeric(dat.output[which(dat.output[,1] != 0),1])

dat.select <- subset(dat.output, dat.output[,1] != 0)
dat.outlier <- as.matrix(dat.select[,-1])

dat.out <- melt(t(dat.outlier))
colnames(dat.out) <- c("SampleID","Gene","Value")
dat.out <- subset(dat.out, dat.out$Value != 0)
dat.out <- dat.out[,-3]

dat.out$SampleID <- as.character(dat.out$SampleID)
dat.out$Gene <- as.character(dat.out$Gene)

### ORDER BY SAMPLEID ---
dat.out <- dat.out[order(dat.out$SampleID, decreasing=FALSE),]

### WRITE OUTPUT ---
file.output <- file.path(dir.output, "expr_outlier_genes.tsv")
write.table(dat.out, file.output, sep="\t", row.names=F, col.names=T, quote=F)

### AGGREGATE OUTLIERS ---
ids <- unique(dat.out$SampleID)

df <- data.frame(SampleID=ids)
df$SampleID <- as.character(df$SampleID)
df$Genes <- ""

for(i in 1:nrow(df)){
	dat.temp <- subset(dat.out, dat.out$SampleID == df$SampleID[i])
	df$Genes[i] <- paste(dat.temp$Gene, collapse=":")
}

### WRITE OUTPUT ---
file.output <- file.path(dir.output, "meso_outliers_samplewise_genelist.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)
