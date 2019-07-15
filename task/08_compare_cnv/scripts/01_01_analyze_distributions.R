### LOAD LIBRAIES ---
library("stringr")
library("ggplot2")
library("gridExtra")
library("bigmemory")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/08_compare_cnv")
dir.data <- file.path(dir.wrk, "data")
dir.output <- file.path(dir.wrk, "output")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")
dir.ov <- file.path("/Data/Raunak/HITnDRIVE/datasets/TCGA_OV/data/cnv/seq_call_refseq_genes_ov")

### DEFINE FILES ---
file.meso <- file.path(dir.data, "meso_cnv_seg_values_calls_parsed.tsv")
file.ov <- file.path(dir.ov, "tcga_ov_seg_values_calls_with_ids_unique.tsv")
file.sampleid <- file.path(dir.ov, "sample_ids.txt")
file.script <- file.path("/Data/Raunak/softwares/bdvtools/array_process/diff_expr.R")

## Remove Genes with Expr = NA in > 25% of the samples ------
remove.na.matrix <- function(dat, cutoff){
	y <- apply(dat, 1, function(x) length(which(is.na(x))))
	del.index <- which(y > (ncol(dat) * cutoff))
	if(length(del.index) != 0){
		dat <- dat[-del.index,]
	} else{
		return(dat)
	}
	return(dat)
}

### LOAD DATA ---
get.distribution <- function(file.dat, dir.plot, cancer_type){
	require("stringr")
	dat <-  read.delim(file.dat, header=T, stringsAsFactors=F, row.names=1)
	colnames(dat) <- str_replace_all(colnames(dat), "[.]", "-")

	genes <- c(sample(rownames(dat), 25))

	file.plot <- file.path(dir.plot, paste("cnv_distribution_", cancer_type, ".pdf", sep=""))
	pdf(file.plot, height=4, width=4)
	for(i in 1:length(genes)){
		gene <- genes[i]
		h <- hist(as.numeric(dat[gene,]), plot=FALSE) 
		plot(h, 
			main=gene, xlim=c(-2,2),
		 	xlab="CNV segment median", ylab="No. of Samples",
		 	cex.main=0.8, cex.lab=0.7, cex.axis=0.7, las=1, tck=-0.03)
	}
	dev.off()
}


#get.distribution(file.dat=file.ov, dir.plot, cancer_type="OV")



### LOAD MESO CNV ----
dat.meso <-  read.delim(file.meso, header=T, stringsAsFactors=F, row.names=1)
colnames(dat.meso) <- str_replace_all(colnames(dat.meso), "[.]", "-")
dat.meso <- remove.na.matrix(dat=dat.meso, cutoff=0.25)

### LOAD OV CNV ----
sampleid <-  read.delim(file.sampleid, header=F, stringsAsFactors=F)$V1
dat.ov <- read.big.matrix(filename=file.ov, sep="\t", header=TRUE, 
				col.names=sampleid, has.row.names=TRUE, ignore.row.names=FALSE,
				type="double", skip=0)

dat.ov <- as.matrix(dat.ov)
dat.ov <- remove.na.matrix(dat=dat.ov, cutoff=0.25)

#shapirotest for nomarality
get.shapirotest <- function(dat){

	func.shapirotest <- function(x){
		st <- shapiro.test(x)
		return(st)
	}

	shapirotest.list <- apply(as.matrix(dat), 1, function(x) func.shapirotest(x))
	pvalue <- do.call(c,lapply(shapirotest.list, function(x) x$p.value))
	statistic <- do.call(c,lapply(shapirotest.list, function(x) x$statistic))
	fdr <- p.adjust(pvalue, method = "BH", n = length(pvalue))

	dat.summary <- data.frame(Gene=rownames(dat), 
								statistic=statistic,
								pvalue=pvalue, 
								fdr=fdr)	
	rownames(dat.summary) <- c(1:nrow(dat))
	dat.summary$Gene <- as.character(dat.summary$Gene)
	dat.summary <- dat.summary[order(dat.summary$pvalue, decreasing=F),]

	dat.summary$NormalityStatus <- ifelse(dat.summary$fdr > 0.1, 1, 0)	
	return(dat.summary)
}

dat.summary.ov <- get.shapirotest(dat.ov)

file.output <- file.path(dir.output, "shapirotest_results_OV.tsv")
write.table(dat.summary.ov, file.output, sep="\t", row.names=F, col.names=T, quote=F)
