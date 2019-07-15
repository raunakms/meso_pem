### Load Libraries -----
library("stringr")

### DEFINE PATH ----
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/cnv")
dir.data <- file.path(dir.wrk, "mes_wdpm")

### Define Files -----
file.dat <- file.path(dir.data, "meso_wdpm_cnv_seg_values_calls.tsv")
file.refseq <- file.path(dir.wrk, "seq_call_refseq_genes_meso/refseq_gene_table.tsv")
file.array <- file.path("/Data/Raunak/softwares/bdvtools/array_process/array_preprocess.R")

### Load Gene Table -----
des <- read.delim(file.refseq, header=T, stringsAsFactors=F)

### Load CNA Data -----
dat <- read.delim(file.dat, header=T, stringsAsFactors=F, row.names=1)
colnames(dat) <- unlist(lapply(str_split(colnames(dat), "[.]"), function(x) x[4]))
dat <- t(dat)

### Write Table ----
file.dat <- file.path(dir.data, "meso_wdpm_cnv_seg_values_calls_transformed.tsv")
write.table(dat, file.dat, sep="\t", row.names=T, col.names=NA, quote=F)


### Reload data file ------
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
colnames(dat)[1] <- "Gene"
dat <- na.omit(dat)

### TRIM VALUES WTIH " " <SPACE> ---
for(i in 1:nrow(dat)){
	x <- dat[i,]
	ind <- which(str_detect(x, " ") == TRUE)

	if(length(ind) != 0){
		y <- unlist(lapply(str_split(x, " "), function(x) x[which(x != "")]))
		dat[i,] <- y
	}
	cat("PROCESSED", i, "OF", nrow(dat), "\n", sep="\t")
}


### TRIM VALUES WTIH , <COMMA> ---
for(i in 1:nrow(dat)){
	x <- dat[i,]
	ind <- which(str_detect(x, ",") == TRUE)

	if(length(ind) != 0){
		y <- unlist(lapply(str_split(x, ", "), function(x) x[1]))
		dat[i,] <- y
	}
	cat("PROCESSED", i, "OF", nrow(dat), "\n", sep="\t")
}


### Write Table ----
file.dat <- file.path(dir.data, "meso_wdpm_cnv_seg_values_calls_transformed2.tsv")
write.table(dat, file.dat, sep="\t", row.names=F, col.names=T, quote=F)


### Reload data file ------
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)


### Process GeneExpression ---
source(file.array)
dat1 <- getunique.gene.expression(dat)

### Write Table -----------------------------------------------------------
file.dat <- file.path(dir.data, "meso_wdpm_cnv_seg_values_calls_parsed.tsv")
write.table(dat1, file.dat, sep="\t", row.names=T, col.names=NA, quote=F)



### Reload data file -------------------------------------------------------
dat <- read.delim(file.dat, header=T, stringsAsFactors=F, row.names=1)


