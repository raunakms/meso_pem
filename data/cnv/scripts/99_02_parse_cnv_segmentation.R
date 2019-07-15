### DEFINE LIBRARIES ---
library("stringr")

### DEFINE PATH ----
dir.wrk <- file.path("/collinsgroup/Raunak/projects/MESO_peritoneal/data/cnv")
dir.data <- file.path(dir.wrk, "seq_call_refseq_genes_meso")
dir.script <- file.path(dir.wrk, "scripts")

### GET FILES ---
file.dat <- file.path(dir.data, "meso_pem_cnv_seg_mean_call_matrix.tsv")
#file.refseq <- file.path(dir.data, "refseq_genelist_for_nexus.txt")
file.array <- file.path("/collinsgroup/Raunak/softwares/bdvtools/array_process/array_preprocess.R")

### LOAD GENE TABLE ---
#genes <- read.delim(file.refseq, header=FALSE, stringsAsFactors=FALSE)$V1

### LOAD CNA DATA -----
dat <- read.delim(file.dat, header=TRUE, stringsAsFactors=FALSE, row.names=1)
y <- colnames(dat)
genes <- unlist(lapply(str_split(y, "[.]"), function(x) x[4]))
colnames(dat) <- genes
dat <- t(dat)

### REPLACE ? TO NA ---
#dat[dat == "?"] <- NA

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

### REPLACE <SPACE> TO <BLANK> ----
for(i in 1:nrow(dat)){
	for(j in 1:ncol(dat)){
		dat[i,j] <- str_replace_all(dat[i,j], " ", "")
		cat("PROCESSED", "i =", i, "j =", j, "\n", sep=" ")
	}
}

### REPLACE ? TO NA ---
#dat[dat == "?"] <- NA
dat[dat == ""] <- NA

### REMOVE NA ---
na.count <- NA
for(i in 1:nrow(dat)){
	x <- dat[i,]
	ind <- which(is.na(x))
	na.count[i] <- length(ind)
}
del.ind <- which(na.count >= 17)
dat <- dat[-del.ind,]

### WRITE TABLE ----
file.dat <- file.path(dir.data, "meso_pem_cnv_seg_mean_call_matrix_refined.tsv")
write.table(dat, file.dat, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)


### RELOAD DATA FILE ------
dat <- read.delim(file.dat, header=TRUE, stringsAsFactors=FALSE)
colnames(dat) <- str_replace_all(colnames(dat), "[.]", "-")

### GET UNIQUE GENE VALUE ---
source(file.array)
dat.cnv <- getunique.gene.expression(dat)
dat.cnv[which(is.na(dat.cnv))] <- 0

### WRITE TABLE ----
file.dat <- file.path(dir.data, "meso_pem_cnv_seg_mean_call_matrix_parsed.tsv")
write.table(dat.cnv, file.dat, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
