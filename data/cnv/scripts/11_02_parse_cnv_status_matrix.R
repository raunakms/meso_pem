### DEFINE LIBRARIES ---
library("stringr")
library("gplots")

### DEFINE PATH ---
dir.wrk <- file.path("/collinsgroup/Raunak/projects/MESO_peritoneal")
dir.data <- file.path(dir.wrk, "data/cnv/seq_call_refseq_genes_meso/")

### DEFINE FILES ---
file.dat <- file.path(dir.data, "cnv_status_pem_nexus_export_transposed.tsv")

### LOAD DATA ---
dat <- read.delim(file.dat, header=TRUE, stringsAsFactors=FALSE)
colnames(dat) <- str_replace_all(colnames(dat), "[.]", "-")

### ADD GENES AS ROWNAMES ---
items <- dat[,1]
genes <- unlist(lapply(str_split(items, "[(]"), function(x) x[1]))
dat$Gene <- genes

#d <- dat[,-1]
#unique(as.character(as.matrix(d)))

### REPLACE DATA ---
dat[dat == ""] <- 0
dat[dat == "CN Gain"] <- 1
dat[dat == "CN Loss"] <- -1
dat[dat == "High Copy Gain"] <- 2
dat[dat == "CN Gain, High Copy Gain"] <- 2
dat[dat == "Homozygous Copy Loss"] <- -2
dat[dat == "CN Loss, Homozygous Copy Loss"] <- -2
dat[dat == "CN Gain, Homozygous Copy Loss"] <- -2
dat[dat == "CN Gain, CN Loss"] <- 0
dat[dat == "CN Loss, High Copy Gain"] <- 2

### ORDER MATRIX ---
dat <- dat[order(dat$Gene, decreasing=FALSE),]
dat <- dat[!duplicated(dat),]

### WRITE OUTPUT ---
file.output <- file.path(dir.data, "cnv_status_pem_matrix.tsv")
write.table(dat, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

### RE-LOAD DATA ---
dat <- read.delim(file.output, header=TRUE, stringsAsFactors=FALSE)
colnames(dat) <- str_replace_all(colnames(dat), "[.]", "-")

d <- dat[,-1]
colnames(d) <- c(1:ncol(d))
z <- apply(d, 1, function(x) paste(x, collapse=":"))

null0 <- paste(rep(0, ncol(d)), collapse=":")
del.index <- which(z == null0)
dat <- dat[-del.index,]

### WRITE OUTPUT ---
file.output <- file.path(dir.data, "cnv_status_pem_matrix.tsv")
write.table(dat, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

### RE-LOAD DATA ---
dat <- read.delim(file.output, header=TRUE, stringsAsFactors=FALSE)
colnames(dat) <- str_replace_all(colnames(dat), "[.]", "-")

z <- table(dat$Gene)
dup.genes <- names(z[which(z > 1)])
dup.index <- which(dat$Gene %in% dup.genes)

list.ugenes <- list()
for(i in 1:length(dup.genes)){
    list.ugenes[[i]] <- which(dat$Gene %in% dup.genes[i])[1]
}
ugenes.index <- unlist(list.ugenes)

del.index <- setdiff(dup.index, ugenes.index)

dat <- dat[-del.index,]

### WRITE OUTPUT ---
file.output <- file.path(dir.data, "cnv_status_pem_matrix.tsv")
write.table(dat, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


### RE-LOAD DATA ---
dat <- read.delim(file.output, header=TRUE, stringsAsFactors=FALSE, row.names=1)
colnames(dat) <- str_replace_all(colnames(dat), "[.]", "-")

### GET SORTED MATRIX ---
mat <- as.matrix(dat)
htm <- heatmap.2(mat)
mat.sorted <- mat[rev(htm$rowInd), htm$colInd]


### WRITE OUTPUT ---
file.output <- file.path(dir.data, "cnv_status_pem_matrix_sorted.tsv")
write.table(mat.sorted, file.output, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
