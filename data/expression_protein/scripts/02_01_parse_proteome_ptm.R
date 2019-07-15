### DEFINE LIBRARIES ----
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/proteome")
dir.data <- file.path(dir.wrk, "data_ptm")
dir.analysis <- file.path(dir.wrk, "processed_data")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")

### DEFINE FILES ---
file.dat <- file.path(dir.data, "Parsed_withGeneName.tsv")
file.des <- file.path(dir.data, "design_table_proteome.tsv")
file.array <- file.path("/Data/Raunak/softwares/bdvtools/array_process/array_preprocess.R")

### LOAD DATA ---
dat <- read.delim(file.dat, header=TRUE, stringsAsFactors=FALSE)
y <- colnames(dat)[115:141]

### LOAD DESIGN TABLE ---
des <- read.delim(file.des, header=T, stringsAsFactors=F)
des <- des[match(y, des$Column.Names.Tbl),]
ids.normal <- sort(des$SequencingID[which(des$SampleType == "Normal")], decreasing=FALSE)
ids.tumor <- sort(des$SequencingID[which(des$SampleType == "Tumor")], decreasing=FALSE)
ids.celline <- sort(des$SequencingID[which(des$SampleType == "Celline")], decreasing=FALSE)


### FILTER DATA COLUMNS ---
dat <- dat[,c(2,3,4,115:141)]
colnames(dat) <- c("Master","Accession","Gene",des$SequencingID)

### RETAIN CALLS WITH MASTER == "Master Protein" 
dat <- subset(dat, dat$Master == "Master Protein")

### STRATIFY EXPRESSION ---
dat.gene <- subset(dat, select="Gene")
dat.normal <- subset(dat, select=ids.normal)
dat.tumor <- subset(dat, select=ids.tumor)
dat.celline <- subset(dat, select=ids.celline)

### MERGE DATA ----
df <- cbind(dat.gene, dat.normal, dat.tumor, dat.celline)

### PROCESS EXPRESSION  ----
source(file.array)
df.expr <- getunique.gene.expression(df)


### REMOVE NA ---
y.na <- rep(0, nrow(df.expr))
for(i in 1:nrow(df.expr)){
	y.na[i] <- length(which(is.na(df.expr[i,1:23]) == TRUE))
}

del.index <- which(y.na == 23)
expr <- df.expr[-del.index,]

### ORDER GENES ---
expr <- expr[order(rownames(expr), decreasing=FALSE),]

### WRITE OUTPUT ---
file.output <- file.path(dir.analysis, "meso_proteome_ptm_all.tsv")
write.table(expr, file.output, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

cmd <- paste("gzip", file.output, sep=" ")
system(cmd)

### RELOAD DATA ---
file.expr <- file.path(dir.analysis, "meso_proteome_ptm_all.tsv.gz")
expr <- read.delim(file.expr, header=T, stringsAsFactors=F, row.names=1)
colnames(expr) <- str_replace_all(colnames(expr), "[.]", "-")

### SELECT ONLY NORMAL, TUMOR, AND CELL LINES SAMPLES ---
expr <- expr[,-24]


### LOG TRANSFORM ---
expr.log <- log2(expr)

### WRITE OUTPUT ---
file.output <- file.path(dir.analysis, "meso_proteome_ptm_all_log2.tsv")
write.table(expr.log, file.output, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
cmd <- paste("gzip", file.output, sep=" ")
system(cmd)


### QUANTILE NORM ---
dqnorm <- getQuantile.normalize(expr.log)

### WRITE OUTPUT ---
file.output <- file.path(dir.analysis, "meso_proteome_ptm_all_log2_dqnorm.tsv")
write.table(dqnorm, file.output, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
cmd <- paste("gzip", file.output, sep=" ")
system(cmd)


### Z-NORMALIZE ---
znorm <- getZscore(dqnorm)

### WRITE OUTPUT ---
file.output <- file.path(dir.analysis, "meso_proteome_ptm_all_log2_dqnorm_znorm.tsv")
write.table(znorm, file.output, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
cmd <- paste("gzip", file.output, sep=" ")
system(cmd)

