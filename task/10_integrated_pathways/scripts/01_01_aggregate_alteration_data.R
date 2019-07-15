### SetDirectories -----
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.task <- file.path(dir.wrk, "task/10_integrated_pathways")
dir.data <- file.path(dir.task, "data")
dir.analysis <- file.path(dir.task, "output")
dir.plot <- file.path(dir.task, "plots")

### Load LIBRARIES ---------------------------------------------------------------
library("stringr")
library("reshape2")

### DEFINe FILES ---
file.snv <- file.path(dir.wrk, "data/mutation/analysis/alterations_samplewise_mutation.tsv")
file.cnv <- file.path(dir.wrk, "data/cnv/analysis/alterations_samplewise_cna.tsv")
file.fus <- file.path(dir.wrk, "data/genefusion/analysis/alterations_samplewise_fusions.tsv")

### LOAD DATA ---
dat.snv <- read.delim(file.snv, header=T, stringsAsFactors=F)[,1:3]
dat.cnv <- read.delim(file.cnv, header=T, stringsAsFactors=F)[,1:3]
dat.fus <- read.delim(file.fus, header=T, stringsAsFactors=F)[,1:3]

### FILTER ---
dat.snv <- subset(dat.snv, dat.snv$Gene != "")
dat.cnv <- subset(dat.cnv, dat.cnv$Gene != "")
dat.fus <- subset(dat.fus, dat.fus$Gene != "")



### DEFINE MATRIX ---
genes <- unique(c(dat.snv$Gene, dat.cnv$Gene, dat.fus$Gene))
sampleids <- unique(c(dat.snv$SampleID, dat.cnv$SampleID, dat.fus$SampleID))

mat <- matrix("", nrow=length(genes), ncol=length(sampleids), dimnames=list(genes, sampleids))
mat.snv <- mat.cnv  <- mat.fus <- mat

### FUNCTION: FILL MATRIX ---
fill.matrix <- function(mat, dat){
	for(i in 1:nrow(dat)){
		x <- dat$Gene[i]
		y <- dat$SampleID[i]

		mat[x,y] <- dat$Status[i]

		#cat("PROCESSED:", i, "SAMPLEID:", y, "GENE:", x, "STATUS:", dat$Status[i], "\n", sep="\t")
	}
	return(mat)
}

### CALL FUNCTION --
mat.snv <- fill.matrix(mat.snv, dat.snv)
mat.cnv <- fill.matrix(mat.cnv, dat.cnv)
mat.fus <- fill.matrix(mat.fus, dat.fus)

### RESHAPE DATA ---
df.snv <- melt(t(mat.snv))
df.cnv <- melt(t(mat.cnv))
df.fus <- melt(t(mat.fus))
colnames(df.snv) <- colnames(df.cnv) <- colnames(df.fus) <- c("SampleID", "Gene", "Status")

### FUNCTION: GET STATUS --- 
get.status <- function(x){
	y <- c(x[3],x[4],x[5])
	y <- subset(y, y != "")
	result <- paste(y, collapse="")
	return(result)
}

### MEGE RESHAPE DATA ---
df <- data.frame(SampleID=df.snv$SampleID, Gene=df.snv$Gene, Status.SNV=df.snv$Status, Status.CNV=df.cnv$Status, Status.FUS=df.fus$Status)
df$SampleID <- as.character(df$SampleID)
df$Gene <- as.character(df$Gene)
df$Status.SNV <- as.character(df$Status.SNV)
df$Status.CNV <- as.character(df$Status.CNV)
df$Status.FUS <- as.character(df$Status.FUS)

### CALL FUNCTION ---
df$Status <- apply(df, 1, get.status)

### TRIM DATA ---
df <- subset(df, df$Status != "")
df <- df[,c(1,2,6)]


### WRITE OUTPUT ---
file.output <- file.path(dir.data, "aggregate_alteration_oncoprint.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)
