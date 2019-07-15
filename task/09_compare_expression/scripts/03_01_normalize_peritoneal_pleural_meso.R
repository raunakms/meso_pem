### LOAD LIBRARIES ----
library("stringr")
library("DESeq")
library("GenABEL")

#### DEFINE PATH ----
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/09_compare_expression")
dir.data <- file.path(dir.wrk, "data")
dir.plot <- file.path(dir.wrk, "plot")
dir.pleu <- file.path("/Data/Raunak/projects/HITnDRIVE/datasets/TCGA_MESO/data/expression/expression_htseqcount/processed_data")
dir.peri <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/expression/main_calls")

#### DEFINE FILES ----
file.pleu <- file.path(dir.pleu, "TCGA-MESO_expression_htseqcount_ensemblid.tsv")
file.peri <- file.path(dir.peri, "MESO_peritoneal_expression_htseqcount_ensemblid.tsv")
file.array <- file.path("/Data/Raunak/softwares/bdvtools/array_process/array_preprocess.R")

### LOAD EXPR DATA 1 ---
dat.pleu <- read.delim(file.pleu, header=T, stringsAsFactors=F, nrow=10)
classes <- sapply(dat.pleu, class)
dat.pleu <- read.delim(file.pleu, header=T, stringsAsFactors=F, colClasses=classes)
colnames(dat.pleu) <- str_replace_all(colnames(dat.pleu), "[.]", "-")
colnames(dat.pleu)[1] <- "ESEMBLE_ID"
dat.pleu$ESEMBLE_ID <- unlist(lapply(str_split(dat.pleu$ESEMBLE_ID, "[.]"), function(x) x[1]))

### LOAD EXPR DATA 2 ---
dat.peri <- read.delim(file.peri, header=T, stringsAsFactors=F)
colnames(dat.peri) <- str_replace_all(colnames(dat.peri), "[.]", "-")
dat.peri <- subset(dat.peri, dat.peri$Type == "protein_coding")

### ANNOTATION ---
dat.annot <- dat.peri[,1:2]

### GET DATA ---
dat.peri <- dat.peri[,-c(2:7)]

#### Subset data by common genes ---
ids <- intersect(dat.pleu$ESEMBLE_ID, dat.peri$ESEMBLE_ID)
dat.pleu <- subset(dat.pleu, dat.pleu$ESEMBLE_ID %in% ids)
dat.peri <- subset(dat.peri, dat.peri$ESEMBLE_ID %in% ids)

dat.pleu <- dat.pleu[match(ids, dat.pleu$ESEMBLE_ID),]
dat.peri <- dat.peri[match(ids, dat.peri$ESEMBLE_ID),]

### MERGE DATA ---
dat <- merge(dat.pleu, dat.peri, by="ESEMBLE_ID")

### Write Output --
file.dat <- file.path(dir.data, "htseqcount_ensemblid_meso_pleural_peritoneal.tsv")
write.table(dat, file.dat, sep="\t", row.names=F, col.names=T, quote=F)

###
cmd <- paste("gzip", file.dat, sep=" ")
system(cmd)

#### LOAD DATA -------------------------------
file.dat <- file.path(dir.data, "htseqcount_ensemblid_meso_pleural_peritoneal.tsv.gz")
dat <- read.delim(file.dat, header=T, stringsAsFactors=F, row.names=1)
colnames(dat) <- str_replace_all(colnames(dat), "[.]", "-")

conds <- factor(rep("Tumor", ncol(dat)))

#### DEseq Normalization ----------------------
cds <- newCountDataSet(dat,conds)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

normalizedCounts <- t( t(counts(cds)) / sizeFactors(cds) )
colnames(normalizedCounts) <- str_replace_all(colnames(normalizedCounts), "[.]", "-")


### Write Output --------------------------------
file.expr <- file.path(dir.data, "htseqcount_ensemblid_meso_pleural_peritoneal_deseqnorm.tsv")
write.table(normalizedCounts, file.expr, sep="\t", row.names=T, col.names=NA, quote=F)

###
cmd <- paste("gzip", file.expr, sep=" ")
system(cmd)

### RE-LOAD DATA ---
file.dat <- file.path(dir.data, "htseqcount_ensemblid_meso_pleural_peritoneal_deseqnorm.tsv.gz")
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
colnames(dat)[1] <- "Gene"

dat.annot <- subset(dat.annot, dat.annot$ESEMBLE_ID %in% dat$Gene)
dat.annot <- dat.annot[match(dat$Gene, dat.annot$ESEMBLE_ID),]
dat$Gene <- dat.annot$Gene_name

### GET EXPRESSION ---
source(file.array)
expr <- getunique.gene.expression(dat)	

dqnorm <- getQuantile.normalize(expr)
gexpr <- log2(dqnorm)
gexpr[gexpr == "-Inf"] <- 0
colnames(gexpr) <- str_replace_all(colnames(gexpr), "[.]", "-")

### write Output ----
file.output <- file.path(dir.data, "expr_combined_meso_pleural_peritoneal.tsv")
write.table(gexpr, file.output, sep="\t", row.names=T, col.names=NA, quote=F)

###
cmd <- paste("gzip", file.output, sep=" ")
system(cmd)

### INVERSE-NORMAL TRANSFORMATION ----
file.array <- file.path("//JBRCSRV009/CollinsGroup/Raunak/softwares/bdvtools/array_process/array_preprocess.R")
dir.w <- file.path("//JBRCSRV009/CollinsGroup/Raunak/projects/MESO_peritoneal/task/09_compare_expression/data")
file.e <- file.path(dir.w, "expr_combined_meso_pleural_peritoneal.tsv.gz")
gexpr <- read.delim(file.e, header=T, stringsAsFactors=F, row.names=1)
colnames(gexpr) <- str_replace_all(colnames(gexpr), "[.]", "-")

source(file.array)
expr.norm <- getInvNormTransform(dat=gexpr)
expr.norm[is.na(expr.norm)] <- 0

### write Output ----
file.output <- file.path(dir.w, "expr_combined_meso_pleural_peritoneal_invnorm.tsv")
write.table(expr.norm, file.output, sep="\t", row.names=T, col.names=NA, quote=F)

###
cmd <- paste("gzip", file.output, sep=" ")
system(cmd)

##### QC ---
genes <- sample(rownames(gexpr), 25)

ids.pleural <- colnames(dat.pleu)[2:ncol(dat.pleu)]
ids.peritoneal <- colnames(dat.peri)[2:ncol(dat.peri)]

pdf(file.path(dir.plot, "histogram.pdf"))
for(i in 1:length(genes)){
	par(mfrow=c(2,1))
	hist(as.numeric(gexpr[genes[i],ids.pleural]), main=genes[i], xlab="Gene Expression", ylab="No. of Samples", cex.axis=0.7, xlim=c(-16,16))

	hist(as.numeric(gexpr[genes[i],ids.peritoneal]), main=genes[i], xlab="Gene Expression", ylab="No. of Samples", cex.axis=0.7, xlim=c(-16,16))
}
dev.off()
