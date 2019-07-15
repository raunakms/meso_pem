### LOAD LIBRARIES ---
library("stringr")

### DEFINE PATH ----
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/cnv")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")

### DEFINE FILES ---
file.dat <- file.path(dir.data, "meso_cnv_calls_nexus.tsv")

### Load CNA DATA ----
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)

### RETAIN ONLY CALLS MAPPED TO chr3p21 ---
dat <- subset(dat, dat$Chr == "chr3")
dat <- dat[which(str_detect(dat$Cytoband, "p21") == TRUE),]
dat <- subset(dat, dat$SampleID != "MESO-20")

### PREPARE DATA ---
ids <- unique(dat$SampleID)
ctr <- 1
list.dat <- list()
for(i in 1:length(ids)){
	dat.temp1 <- subset(dat, dat$SampleID == ids[i])
	event <- unique(dat.temp1$Event)

	for(j in 1:length(event)){
		dat.temp2 <- subset(dat.temp1, dat.temp1$Event == event[j])

		g <- paste(unique(unlist(str_split(dat.temp2$Gene.Symbols, ","))), collapse=":")
		list.dat[[ctr]] <- data.frame(SampleID=ids[i], Event=event[j], Genes=g)
		ctr <- ctr + 1
	}
}

dat2 <- do.call(rbind.data.frame, list.dat)
dat2$SampleID <- as.character(dat2$SampleID)
dat2$Event <- as.character(dat2$Event)
dat2$Genes <- as.character(dat2$Genes)

dat2$Status <- ifelse( dat2$Event == "CN Gain", 1, ifelse(dat2$Event == "CN Loss", -1, -2) )

### GET GLOBAL GENE FREQ ---
gg <- unlist(str_split(dat2$Genes, ":"))
df.gg <- data.frame(Gene=names(table(gg)), Freq=as.numeric(table(gg)))
df.gg$Gene <- as.character(df.gg$Gene)
df.gg <- df.gg[order(df.gg$Freq, decreasing=TRUE),]
df.gg <- subset(df.gg, df.gg$Freq >= 3)

### DELETE UNWANTED GENES ---
del.index1 <- which(str_detect(df.gg$Gene, "-AS1") == TRUE)
del.index2 <- which(str_detect(df.gg$Gene, "MIR") == TRUE)
del.index3 <- which(str_detect(df.gg$Gene, "LOC") == TRUE)
del.index4 <- which(str_detect(df.gg$Gene, "LINC") == TRUE)
del.index5 <- which(str_detect(df.gg$Gene, "orf") == TRUE)

del.index <- unique(c(del.index1,del.index2,del.index3,del.index4,del.index5))
df.gg <- df.gg[-del.index,]

### REPLACE GENES ---
for(i in 1:nrow(dat2)){
	dat2$Genes[i] <- paste(intersect(str_split(dat2$Genes[i], ":")[[1]], df.gg$Gene), collapse=":")
}

### CREATE BI-PARTITE GRAPH ---
ids <- unique(dat2$SampleID)
bip <- matrix(0, nrow=length(ids), ncol=length(df.gg$Gene), dimnames=list(ids, df.gg$Gene))

for(i in 1:nrow(dat2)){
	id <- dat2$SampleID[i]
	status <- dat2$Status[i]
	genes <- str_split(dat2$Genes[i], ":")[[1]]
	index.g <- which(colnames(bip) %in% genes)
	index.s <- which(rownames(bip) == id)
	bip[index.s, index.g] <- status
}

bip <- t(bip)

### GET MART ANNOTATION FILE ---
file.annot <- file.path(dir.analysis, "genes_mart_export.tsv")
annot <- read.delim(file.annot, header=T, stringsAsFactors=F)
annot <- annot[order(annot$Start, decreasing=FALSE),]

genes <- intersect(rownames(bip), annot$Gene)

bip <- subset(bip, rownames(bip) %in% genes)
annot <- subset(annot, annot$Gene %in% genes)

bip <- bip[match(annot$Gene, rownames(bip)),]

df <- cbind(annot, bip)

### OUTPUT DATA ---
file.output <- file.path(dir.analysis, "MESO-PeM_chr3p21_status.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)
