### DEFINE PATH ---------------------------------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data")
dir.cna <- file.path(dir.wrk, "cnv/analysis")
dir.exp <- file.path(dir.wrk, "expression/analysis")

### DEFINE FILES --------------------------------------------------------------------
file.cna <- file.path(dir.cna, "meso_cnv_genes.tsv")
file.exp <- file.path(dir.exp, "meso_peritoneal_gene_expression_znorm_extended.tsv")

### Load Libraries ----------------------------------------------------------
library("stringr")

### LOAD FILES -----------------------------------------------------------------------
dat.cna <- read.delim(file.cna, header=T, stringsAsFactors=F)
dat.cna <- dat.cna[,c(1,10,5,9)]
dat.cna <- subset(dat.cna, dat.cna$SampleID != "MESO-20")

dat.exp <- read.delim(file.exp, header=T, stringsAsFactors=F, row.names=1)
colnames(dat.exp) <- str_replace(colnames(dat.exp), "[.]", "-")

### Separate Data by Overlapped SampleIDs ---------------------------------------------
ids.cna <- unique(dat.cna$SampleID)
ids.exp <- colnames(dat.exp)

ids.cna_exp <- intersect(ids.cna, ids.exp)
ids.cna_only <- setdiff(ids.cna, ids.exp)
ids.exp_only <- setdiff(ids.exp, ids.cna)

dat.cna_both <- subset(dat.cna, dat.cna$SampleID %in% ids.cna_exp)
dat.cna_only <- subset(dat.cna, dat.cna$SampleID %in% ids.cna_only)

### FILL-UP GENE-EXPRESSION PROFILE VALUES FOR CNA GENES --------------------------------
dat.cna_both$Expr <- NA

for(i in 1:nrow(dat.cna_both)){
	x <- dat.cna_both$Gene[i]
	y <- dat.cna_both$SampleID[i]
	
	exp.val <- as.numeric(dat.exp[x,y])
	
	if(!is.na(exp.val)){
		dat.cna_both$Expr[i] <- exp.val
	}
}

### Find Consistent Direction of Expression -----------------------------------------------
dat.cna_both <- subset(dat.cna_both, !is.na(dat.cna_both$Expr))

dat.cna_both$Cor.Status <- 0
for(i in 1:nrow(dat.cna_both)){
	if((dat.cna_both$Probe.Median[i] > 0) & (dat.cna_both$Expr[i] > 0)){
		dat.cna_both$Cor.Status[i] <- 1
	}else if((dat.cna_both$Probe.Median[i] < 0) & (dat.cna_both$Expr[i] < 0)){
		dat.cna_both$Cor.Status[i] <- 1
	}else{
		dat.cna_both$Cor.Status[i] <- 0
	}
	#cat("PROCESSED:", i, "OF", nrow(dat.cna_both), "\n", sep="\t")	
}

### Retain rows that correlate ---------------------------------------------------------------
dat.cna_both <- subset(dat.cna_both, dat.cna_both$Cor.Status == 1)

file.output <- file.path(dir.cna, "correlation_cna_expr.tsv")
write.table(dat.cna_both, file.output, sep="\t", row.names=F, col.names=T, quote=F)

dat.cna_both_amp <- subset(dat.cna_both, dat.cna_both$Status %in% c("CN Gain","High Copy Gain"))
dat.cna_both_del <- subset(dat.cna_both, dat.cna_both$Status %in% c("CN Loss","Homozygous Copy Loss"))


### Stat: CNA Genes ------------------------------------------------------------
dat.cna_both <- dat.cna_both[,-6]
dat.cna_only$Expr <- NA
dat.cna_both_only <- rbind(dat.cna_both, dat.cna_only)

dstat2 <- data.frame(Gene=names(table(dat.cna_both_only$Gene)), Freq=as.numeric(table(dat.cna_both_only$Gene)))
dstat2 <- dstat2[order(dstat2$Freq, decreasing=T),]

subset(dat.cna_both_only, dat.cna_both_only$Gene == "BAP1")


### Prepare Data for HIT'nDRIVE -------------------------------------------------
dat <- dat.cna_both_only
dat$CNA.Status <- dat$Status
dat$Status <- ifelse((dat$CNA.Status == "CN Gain") | (dat$CNA.Status == "High Copy Gain"), "AMP;", "DEL;")

df.alt <- dat[,c(1,2,3,6)]

file.output <- file.path(dir.cna, "alterations_samplewise_cna.tsv")
write.table(df.alt, file.output, sep="\t", row.names=F, col.names=T, quote=F)

### (1) Stat: CNV --------------------------------------------------------------
file.df <- file.path(dir.cna, "cnv_samplewise_stats.tsv")
df <- read.delim(file.df, header=T, stringsAsFactors=F)[-20,]
dstat1 <- data.frame(SampleID=names(table(df.alt$SampleID)), 
					FreqGeneFiltered=as.numeric(table(df.alt$SampleID)))
dstat1$SampleID <- as.character(dstat1$SampleID)

df$FreqGeneFiltered <- dstat1$FreqGeneFiltered

file.output <- file.path(dir.cna, "cnv_samplewise_stats.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)


