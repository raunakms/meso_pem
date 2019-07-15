### LOAD LIBRAIES ---
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/08_compare_cnv")
dir.data <- file.path(dir.wrk, "data")
dir.output <- file.path(dir.wrk, "output")
dir.result <- file.path(dir.wrk, "results")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")
dir.batch <- file.path(dir.result, "CNV_MESO_OV_WCOX")

### DEFINE FILES ---
file.nexus <- file.path(dir.data, "nexus_comparision_meso_ov_parsed.tsv")
file.features <- file.path(dir.output, "cnv_meso_ov_diff_feature_genelist_wcox_nexus.txt")
file.annot <- file.path(dir.output, "cnv_meso_ov_diff_summary_table_by_chromosome.tsv")
files.dat <- list.files(dir.batch, pattern="performance_summary", full.names=FALSE)
genes <- unlist(lapply(str_split(str_replace_all(files.dat, "[.]", "_"), "_"), function(x) x[3]))

### CREATE MATRIX TO FILL ---
items.col <- c("AC","AUC","MCC","SN","SP","PPV","NPV","FPR","FNR","LRP","LRN")
mat <- matrix(NA, nrow=length(genes), ncol=11, dimnames=list(genes, items.col))

for(i in 1:length(files.dat)){

	# LOAD FILE ---
	file.dat <- file.path(dir.batch, files.dat[i])
	dat <- read.delim(file.dat, header=T, stringsAsFactors=F)

	# FILL MATRIX --- 
	mat[i,] <- dat$Mean

	cat("PROCESSED:", i, "OF", length(files.dat), "\n", sep="\t")
}

### WRITE OUTPUT ---
file.output <- file.path(dir.output, "cnv_meso_ov_diff_feature_performance_summary_allgenes.tsv")
write.table(mat, file.output, sep="\t", row.names=T, col.names=NA, quote=F) 


### LOAD FEATURES ---
#features <- read.table(file.features, header=F, stringsAsFactors=F)$V1

### SUBSET DATA ---
#mat.select <- subset(mat, rownames(mat) %in% features)

### WRITE OUTPUT ---
#file.output <- file.path(dir.output, "cnv_meso_ov_diff_feature_performance_summary.tsv")
#write.table(mat.select , file.output, sep="\t", row.names=T, col.names=NA, quote=F) 

### RE-LOAD MATRIX ----
#file.mat <- file.path(dir.output, "cnv_meso_ov_diff_feature_performance_summary.tsv")
#mat <- read.delim(file.mat, header=T, stringsAsFactors=F, row.names=1)

### ORDER MATRIX BY MCC --
#mat <- mat[order(mat$MCC, decreasing=TRUE),]

### FILTER MATRIX BY MCC >= 0.9 ---
#mat <- mat[which(mat$MCC >= 0.9),]

### LOAD ANNOTATION DATA ---
#annot <- read.delim(file.annot, header=T, stringsAsFactors=F)
#annot <- subset(annot, annot$Gene %in% rownames(mat))
#annot <- annot[match(rownames(mat), annot$Gene),]

### MERGE DATA ---
#df <- cbind(annot, mat)

### LOAD NESUS DATA ---
#dat.nexus <- read.delim(file.nexus, header=T, stringsAsFactors=F)

### FILL df DATA FRAME ---
#df$CNEvent <- ""
#df$Percent.Freq.MESO <- 0
#df$Percent.Freq.OV <- 0

#for(i in 1:nrow(df)){
#	gene <- df$Gene[i]
#	index <- which(str_detect(dat.nexus$Gene.Symbols, gene) == TRUE)
#	dtemp <- dat.nexus[index,]
#
#	df$CNEvent[i] <- paste(unique(dtemp$Event), collapse=":")
#	df$Percent.Freq.MESO[i] <- mean(dtemp$Percentage.Freq.in.Peritoneal)
#	df$Percent.Freq.OV[i] <- mean(dtemp$Percentage.Freq.in.Ovarian)
#}

### WRITE OUTPUT ---
#file.output <- file.path(dir.output, "cnv_meso_ov_diff_feature_performance_summary_compilation.tsv")
#write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F) 




### RELOAD DATA ---
#file.dat <- file.path(dir.output, "cnv_meso_ov_diff_feature_performance_summary_compilation.tsv")
#dat <- read.delim(file.dat, header=T, stringsAsFactors=F)

### LOAD FEATURE DATA ---
file.features <- file.path(dir.output, "cnv_meso_ov_diff_summary_table_by_chromosome_wcox_nexus.tsv")
dat.features <- read.delim(file.features, header=T, stringsAsFactors=F)

### TRIM DATA ---
mat.select <- subset(mat, rownames(mat) %in% dat.features$Gene)
mat.select <- mat.select[match(dat.features$Gene, rownames(mat.select)),]

### WRITE OUTPUT ---
file.output <- file.path(dir.output, "cnv_meso_ov_diff_feature_performance_summary.tsv")
write.table(mat.select, file.output, sep="\t", row.names=T, col.names=NA, quote=F) 


mat.pass1 <- subset(mat.select, mat.select[,2] > 0.5)
mat.pass2 <- subset(mat.pass1, mat.pass1[,3] > 0.70)


### WRITE OUTPUT ---
file.output <- file.path(dir.output, "cnv_meso_ov_diff_feature_performance_summary_finalpass30.tsv")
write.table(mat.pass2, file.output, sep="\t", row.names=T, col.names=NA, quote=F) 



### SUBSET ANNOTATION ----
file.dat <- file.path(dir.output, "cnv_meso_ov_diff_summary_table_by_chromosome_wcox_nexus_excel_export.tsv")
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
df <- subset(dat, dat$Gene %in% rownames(mat.pass2))
file.output <- file.path(dir.output, "cnv_meso_ov_diff_summary_table_by_chromosome_wcox_nexus_excel_export_finalpass30.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F) 
