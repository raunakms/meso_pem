### DEFINE LIBRARIES --
library("stringr")

### DEFINe PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/mutation")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")

### DEFINE FILES ---
file.maf <- file.path(dir.data, "VPCLAGA.MESO_PeM.nonsilent_mutation.maf")

# LOAD DATA ---
dat.maf <- read.delim(file.maf, header=TRUE, stringsAsFactors=FALSE)

### TRIM DATA ---
dat <- subset(dat.maf, select=c("Tumor_Sample_Barcode","Hugo_Symbol","Variant_Type"))
dat$Variant_Type[which(dat$Variant_Type %in% c("INS","DEL"))] <- "INDEL"
dat$Variant_Type <- paste(dat$Variant_Type, ";", sep="")

### DISCARD REDUNDANT DATA ---
df <- dat[,1:2]
df <- df[!duplicated(df),]

df$Status.INDEL <- df$Status.SNV <- ""

### ADD SNV DATA ---
for(i in 1:nrow(dat)){
	id <- dat$Tumor_Sample_Barcode[i]
	gene <- dat$Hugo_Symbol[i]
	s <- dat$Variant_Type[i]

	if(s == "SNP;"){
		df$Status.SNV[which((df$Tumor_Sample_Barcode == id) & (df$Hugo_Symbol == gene))] <- s
	}
}

### ADD INDEL DATA ---
for(i in 1:nrow(dat)){
	id <- dat$Tumor_Sample_Barcode[i]
	gene <- dat$Hugo_Symbol[i]
	s <- dat$Variant_Type[i]

	if(s == "INDEL;"){
		df$Status.INDEL[which((df$Tumor_Sample_Barcode == id) & (df$Hugo_Symbol == gene))] <- s
	}
}

### MERGE ALTERATIONS ---
remove.null <- function(y){ return(subset(y, y != "")) }
df$Status <- apply(df, 1, function(x) paste(remove.null(c(x[3],x[4])), collapse=""))
colnames(df) <- c("SampleID","Gene","V1","V2", "Status")

### TRIM DATA ---
df <- subset(df, select=c("SampleID","Gene","Status"))

### WRITE OUTPUT ---
file.output <- file.path(dir.analysis, "snv_alterations_status.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)

### LOAD DRIVER GENES ---
dir.task <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/04_driver_stats/output")
file.driver <- file.path(dir.task, "meso_driver_genelist.txt")
genes.driver <- read.table(file.driver, header=F, stringsAsFactors=F)$V1

### RETAIN ONLY DRIVER GENES ---
df <- subset(df, df$Gene %in% genes.driver)

### OUTPUT DATA ---
file.output <- file.path(dir.analysis, "MESO-PeM_drivergene_snv_status.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)
