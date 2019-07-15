### DEFINE LIBRARIES ---
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.snv <- file.path(dir.wrk, "data/mutation/analysis")
dir.cnv <- file.path(dir.wrk, "data/cnv/analysis")
dir.fus <- file.path(dir.wrk, "data/genefusion/analysis")
dir.task <- file.path(dir.wrk, "task/10_integrated_pathways")
dir.data <- file.path(dir.task, "data")

### DEFINE FILES ---
file.snv <- file.path(dir.snv, "snv_alterations_status.tsv") 
file.cnv <- file.path(dir.cnv, "meso_cnv_genes.tsv") 
file.fus <- file.path(dir.fus, "fusion_alteration_status.tsv") 

### LOAD SNV ----
dat.snv <- read.delim(file.snv, header=TRUE, stringsAsFactors=FALSE)

### LOAD CNV ----
dat.cnv <- read.delim(file.cnv, header=TRUE, stringsAsFactors=FALSE)
dat.cnv <- subset(dat.cnv, select=c("SampleID","Gene","Status"))
dat.cnv <- dat.cnv[!duplicated(dat.cnv),]

dat.cnv$Type <- ""
dat.cnv$Type[which(dat.cnv$Status %in% c("CN Gain","High Copy Gain"))]  <- "AMP;"
dat.cnv$Type[which(dat.cnv$Status %in% c("CN Loss","Homozygous Copy Loss"))]  <- "DEL;"

dat.cnv <- dat.cnv[,-3]
colnames(dat.cnv)[3] <- "Status"

### LOAD SNV ----
dat.fus <- read.delim(file.fus, header=TRUE, stringsAsFactors=FALSE)
dat.fus$Status <- "FUSION;"

### AGGREGATE DATA ---
dat <- rbind(dat.snv, dat.cnv, dat.fus)
dat <- dat[order(dat$SampleID, decreasing=FALSE),]

### WRTIE OUTPUT ---
file.output <- file.path(dir.data, "meso_alterations_samplewise_status_all.tsv")
write.table(dat, file.output, sep="\t", row.names=F, col.names=T, quote=F)

### LOAD PATHWAY GENE FILE ---
file.sig <- file.path(dir.data, "pathway_gene_relationship.tsv")
dat.sig <- read.delim(file.sig, header=TRUE, stringsAsFactors=FALSE)
dat.sig$Group <- apply(dat.sig, 1, function(x) paste(x[1], x[2], sep=":"))

### SUBSET DATA ---
df <- subset(dat, dat$Gene %in% dat.sig$Gene)

### ADD GROUP LABEL ---
df$Group <- ""
for(i in 1:nrow(df)){
	df$Group[i] <- dat.sig$Group[which(dat.sig$Gene == df$Gene[i])]
}

### TRIM DATA ---
df <- subset(df, select=c("SampleID","Group","Status"))
colnames(df)[2] <- "Gene"

df$Status <- str_replace(df$Status, "SNP;INDEL;", "MUT;")
df$Status <- str_replace(df$Status, "SNP;", "MUT;")
df$Status <- str_replace(df$Status, "INDEL;", "IND;")


### WRTIE OUTPUT ---
file.output <- file.path(dir.data, "aggregate_alteration_oncoprint_by_pathway_input.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)
