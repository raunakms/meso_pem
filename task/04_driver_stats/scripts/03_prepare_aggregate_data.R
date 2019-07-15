### DEFINE LIBRARIES ---
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.snv <- file.path(dir.wrk, "data/mutation/analysis")
dir.cnv <- file.path(dir.wrk, "data/cnv/analysis")
dir.fus <- file.path(dir.wrk, "data/genefusion/analysis")
dir.task <- file.path(dir.wrk, "task/04_driver_stats")
dir.output <- file.path(dir.wrk, "output")

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

### WRTIE OUTPUT ---
file.output <- file.path(dir.output, "meso_alterations_samplewise_status_all.tsv")
write.table(lpar, file.output, sep="\t", row.names=F, col.names=T, quote=F)
