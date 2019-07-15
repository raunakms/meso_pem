### Load Libraries ----------------------------------------------------------
library("stringr")
library("GenABEL")

### DEFINE PATH ----
dir.wrk <- file.path("//jbrcsrv009/CollinsGroup/Raunak/projects/MESO_peritoneal/data/proteome")
dir.data <- file.path(dir.wrk, "data_proteome_discover")
dir.analysis <- file.path(dir.wrk, "processed_data")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")

### DEFINE FILE ----
file.expr <- file.path(dir.analysis, "meso_proteome_proteomediscover_all_log2.tsv.gz")
file.array <- file.path("//JBRCSRV009/CollinsGroup/Raunak/softwares/bdvtools/array_process/array_preprocess.R")

### LOAD DATA ---
expr <- read.delim(file.expr, header=T, stringsAsFactors=F, row.names=1)
colnames(expr) <- str_replace_all(colnames(expr), "[.]", "-")

### INVERSE-NORMAL TRANSFORMATION ----
source(file.array)
expr.norm <- getInvNormTransform(dat=expr)

### WRITE OUTPUT ---
file.output <- file.path(dir.analysis, "meso_proteome_proteomediscover_all_invnorm.tsv")
write.table(expr.norm, file.output, sep="\t", row.names=T, col.names=NA, quote=F)

### COMPRESS ---
cmd <- paste("gzip", file.output, sep=" ")
system(cmd)
