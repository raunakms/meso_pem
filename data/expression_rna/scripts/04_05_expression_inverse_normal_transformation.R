### Load Libraries ----------------------------------------------------------
library("stringr")
library("GenABEL")

### DEFINE PATH ----
dir.wrk <- file.path("//JBRCSRV009/CollinsGroup/Raunak/projects/MESO_peritoneal/data/expression")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")

### DEFINE FILE ----
file.expr <- file.path(dir.analysis, "meso_peritoneal_gene_expression_proteincoding.tsv.gz")
file.array <- file.path("//JBRCSRV009/CollinsGroup/Raunak/softwares/bdvtools/array_process/array_preprocess.R")

### LOAD DATA ---
expr <- read.delim(file.expr, header=T, stringsAsFactors=F, row.names=1)
colnames(expr) <- str_replace_all(colnames(expr), "[.]", "-")

### INVERSE-NORMAL TRANSFORMATION ----
source(file.array)
expr.norm <- getInvNormTransform(dat=expr)

### WRITE OUTPUT ---
file.output <- file.path(dir.analysis, "meso_peritoneal_gene_expression_invnorm_proteincoding.tsv")
write.table(expr.norm, file.output, sep="\t", row.names=T, col.names=NA, quote=F)

### COMPRESS ---
cmd <- paste("gzip", file.output, sep=" ")
system(cmd)
