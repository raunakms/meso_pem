### Load Libraries ----------------------------------------------------------
library("stringr")

### SetDirectories ----------------------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/expression")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")

### Define Files ------------------------------------------------------------
file.dat <- file.path(dir.data, "Gene_expression_HTSeq_DESeq-normalized.tsv.gz")
file.array <- file.path("/Data/Raunak/softwares/bdvtools/array_process/array_preprocess.R")

### Load Expression Data -----------------------------------------------------
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
colnames(dat)[8:ncol(dat)] <- str_replace(colnames(dat)[8:ncol(dat)], "[.]", "-")

dat <- subset(dat, dat$Type == "protein_coding")
dat <- dat[,c(2,8:ncol(dat))]
colnames(dat)[1] <- "Gene"

### Process GeneExpression ----------------------------------------------------
source(file.array)
expr <- getunique.gene.expression(dat)

dqnorm <- getQuantile.normalize(expr)
colnames(dqnorm) <- colnames(expr)
file.output <- file.path(dir.analysis, "meso_peritoneal_gene_expression_proteincoding_qnorm.tsv")
write.table(dqnorm, file.output, sep="\t", row.names=T, col.names=NA, quote=F)

gexpr <- log2(dqnorm)
gexpr[gexpr == "-Inf"] <- 0

file.output <- file.path(dir.analysis, "meso_peritoneal_gene_expression_proteincoding.tsv")
write.table(gexpr, file.output, sep="\t", row.names=T, col.names=NA, quote=F)

exprz <- getZscore(gexpr)
exprz[exprz == "NaN"] <- 0

file.output <- file.path(dir.analysis, "meso_peritoneal_gene_expression_znorm_proteincoding.tsv")
write.table(exprz, file.output, sep="\t", row.names=T, col.names=NA, quote=F)


###
file.output <- file.path(dir.analysis, "meso_peritoneal_gene_expression_extended.tsv")
write.table(gexpr, file.output, sep="\t", row.names=T, col.names=NA, quote=F)

file.output <- file.path(dir.analysis, "meso_peritoneal_gene_expression_znorm_extended.tsv")
write.table(exprz, file.output, sep="\t", row.names=T, col.names=NA, quote=F)

