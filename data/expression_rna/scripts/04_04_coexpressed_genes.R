### SetDirectories ----------------------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/expression")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")

### Define Files ------------------------------------------------------------
file.expr <- file.path(dir.analysis, "meso_peritoneal_gene_expression_znorm_proteincoding.tsv")

### Load Libraries ----------------------------------------------------------
library("stringr")
library("ggplot2")

### Load Expression Data -----------------------------------------------------
expr <- read.delim(file.expr, header=T, row.names=1,  stringsAsFactors=F)
colnames(expr) <- str_replace(colnames(expr), "[.]", "-")

### Remove Genes with Expr = 0 in >4 samples  ---------------------------------
y <- apply(expr, 1, function(x) length(which(x == 0)))
del.index <- which(y > 4)
expr <- expr[-del.index,]

### Compute Pearson Correlation between genes  ---------------------------------
dat.pcc <- cor(t(expr), method = "pearson")

dat.pcc[(dat.pcc < 0.8) & (dat.pcc > -0.8)] <- 0

for(i in 1:nrow(dat.pcc)){
	df <- data.frame(Gene=colnames(dat.pcc)[-i], pcc=as.numeric(dat.pcc[i,-i]))
	df$Gene <- as.character(df$Gene)
	df <- subset(df, df$pcc != 0)
}
