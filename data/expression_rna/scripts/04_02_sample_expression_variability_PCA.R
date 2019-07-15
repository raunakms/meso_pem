### SetDirectories ----------------------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/expression")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")

### Define Files ------------------------------------------------------------
file.expr <- file.path(dir.analysis, "meso_peritoneal_gene_expression_proteincoding.tsv")

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

### Principal Component Analysis ----------------------------------------------
dat.pca <- prcomp(t(expr))

dat.eigenvectors <- as.data.frame(dat.pca$rotation[,1:3])
dat.eigenvectors$max.pc <- apply(dat.eigenvectors, 1, max)

cutoff <- sd(dat.eigenvectors$max.pc) * 2
dat.eigenvectors$var.select <- ifelse(dat.eigenvectors$max.pc > cutoff, 1, 0)
genes.select <- rownames(dat.eigenvectors)[which(dat.eigenvectors$var.select == 1)]

#### Select expression matrix with subset gene features ------------------------
expr.subset <- subset(expr, rownames(expr) %in% genes.select)

dat.pca.new <- prcomp(expr.subset)
dat.eigenvectors <- as.data.frame(dat.pca.new $rotation[,1:2])


###Plot
s <- unlist(lapply(str_split(rownames(dat.eigenvectors), "-"), function(x) x[2]))
file.plot <- file.path(dir.plot, "PCA_variable_expression.pdf")
pdf(file.plot, height=4, width=4)
plot(x=dat.eigenvectors$PC1, y=dat.eigenvectors$PC2,
		xlim=c(0.20,0.29), 
		ylim=c(-0.4,1.0),
		type="p", pch=1, col="black",
		cex=1, cex.axis=0.5, cex.lab=0.5, cex.main=0.5,
		las=1, tck = -.03, xlab="PC1", ylab="PC2",
		main="PCA: MESO Samples Using Variable Expression Genes")
text(x=dat.eigenvectors$PC1, y=dat.eigenvectors$PC2, labels=s, cex=0.5, pos=1, offset=0.5)
dev.off()


