#### DEFINE PATH ------------------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/05_consensus_clustering")
dir.data <- file.path(dir.wrk, "data")
dir.plot <- file.path(dir.wrk, "plot")
dir.output <- file.path(dir.wrk, "output")
dir.expr <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/expression/main_calls")

### Define Files ------------------------------------------------------------
file.dat <- file.path(dir.expr, "Gene_expression_HTSeq_DESeq-normalized.tsv.gz")
file.array <- file.path("/Data/Raunak/softwares/bdvtools/array_process/array_preprocess.R")

### Load Libraries ----------------------------------------------------------
library("stringr")
library("proxy")
library("ConsensusClusterPlus")

### Load Expression Data -----------------------------------------------------
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
colnames(dat)[8:ncol(dat)] <- str_replace(colnames(dat)[8:ncol(dat)], "[.]", "-")

### SELECT FEATURES ----------------------------------------------------------
features <- c("miRNA","misc_RNA","Mt_rRNA","Mt_tRNA","rRNA","snoRNA","snRNA")

dat <- subset(dat, dat$Type %in% features)
dat <- dat[,c(2,8:ncol(dat))]
colnames(dat)[1] <- "Gene"

### Process GeneExpression ----------------------------------------------------
source(file.array)
expr <- getunique.gene.expression(dat)
dqnorm <- getQuantile.normalize(expr)
gexpr <- log2(dqnorm)
gexpr[gexpr == "-Inf"] <- 0
colnames(gexpr) <- str_replace(colnames(gexpr), "[.]", "-")
expr <- gexpr


### Remove Genes with Expr = 0 or NA  ---------------------------------
y1 <- apply(expr, 1, function(x) length(which(x == 0)))
y2 <- apply(expr, 1, function(x) length(which(is.na(x))))
del.index1 <- which(y1 >= 10)
del.index2 <- which(y2 >= 1)
del.index <- unique(c(del.index1, del.index2))
expr <- expr[-del.index,]


### Principal Component Analysis --------------------------------------
dat.pca <- prcomp(t(expr))

dat.eigenvectors <- as.data.frame(dat.pca$rotation[,1:3])
dat.eigenvectors$max.pc <- apply(dat.eigenvectors, 1, max)

cutoff <- sd(dat.eigenvectors$max.pc) * 1.5
dat.eigenvectors$var.select <- ifelse(dat.eigenvectors$max.pc > cutoff, 1, 0)
genes.select <- rownames(dat.eigenvectors)[which(dat.eigenvectors$var.select == 1)]


#### Select expression matrix with subset gene features ----------------
expr.var <- subset(expr, rownames(expr) %in% genes.select)

### Consensus Cluster ---------------------------------------------------
dir.batch <- file.path(dir.plot, "04_MESO_consensus_cluster_expression_snrna")
dm <- ConsensusClusterPlus(d=as.matrix(expr.var), maxK=10, reps=10000, pItem=0.8, pFeature=1,
							clusterAlg="hc", distance="pearson", seed=12345,
							innerLinkage="average", finalLinkage="average",
							title=dir.batch, plot="pdf", verbose=TRUE)


														
### Get Cluster Sample Class ---------------------------------------------
list.dclass <- list()
ctr <- 1
for(i in 2:10){
	list.dclass[[ctr]] <- data.frame(k=i,
								SampleID=names(dm[[i]][["consensusClass"]]),
								Class=as.numeric(dm[[i]][["consensusClass"]]),
								Color=dm[[i]][["clrs"]][[1]])
	ctr <- ctr + 1
}
dclass <- do.call(rbind.data.frame, list.dclass)


#### Since, optimum k = 2 -----------------------------------------------
dclass <- subset(dclass, dclass$k == 2)
dclass <- dclass[order(dclass$Class, decreasing=F),]
write.table(dclass, file.path(dir.output, "meso_consensus_cluster_expression_snrna_cluster.tsv"), sep="\t", row.names=F, col.names=T, quote=F)


### Measure Correlation and Distance -------------------------
dat.pcc <- cor(expr.var, method="pearson")
dat.dist <- pr_simil2dist(dat.pcc)

### Cluster data ----------------------------------------------
hc <- hclust(dist(dat.dist, method = "euclidean"), method = "average")


### PLOT -----------------------------------------------
file.plot <- file.path(dir.plot, "meso_dendogram_expression_snrna.pdf")
pdf(file.plot, height=4, width=5)
	plot(hc, cex = 0.7, hang=0.1, xlab="", ylab="", main="", axes=FALSE)
dev.off()


