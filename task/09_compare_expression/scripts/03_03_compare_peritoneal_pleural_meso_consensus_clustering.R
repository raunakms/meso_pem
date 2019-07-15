### LOAD LIBRARIES ----
library("stringr")
library("proxy")
library("ConsensusClusterPlus")

#### DEFINE PATH ----
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/09_compare_expression")
dir.data <- file.path(dir.wrk, "data")
dir.plot <- file.path(dir.wrk, "plot")
dir.clust <- file.path(dir.wrk, "consensus_cluster")

### DEFINE FILE ---
file.dat <- file.path(dir.data, "expr_combined_meso_pleural_peritoneal.tsv.gz")
file.des <- file.path(dir.data, "expr_combined_meso_pleural_peritoneal.txt")
file.diffexpr <- file.path("/Data/Raunak/softwares/bdvtools/array_process/wilcox_rank_test.R")

## Remove Genes with Expr = 0 in > 50% of the samples ------
remove.0.matrix <- function(dat, cutoff){
	y <- apply(dat, 1, function(x) length(which(x == 0)))
	del.index <- which(y > (ncol(dat) * cutoff))
	if(length(del.index) != 0){
		dat <- dat[-del.index,]
	} else{
		return(dat)
	}
	return(dat)
}

### LOAD DATA ---
dat <- read.delim(file.dat, header=T, stringsAsFactors=F, row.names=1)
colnames(dat) <- str_replace_all(colnames(dat), "[.]", "-")
dat <- remove.0.matrix(dat=dat, cutoff=0.25)

### SAMPLE CLASS ---
sample.class <- read.delim(file.des, header=F, stringsAsFactors=F)$V1
des <- data.frame(SampleID=colnames(dat), Class=sample.class)
des$SampleID <- as.character(des$SampleID)
des$Class <- as.character(des$Class)

sampleids.pleu <- des$SampleID[which(des$Class == "MESO-PLEURAL")]
sampleids.peri <- des$SampleID[which(des$Class == "MESO-PERITONEAL")]


### Principal Component Analysis --------------------------------------
expr <- dat
dat.pca <- prcomp(t(expr))

dat.eigenvectors <- as.data.frame(dat.pca$rotation[,1:3])
dat.eigenvectors$max.pc <- apply(dat.eigenvectors, 1, max)

cutoff <- sd(dat.eigenvectors$max.pc) * 1
dat.eigenvectors$var.select <- ifelse(dat.eigenvectors$max.pc > cutoff, 1, 0)
genes.select <- rownames(dat.eigenvectors)[which(dat.eigenvectors$var.select == 1)]

#> length(genes.select)
#[1] 5552
#> dim(expr)
#[1] 17667   101

#### Select expression matrix with subset gene features ----------------
expr.var <- subset(expr, rownames(expr) %in% genes.select)

### Consensus Cluster ---------------------------------------------------
dir.batch <- file.path(dir.clust, "01_MESO_expr_peritoneal_pleural")
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

#### Since, optimum k = 3 -----------------------------------------------
dclass <- subset(dclass, dclass$k == 3)


### Measure Correlation and Distance -------------------------
dat.pcc <- cor(expr.var, method="pearson")
dat.dist <- pr_simil2dist(dat.pcc)

### Cluster data ----------------------------------------------
hc <- hclust(dist(dat.dist, method = "euclidean"), method = "average")


### PLOT -----------------------------------------------
file.plot <- file.path(dir.plot, "meso_dendogram_expression_peri_pleu.pdf")
pdf(file.plot, height=4, width=10)
	plot(hc, cex = 0.5, hang=0.1, xlab="", ylab="", main="", axes=FALSE)
dev.off()

