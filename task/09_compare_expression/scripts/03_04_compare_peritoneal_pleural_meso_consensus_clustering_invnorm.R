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
file.dat <- file.path(dir.data, "expr_combined_meso_pleural_peritoneal_invnorm.tsv.gz")
file.cls <- file.path(dir.data, "expr_combined_meso_pleural_peritoneal.txt")
file.diffexpr <- file.path("/Data/Raunak/softwares/bdvtools/array_process/wilcox_rank_test.R")
file.des <- file.path(dir.data, "design_table.tsv")

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
sample.class <- read.delim(file.cls, header=F, stringsAsFactors=F)$V1
cls <- data.frame(SampleID=colnames(dat), Class=sample.class)
cls$SampleID <- as.character(cls$SampleID)
cls$Class <- as.character(cls$Class)

sampleids.pleu <- cls$SampleID[which(cls$Class == "MESO-PLEURAL")]
sampleids.peri <- cls$SampleID[which(cls$Class == "MESO-PERITONEAL")]

### DESIGN TABLE ---
des <- read.delim(file.des, header=TRUE, stringsAsFactors=FALSE)
des <- subset(des, des$SampleType == "Tumor")

### Principal Component Analysis --------------------------------------
expr <- dat
dat.pca <- prcomp(t(expr))

dat.eigenvectors <- as.data.frame(dat.pca$rotation[,1:3])
dat.eigenvectors$max.pc <- apply(dat.eigenvectors, 1, max)

cutoff <- sd(dat.eigenvectors$max.pc) * 1
dat.eigenvectors$var.select <- ifelse(dat.eigenvectors$max.pc > cutoff, 1, 0)
genes.select <- rownames(dat.eigenvectors)[which(dat.eigenvectors$var.select == 1)]

#> length(genes.select)
#[1] 10377
#> dim(expr)
#[1] 17667   101

#### Select expression matrix with subset gene features ----------------
expr.var <- subset(expr, rownames(expr) %in% genes.select)

### Consensus Cluster ---------------------------------------------------
dir.batch <- file.path(dir.clust, "02_MESO_expr_peritoneal_pleural_invnorm")
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
								Class=paste("Class_", as.numeric(dm[[i]][["consensusClass"]]), sep=""),
								Color=dm[[i]][["clrs"]][[1]])
	ctr <- ctr + 1
}
dclass <- do.call(rbind.data.frame, list.dclass)

#### Since, optimum k = 4 -----------------------------------------------
dclass <- subset(dclass, dclass$k == 6)
dclass$SampleID <- as.character(dclass$SampleID)
dclass$Class <- as.character(dclass$Class)
dclass$Subtype <- ""

for(i in 1:86){
	dclass$Subtype[i] <- des$Subtype[which(des$SampleID %in% dclass$SampleID[i])]
}


#> table(dclass$Subtype, dclass$Class)
#
#                    Class_1 Class_2 Class_3 Class_4 Class_5 Class_6
#                          0       0       0       1      15       0
#  Biphasic                7       1       0      13       0       0
#  DiffuseMaglignant       2       1       0       2       0       0
#  Epithelioid            32       5       8       9       2       1
#  Sarcomatoid             0       0       0       2       0       0


write.table(dclass, file.path(dir.clust, "02_MESO_expr_peritoneal_pleural_invnorm_cluster.tsv"), sep="\t", row.names=F, col.names=T, quote=F)


### Measure Correlation and Distance -------------------------
dat.pcc <- cor(expr.var, method="pearson")
dat.dist <- pr_simil2dist(dat.pcc)

### Cluster data ----------------------------------------------
hc <- hclust(dist(dat.dist, method = "euclidean"), method = "average")


### PLOT -----------------------------------------------
file.plot <- file.path(dir.plot, "meso_dendogram_expression_peri_pleu_invnorm.pdf")
pdf(file.plot, height=4, width=10)
	plot(hc, cex = 0.5, hang=0.1, xlab="", ylab="", main="", axes=FALSE)
dev.off()

