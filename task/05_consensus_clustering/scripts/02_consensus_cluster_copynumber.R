#### LOAD LIBRARIES ---------------------------------------------------
library("stringr")
library("proxy")
library("dendextend")
library("ConsensusClusterPlus")

#### DEFINE PATH ------------------------------------------------------
#dir.wrk <- file.path("//jbrcsrv009/CollinsGroup/Raunak/projects/MESO_peritoneal/task/05_consensus_clustering")
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/05_consensus_clustering")
dir.data <- file.path(dir.wrk, "data")
dir.plot <- file.path(dir.wrk, "plot")
dir.output <- file.path(dir.wrk, "output")

#### DEFINE FILE ------------------------------------------------------
file.dat <- file.path(dir.data, "meso_cnv_seg_values_calls_parsed.tsv")

#### Load Expression Data ---------------------------------------------
dat <- read.delim(file.dat, header=T, row.names=1)
colnames(dat) <- str_replace(colnames(dat), "[.]", "-")

### Remove Genes with CNA = NA  ---------------------------------
y <- apply(dat, 1, function(x) length(which(is.na(x))))
del.index <- which(y >= 1)
dat <- dat[-del.index,]

### Principal Component Analysis --------------------------------------
#expr <- na.omit(expr)
dat.pca <- prcomp(t(dat))

dat.eigenvectors <- as.data.frame(dat.pca$rotation[,1:3])
dat.eigenvectors$max.pc <- apply(dat.eigenvectors, 1, max)

cutoff <- sd(dat.eigenvectors$max.pc) * 2
dat.eigenvectors$var.select <- ifelse(dat.eigenvectors$max.pc > cutoff, 1, 0)
genes.select <- rownames(dat.eigenvectors)[which(dat.eigenvectors$var.select == 1)]


#### Select expression matrix with subset gene features ----------------
dat.var <- subset(dat, rownames(dat) %in% genes.select)

### Consensus Cluster ---------------------------------------------------
#dir.batch <- file.path(dir.plot, "02_MESO_consensus_cluster_cnv")
dir.batch <- file.path(dir.plot, "test")
dm <- ConsensusClusterPlus(d=as.matrix(dat.var), maxK=10, reps=10000, pItem=0.8, pFeature=1,
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
dclass$SampleID <- as.character(dclass$SampleID)
dclass$Class <- as.character(dclass$Class)
dclass$Color <- as.character(dclass$Color)

write.table(dclass, file.path(dir.plot, "02_MESO_consensus_cluster_cnv/02_MESO_consensus_cluster_cnv.tsv"), sep="\t", row.names=F, col.names=T, quote=F)


#### FUNCTION: ---
plot_dendo <- function(dat_class, k){
	# SUBSET BY k class ---
	df <- dat_class[which(dat_class$k == k),]

	# GET CLUSTER ---
	hc <- dm[[k]]$consensusTree
	ids <- names(dm[[k]]$consensusClass)
	color.ids <- dm[[k]]$clrs[[1]]
	#hc <- hclust(dist(dat.dist, method = "euclidean"), method = "average")
	dend <- as.dendrogram(hc)
	labels(dend) <- ids[hc$order]
	labels_colors(dend) <- color.ids[hc$order]

	# PLOT ----
	p <- plot(dend, cex = 0.7, xlab="", ylab="", main=paste("k =", k, sep=""), axes=FALSE)
	
	return(p)
}


### CALL FUNCTION ---
file.plot <- file.path(dir.batch, "meso_dendogram_cnv.pdf")
pdf(file.plot, height=3, width=5)
	plot_dendo(dat_class=dclass, k=2)
	plot_dendo(dat_class=dclass, k=3)
	plot_dendo(dat_class=dclass, k=4)
	plot_dendo(dat_class=dclass, k=5)
	plot_dendo(dat_class=dclass, k=6)
dev.off()
