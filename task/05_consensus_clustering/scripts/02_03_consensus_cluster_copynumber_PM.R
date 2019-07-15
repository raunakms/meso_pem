#### LOAD LIBRARIES ---------------------------------------------------
library("stringr")
library("dendextend")
library("ConsensusClusterPlus")

#### DEFINE PATH ------------------------------------------------------
dir.wrk <- file.path("//jbrcsrv009/CollinsGroup/Raunak/projects/MESO_peritoneal/task/05_consensus_clustering")
#dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/05_consensus_clustering")
dir.data <- file.path(dir.wrk, "data")
dir.plot <- file.path(dir.wrk, "plot")
dir.output <- file.path(dir.wrk, "output")
#dir.cnv <- file.path("/Data/Raunak/projects/HITnDRIVE/datasets/TCGA_MESO/data/cnv/seq_call_refseq_genes_meso")
dir.cnv <- file.path("//jbrcsrv009/CollinsGroup/Raunak/projects/HITnDRIVE/datasets/TCGA_MESO/data/cnv/seq_call_refseq_genes_meso")

#### DEFINE FILE ------------------------------------------------------
file.dat <- file.path(dir.cnv, "tcga_meso_cnv_seg_values_calls_parsed.tsv.gz")

## Remove Genes with Expr = NA in > 25% of the samples ------
remove.na.matrix <- function(dat, cutoff){
	y <- apply(dat, 1, function(x) length(which(is.na(x))))
	del.index <- which(y > (ncol(dat) * cutoff))
	if(length(del.index) != 0){
		dat <- dat[-del.index,]
	} else{
		return(dat)
	}
	return(dat)
}

#### Load Expression Data ---------------------------------------------
dat <- read.delim(file.dat, header=T, row.names=1)
colnames(dat) <- str_replace(colnames(dat), "[.]", "-")
dat <- remove.na.matrix(dat=dat, cutoff=0.25)
dat[is.na(dat)] <- 0

### Principal Component Analysis --------------------------------------
dat.pca <- prcomp(t(dat))

dat.eigenvectors <- as.data.frame(dat.pca$rotation[,1:3])
dat.eigenvectors$max.pc <- apply(dat.eigenvectors, 1, max)

cutoff <- sd(dat.eigenvectors$max.pc) * 2
dat.eigenvectors$var.select <- ifelse(dat.eigenvectors$max.pc > cutoff, 1, 0)
genes.select <- rownames(dat.eigenvectors)[which(dat.eigenvectors$var.select == 1)]


#### Select expression matrix with subset gene features ----------------
dat.var <- subset(dat, rownames(dat) %in% genes.select)

### Consensus Cluster ---------------------------------------------------
dir.batch <- file.path(dir.plot, "02_03_MESO_consensus_cluster_cnv_PM")
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

write.table(dclass, file.path(dir.batch, "cluster_groups.tsv"), sep="\t", row.names=F, col.names=T, quote=F)


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
