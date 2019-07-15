### DEFINE LIBRARIES ---
library("stringr")
library("gplots")
library("gridExtra")
library("reshape2")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.task <- file.path(dir.wrk, "task/10_integrated_pathways")
dir.data <- file.path(dir.task, "data")
dir.analysis <- file.path(dir.task, "output")
dir.plot <- file.path(dir.task, "plots")
dir.expr <- file.path(dir.wrk, "data/expression/analysis")
dir.des <- file.path(dir.wrk, "data/annotation")

### DEFINE FILES ---
file.des <- file.path(dir.des, "design_table_3p21genes.tsv")
file.expr <- file.path(dir.expr, "meso_peritoneal_gene_expression_invnorm_proteincoding.tsv.gz")
file.diffexpr <- file.path("/Data/Raunak/softwares/bdvtools/array_process/wilcox_rank_test.R")

### FUNCTION: get.heatmap ---
get.heatmap <- function(expr.sig, file.plot, plot.height, plot.width, label.groups){
	require("gplots")
	require("RColorBrewer")

	# MANAGE GROUPS ---
	jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))
	#label.groups <- c(rep("#fffc00", length(grp0)), rep("blue", length(grp1)))

	# PLOT ---
	pdf(file.plot, height=plot.height, width=plot.width)
	heatmap.2(as.matrix(t(expr.sig)), 
          col = rev(jColFun(1024)),
          Colv=TRUE, Rowv=TRUE, 
          dendrogram="both", trace="none",  scale="none",
          cexCol=0.8, cexRow=0.8, symbreaks=TRUE, margin=c(5,5),
          hclustfun = function(x) hclust(x, method = "ward.D2"), 
          distfun = function(x) dist(x, method = "euclidean"),
          colsep=c(1:500), rowsep=c(1:500),
          sepcolor="black", sepwidth=c(0.0005,0.0005), 
		  ColSideColors=label.groups,  
		  xlab="", ylab="",
          key="TRUE", keysize=1, density.info="none", symkey=0,
		  key.title = NA, key.xlab = NA, key.ylab = NA,
		  key.par=list(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0)))
	#legend("topright", legend=c("BAP1-Wt","BAP1-Del"), 
	#	  fill=c("#fffc00","blue"), border=TRUE, bty="n", 
	#	  x.intersp = 1, y.intersp = 1, cex=1)
	dev.off()
}


### GENES ---
#genes.sig <- c("LAG3","CD274","PDCD1LG2","TNFRSF4","TNFRSF9","PDCD1","CTLA4","ICOS","CD247","TLR9","SAGE1")
genes.sig <- c("LAG3","CD274","PDCD1LG2","PDCD1","CTLA4","ICOS","CD247")

### LOAD DATA ---
expr <- read.delim(file.expr, header=TRUE, stringsAsFactors=FALSE, row.names=1)
colnames(expr) <- str_replace_all(colnames(expr), "[.]", "-")
expr[is.na(expr)] <- 0


### LOAD DESIGN TABLE ---
des <- read.delim(file.des, header=TRUE, stringsAsFactors=FALSE)
grp0 <- intersect(colnames(expr), des$Sample.ID[which(des$Group.3p21 == 0)])
grp1 <- intersect(colnames(expr), des$Sample.ID[which(des$Group.3p21 == 1)])


### ARRANGE SAMPLES ---
#expr <- subset(expr, select=c(grp0, grp1))

expr0 <- subset(expr, select=grp0)
expr1 <- subset(expr, select=grp1)
expr.sig0 <- subset(expr0, rownames(expr0) %in% genes.sig)
expr.sig1 <- subset(expr1, rownames(expr1) %in% genes.sig)

expr.sig0 <- expr.sig0[match(genes.sig, row.names(expr.sig0)),]
expr.sig1 <- expr.sig1[match(genes.sig, row.names(expr.sig1)),]


file.plot0 <- file.path(dir.plot, "heatmap_ImmuneCheckpointM_mRNA_grp0.pdf")
get.heatmap(expr.sig0, file.plot0, plot.height=4, plot.width=3, label.groups=rep("#fffc00", length(grp0)))

file.plot1 <- file.path(dir.plot, "heatmap_ImmuneCheckpointM_mRNA_grp1.pdf")
get.heatmap(expr.sig1, file.plot1, plot.height=4, plot.width=3, label.groups=rep("blue", length(grp1)))


expr01 <- subset(expr, select=c(grp0,grp1))
expr.sig01 <- subset(expr01, rownames(expr01) %in% genes.sig)
expr.sig01 <- expr.sig01[match(genes.sig, row.names(expr.sig01)),]
file.plot01 <- file.path(dir.plot, "heatmap_ImmuneCheckpointM_mRNA.pdf")
get.heatmap(t(expr.sig01), file.plot01, plot.height=4, plot.width=3, label.groups=c(rep("#fffc00", length(grp0)), rep("blue", length(grp1))) )


