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
file.sig <- file.path(dir.data, "DNA_REPAIR_genelist.txt")
file.expr <- file.path(dir.expr, "meso_peritoneal_gene_expression_invnorm_proteincoding.tsv.gz")
file.diffexpr <- file.path("/Data/Raunak/softwares/bdvtools/array_process/wilcox_rank_test.R")

### LOAD SIGNATURE GENES ---
genes.sig<- read.delim(file.sig, header=FALSE, stringsAsFactors=FALSE)$V1

### LOAD DATA ---
expr <- read.delim(file.expr, header=TRUE, stringsAsFactors=FALSE, row.names=1)
colnames(expr) <- str_replace_all(colnames(expr), "[.]", "-")
expr[is.na(expr)] <- 0

### LOAD DESIGN TABLE ---
des <- read.delim(file.des, header=TRUE, stringsAsFactors=FALSE)
grp0 <- intersect(colnames(expr), des$Sample.ID[which(des$Group.3p21 == 0)])
grp1 <- intersect(colnames(expr), des$Sample.ID[which(des$Group.3p21 == 1)])


### FUNCTION: get.heatmap ---
get.heatmap <- function(expr.sig, file.plot, plot.height, plot.width){
	require("gplots")
	require("RColorBrewer")

	# MANAGE GROUPS ---
	jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))
	label.groups <- c(rep("#fffc00", length(grp0)), rep("blue", length(grp1)))

	# PLOT ---
	pdf(file.plot, height=plot.height, width=plot.width)
	heatmap.2(as.matrix(t(expr.sig)), 
          col = rev(jColFun(1024)),
          Colv=TRUE, Rowv=TRUE, 
          dendrogram="both", trace="none",  scale="none",
          cexCol=0.8, cexRow=0.8, symbreaks=TRUE, margin=c(20,5),
          hclustfun = function(x) hclust(x, method = "ward.D2"), 
          distfun = function(x) dist(x, method = "euclidean"),
          colsep=c(1:500), rowsep=c(1:500),
          sepcolor="black", sepwidth=c(0.0005,0.0005), 
		  RowSideColors=label.groups,  xlab="", ylab="",
          key="TRUE", keysize=1, density.info="none", symkey=0,
		  key.title = NA, key.xlab = NA, key.ylab = NA,
		  key.par=list(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0)))
	#legend("topright", legend=c("BAP1-Wt","BAP1-Del"), 
	#	  fill=c("#fffc00","blue"), border=TRUE, bty="n", 
	#	  x.intersp = 1, y.intersp = 1, cex=1)
	dev.off()
}


### ARRANGE SAMPLES ---
expr <- subset(expr, select=c(grp0, grp1))
expr.sig <- subset(expr, rownames(expr) %in% genes.sig)
file.plot <- file.path(dir.plot, "heatmap_DNARepair_mRNA.pdf")
get.heatmap(expr.sig, file.plot, plot.height=5, plot.width=7)

