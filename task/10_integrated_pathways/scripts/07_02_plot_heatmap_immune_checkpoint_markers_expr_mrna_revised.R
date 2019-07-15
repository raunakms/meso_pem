### DEFINE LIBRARIES ---
library("stringr")
library("gplots")
library("gridExtra")
library("reshape2")
library("ggplot2")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.task <- file.path(dir.wrk, "task/10_integrated_pathways")
dir.data <- file.path(dir.task, "data")
dir.analysis <- file.path(dir.task, "output")
dir.plot <- file.path(dir.task, "plots")
dir.expr <- file.path(dir.wrk, "data/expression/analysis")
dir.des <- file.path(dir.wrk, "data/annotation")
dir.diff <- file.path(dir.wrk, "task/17_expression_3p21/output")

### DEFINE FILES ---
file.des <- file.path(dir.des, "design_table_3p21genes.tsv")
file.expr <- file.path(dir.expr, "meso_peritoneal_gene_expression_invnorm_proteincoding.tsv.gz")
file.diffexpr <- file.path(dir.diff, "diffexpr_3p21gene_status_wilcoxranktest_mrna.tsv")

### FUNCTION: get.heatmap ---
get.heatmap <- function(expr.sig, file.plot, plot.height, plot.width, label.groups){
	require("gplots")
	require("RColorBrewer")

	# MANAGE GROUPS ---
	jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))
	#label.groups <- c(rep("#fffc00", length(grp0)), rep("blue", length(grp1)))

	# PLOT ---
	pdf(file.plot, height=plot.height, width=plot.width)
	heatmap.2(as.matrix(expr.sig), 
          col = rev(jColFun(1024)),
          Colv=FALSE, Rowv=TRUE, 
          dendrogram="row", trace="none",  scale="none",
          cexCol=0.8, cexRow=0.8, symbreaks=TRUE, margin=c(12,8),
          hclustfun = function(x) hclust(x, method = "ward.D2"), 
          distfun = function(x) dist(x, method = "euclidean"),
          colsep=c(1:500), rowsep=c(1:500),
          sepcolor="black", sepwidth=c(0.0005,0.0005), 
		  ColSideColors=label.groups,  
		  xlab="", ylab="",
          key="TRUE", keysize=1, density.info="none", symkey=0,
		  key.title = NA, key.xlab = NA, key.ylab = NA,
		  key.par=list(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0)))
	dev.off()
}

#TUMOR					T-CELL
#"CD274","PDCD1LG2"		"PDCD1"
#"CD80","CD86"			"CD28","CTLA4"
#"ICOSLG"				"ICOS"
#"TNFRSF14"				"BTLA"
#"HLA-A","HLA-B","HLA-C","LAG3"
#"LGALS9"				"HAVCR2"
#"TNFSF9"				"TNFRSF9"
#"TNFRSF4"				"TNFSF4"
#"CD70"					"CD27"
#"CD40"					"CD40LG"


### GENES ---
genes.sig <- c("CD274","PDCD1LG2","PDCD1","CD80","CD86","CD28","CTLA4","ICOSLG","ICOS","TNFRSF14","BTLA","HLA-A","HLA-B","HLA-C","LAG3","LGALS9","HAVCR2")

ids <- c("MESO-08","MESO-18E","MESO-18A","MESO-19","MESO-13","MESO-12","MESO-11",
		"MESO-14","MESO-17","MESO-09","MESO-02","MESO-10","MESO-07","MESO-06","MESO-05")

### LOAD DATA ---
expr <- read.delim(file.expr, header=TRUE, stringsAsFactors=FALSE, row.names=1)
colnames(expr) <- str_replace_all(colnames(expr), "[.]", "-")
expr[is.na(expr)] <- 0

### LOAD DESIGN TABLE ---
des <- read.delim(file.des, header=TRUE, stringsAsFactors=FALSE)
grp0 <- intersect(colnames(expr), des$Sample.ID[which(des$Group.3p21 == 0)])
grp1 <- intersect(colnames(expr), des$Sample.ID[which(des$Group.3p21 == 1)])


### GET HEATMAP ---
expr01 <- subset(expr, select=c(grp0,grp1))
expr.sig01 <- subset(expr01, rownames(expr01) %in% genes.sig)
expr.sig01 <- subset(expr.sig01, select=ids)
#expr.sig01 <- expr.sig01[match(genes.sig, row.names(expr.sig01)),]
file.plot01 <- file.path(dir.plot, "heatmap_ImmuneCheckpointM_mRNA_revised.pdf")
get.heatmap(expr.sig01, file.plot01, plot.height=4, plot.width=3.5, label.groups=c(rep("#fffc00", length(grp0)), rep("blue", length(grp1))) )


### LOAD DIFF EXPR DATA ---
dat <- read.delim(file.diffexpr, header=TRUE, stringsAsFactors=FALSE)
dat <- subset(dat, dat$Gene %in% genes.sig)
dat <- subset(dat, select=c("Gene","pvalue"))
dat$nlogp <- -log10(dat$pvalue)

g <- c("ICOSLG","LAG3","CD274","PDCD1LG2","CD86","HAVCR2","PDCD1","CTLA4","ICOS","CD28","CD80","BTLA","TNFRSF14","HLA-C","LGALS9","HLA-A","HLA-B")
dat <- dat[match(g, dat$Gene),]

get.barplot <- function(dat, g){
	dat$Gene <- factor(dat$Gene, levels=g)

	# GENERATE PLOT ---
	p <- ggplot(dat, aes(x=Gene, y=nlogp)) + 
			geom_bar(fill="#000000", stat = "identity", color=NA, width=0.6) + 
			coord_flip() +
			theme(
				axis.text.x = element_text(size = 5, color="black"),
				axis.text.y = element_text(size = 5, color="black"),
				axis.title = element_text(size = 5, color="black"),
				plot.title = element_text(size = 8, color="black"),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.ticks = element_line(size=0.4, color="black"),	
				strip.text = element_text(size=6, color="black"),
				strip.background = element_rect(fill="white", color="white"),
				panel.background = element_rect(fill="white", color="black"),
				legend.text = element_text(size = 4.5, color="black"),
				legend.title = element_blank(),
				legend.key.size = unit(0.3, "cm"),
				legend.position = "none") +
			ylab("") +
			xlab("") + 
			ggtitle("") 
	return(p)	
}


# PLOT ---
file.plot <- file.path(dir.plot, "barplot_ImmuneCheckpointM_mRNA.pdf")
pdf(file.plot, width=1.5, height=2)
	p <- get.barplot(dat, g)
	grid.arrange(p, ncol=1, nrow=1)
dev.off()
