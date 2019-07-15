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
file.sig <- file.path(dir.data, "genesets_Chromatin_Modification.tsv")
file.expr <- file.path(dir.expr, "meso_peritoneal_gene_expression_invnorm_proteincoding.tsv.gz")
file.diffexpr <- file.path("/Data/Raunak/softwares/bdvtools/array_process/wilcox_rank_test.R")

### LOAD SIGNATURE GENES ---
dat.gmt <- read.delim(file.sig, header=TRUE, stringsAsFactors=FALSE)

### LOAD DATA ---
expr <- read.delim(file.expr, header=TRUE, stringsAsFactors=FALSE, row.names=1)
colnames(expr) <- str_replace_all(colnames(expr), "[.]", "-")
expr[is.na(expr)] <- 0

### LOAD DESIGN TABLE ---
des <- read.delim(file.des, header=TRUE, stringsAsFactors=FALSE)
grp0 <- intersect(colnames(expr), des$Sample.ID[which(des$Group.3p21 == 0)])
grp1 <- intersect(colnames(expr), des$Sample.ID[which(des$Group.3p21 == 1)])

### ARRANGE SAMPLES ---
expr <- subset(expr, select=c(grp0, grp1))




### FUNCTION: parseGMTdata ---
parseGMTdata <- function(dat){
	require("stringr")
	list.dat <- list()
	for(i in 1:nrow(dat)){
		group <- dat$Category[i]
		genes <- str_split(dat$Genesets[i], ":")[[1]]

		list.dat[[i]] <- data.frame(Group=group, Gene=genes)
	}

	df <- do.call(rbind.data.frame, list.dat)
	colnames(df) <- c("Group","Gene")
	df$Group <- as.character(df$Group)
	df$Gene <- as.character(df$Gene)	
	return(df)
}

### FUNCTION: get.exprsig ----
get.exprsig <- function(expr, dat, group){
	# SUBSET BY GROUP ---
	df <- subset(dat, dat$Group %in% group)
	genes.sig <- df$Gene

	# SUBSET EXPRESSION ---
	expr.sig <- subset(expr, rownames(expr) %in% genes.sig)

	return(expr.sig)
}

### FUNCTION: ----
dataByGroup <- function(expr, df){
	require("reshape2")

	# SUBSET BY GENES ---
	#df <- subset(df, df$pvalue < 0.1)
	gexpr <- subset(expr, rownames(expr) %in% df$Gene)

	# MELT DATA ---
	dm <- melt(as.matrix(t(gexpr)))
	colnames(dm) <- c("SampleID","Gene","Value")
	dm$SampleID <- as.character(dm$SampleID)
	dm$Gene <- factor(dm$Gene, levels=df$Gene)
	dm$Group <- ""

	dm$Group[which(dm$SampleID %in% colnames(gexpr)[1:7])] <- "BAP1-Intact"
	dm$Group[which(dm$SampleID %in% colnames(gexpr)[8:15])] <- "BAP1-Del"

	dm$Group <- factor(dm$Group, c("BAP1-Intact","BAP1-Del"))

	return(dm)
}

### FUNCTION: ----
get.boxplot <- function(df, plotitle, nrow, ncol){
	require("ggplot2")
	cbPalette <- c("#fffc00","blue")

	# PLOT ---
	p <- ggplot(df, aes(x=Group, y=Value)) +
			geom_boxplot(aes(fill=Group, alpha=0.3), lwd=0.3, color="black", outlier.size=0.1, outlier.alpha=0.3, notch=FALSE) +
			scale_fill_manual(values=cbPalette) +
			#geom_jitter(width=0.1, color="black", alpha=0.5, size=0.5) +
			facet_wrap(~ Gene, nrow=nrow, ncol=ncol) +
			theme(
				axis.text.x = element_text(size=7, color="black", angle=90, hjust=0, vjust=0),
				axis.text.y = element_text(size=7, color="black"),
				axis.title = element_text(size=7, color="black"),
				plot.title = element_text(size=10, color="black", hjust=0),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),			
				axis.ticks = element_line(size=0.4, color="black"),
				strip.text = element_text(size=6, color="black"),
				strip.background = element_rect(fill="white", color="white"),
				panel.background = element_rect(fill="white", color="black"),
				legend.position="none") + 
			ylab("Expression") +
			xlab("") + ggtitle(plotitle) 

	return(p)	
}

### FUNCTION: main.func ---
main.func <- function(file.diffexpr, expr, dat, grp0, grp1, path.group, nrow, ncol){
	source(file.diffexpr)
	expr.sig <- get.exprsig(expr, dat, group=path.group)
	df <- get.wilcox.rank.test(expr.sig, class1=grp0, class2=grp1)
	dm <- dataByGroup(expr, df)
	p <- get.boxplot(df=dm, plotitle=path.group, nrow=nrow, ncol=ncol)
	return(p)
}


###
> data.frame(Category=dat$Category, Genes=unlist(lapply(str_split(dat$Genesets, ":"), function(x) length(x))))
                                                 Category Genes
1                   GO_ATP_DEPENDENT_CHROMATIN_REMODELING    74
2                    GO_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY   177
3                                GO_CHROMATIN_DISASSEMBLY    17
4                               GO_CHROMATIN_MODIFICATION   539
5                                 GO_CHROMATIN_REMODELING   150
6                      GO_COVALENT_CHROMATIN_MODIFICATION   345
7    GO_DNA_REPLICATION_DEPENDENT_NUCLEOSOME_ORGANIZATION    32
8  GO_DNA_REPLICATION_INDEPENDENT_NUCLEOSOME_ORGANIZATION    53
9        GO_POSITIVE_REGULATION_OF_CHROMATIN_MODIFICATION    85
10                     GO_REGULATION_OF_CHROMATIN_BINDING    17
11       GO_NEGATIVE_REGULATION_OF_CHROMATIN_MODIFICATION    45
12                   GO_REGULATION_OF_CHROMATIN_SILENCING    21


### GROUP: DNA_Repair 1 ---
dat.agg <- parseGMTdata(dat=dat.gmt[-c(2,4,5,6),])
list.plot <- list()
	#for(k in 1:2){
for(k in 1:length(unique(dat.agg$Group))){
	dat <- subset(dat.agg, dat.agg$Group == unique(dat.agg$Group)[k])
	list.plot[[k]] <- main.func(file.diffexpr, expr, dat, grp0, grp1, path.group=unique(dat.agg$Group)[k], nrow=8, ncol=11)
}


file.plot <- file.path(dir.plot, "boxplot_mRNA_ChromatinModification_genesets_01.pdf")
pdf(file.plot, height=7, width=6.8)
	for(k in 1:length(list.plot)){
		grid.arrange(list.plot[[k]], ncol=1, nrow=1)
	}
dev.off()





### GROUP: DNA_Repair 2 ---
dat.agg <- parseGMTdata(dat=dat.gmt[5,])
list.plot <- list()
	#for(k in 1:2){
for(k in 1:length(unique(dat.agg$Group))){
	dat <- subset(dat.agg, dat.agg$Group == unique(dat.agg$Group)[k])
	list.plot[[k]] <- main.func(file.diffexpr, expr, dat, grp0, grp1, path.group=unique(dat.agg$Group)[k], nrow=10, ncol=15)
}


file.plot <- file.path(dir.plot, "boxplot_mRNA_ChromatinModification_genesets_02.pdf")
pdf(file.plot, height=10, width=8)
	for(k in 1:length(list.plot)){
		grid.arrange(list.plot[[k]], ncol=1, nrow=1)
	}
dev.off()




### COSTOM GENES ---
dat.agg <- data.frame(Group="CHROMATIN_REMODELING", Gene=c("BAP1","PBRM1","SETD2","SMARCC1","EP300"))
dat.agg$Group <- as.character(dat.agg$Group)
dat.agg$Gene <- as.character(dat.agg$Gene)

list.plot <- list()
for(k in 1:length(unique(dat.agg$Group))){
	dat <- subset(dat.agg, dat.agg$Group == unique(dat.agg$Group)[k])
	list.plot[[k]] <- main.func(file.diffexpr, expr, dat, grp0, grp1, path.group=unique(dat.agg$Group)[k], nrow=1, ncol=5)
}


file.plot <- file.path(dir.plot, "boxplot_mRNA_ChromatinRemodeling_selected.pdf")
pdf(file.plot, height=2, width=3)
	for(k in 1:length(list.plot)){
		grid.arrange(list.plot[[k]], ncol=1, nrow=1)
	}
dev.off()

