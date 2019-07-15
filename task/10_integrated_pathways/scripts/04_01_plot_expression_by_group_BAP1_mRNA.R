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
file.sig <- file.path(dir.data, "pathway_gene_relationship.tsv")
file.expr <- file.path(dir.expr, "meso_peritoneal_gene_expression_invnorm_proteincoding.tsv.gz")
file.diffexpr <- file.path("/Data/Raunak/softwares/bdvtools/array_process/wilcox_rank_test.R")

### LOAD SIGNATURE GENES ---
dat <- read.delim(file.sig, header=TRUE, stringsAsFactors=FALSE)
colnames(dat) <- c("Group","Gene")

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

#> unique(dat$Group)
# [1] "DNA_Repair"
# [2] "CellCycle:DNA_Repair"
# [3] "Chromatin_Remodeling"
# [4] "SWI/SNF_Complex:Chromatin_Remodeling"
# [5] "MAPK_Signaling"
# [6] "PI3K_Signaling"
# [7] "TOR_Signaling"
# [8] "Wnt_Signaling"
# [9] "Hippo_Signaling"
#[10] "NADPH_Metabolism"
#[11] "Uncategorized"



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



### GROUP: DNA_Repair ---
source(file.diffexpr)
expr.sig <- get.exprsig(expr, dat, group=c("DNA_Repair","CellCycle:DNA_Repair"))
df <- get.wilcox.rank.test(expr.sig, class1=grp0, class2=grp1)
dm <- dataByGroup(expr, df)
p <- get.boxplot(df=dm, plotitle="DNA_Repair", nrow=3, ncol=11)

file.plot <- file.path(dir.plot, "boxplot_mRNA_DNARepair.pdf")
pdf(file.plot, height=4, width=6.8)
	grid.arrange(p, ncol=1, nrow=1)
dev.off()


### GROUP: Chromatin_Remodeling ---
expr.sig <- get.exprsig(expr, dat, group=c("Chromatin_Remodeling","SWI/SNF_Complex:Chromatin_Remodeling"))
df <- get.wilcox.rank.test(expr.sig, class1=grp0, class2=grp1)
dm <- dataByGroup(expr, df)
p <- get.boxplot(df=dm, plotitle="Chromatin_Remodeling", nrow=2, ncol=11)

file.plot <- file.path(dir.plot, "boxplot_mRNA_ChromatinRemodeling.pdf")
pdf(file.plot, height=3, width=6.8)
	grid.arrange(p, ncol=1, nrow=1)
dev.off()


### GROUP: MAPK_Signaling ---
expr.sig <- get.exprsig(expr, dat, group="MAPK_Signaling")
df <- get.wilcox.rank.test(expr.sig, class1=grp0, class2=grp1)
dm <- dataByGroup(expr, df)
p <- get.boxplot(df=dm, plotitle="MAPK_Signaling", nrow=1, ncol=11)

file.plot <- file.path(dir.plot, "boxplot_mRNA_MAPKSignaling.pdf")
pdf(file.plot, height=2, width=6.8)
	grid.arrange(p, ncol=1, nrow=1)
dev.off()


### GROUP: PI3K_Signaling ---
expr.sig <- get.exprsig(expr, dat, group="PI3K_Signaling")
df <- get.wilcox.rank.test(expr.sig, class1=grp0, class2=grp1)
dm <- dataByGroup(expr, df)
p <- get.boxplot(df=dm, plotitle="PI3K_Signaling", nrow=1, ncol=11)

file.plot <- file.path(dir.plot, "boxplot_mRNA_PI3KSignaling.pdf")
pdf(file.plot, height=2, width=4)
	grid.arrange(p, ncol=1, nrow=1)
dev.off()


### GROUP: TOR_Signaling ---
expr.sig <- get.exprsig(expr, dat, group="TOR_Signaling")
df <- get.wilcox.rank.test(expr.sig, class1=grp0, class2=grp1)
dm <- dataByGroup(expr, df)
p <- get.boxplot(df=dm, plotitle="TOR_Signaling", nrow=1, ncol=11)

file.plot <- file.path(dir.plot, "boxplot_mRNA_TORSignaling.pdf")
pdf(file.plot, height=2, width=4)
	grid.arrange(p, ncol=1, nrow=1)
dev.off()


### GROUP: Uncategorized ---
expr.sig <- get.exprsig(expr, dat, group=c("Hippo_Signaling","NADPH_Metabolism","Uncategorized"))
df <- get.wilcox.rank.test(expr.sig, class1=grp0, class2=grp1)
dm <- dataByGroup(expr, df)
p <- get.boxplot(df=dm, plotitle="Uncategorized", nrow=2, ncol=11)

file.plot <- file.path(dir.plot, "boxplot_mRNA_Uncategorized.pdf")
pdf(file.plot, height=3, width=6.8)
	grid.arrange(p, ncol=1, nrow=1)
dev.off()


