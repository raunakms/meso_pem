### LOAD LIBRARIES ---
library("stringr")
library("reshape2")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.mrna <- file.path(dir.wrk, "data/expression/analysis")
dir.prot <- file.path(dir.wrk, "data/proteome/processed_data")
dir.analysis <- file.path(dir.wrk, "task/07_compare_mRNA_protein")
dir.output <- file.path(dir.analysis, "output")
dir.plot <- file.path(dir.analysis, "plot")
dir.script <- file.path(dir.analysis, "scripts")

### DEFINE FILES ---
file.mrna <- file.path(dir.mrna, "meso_peritoneal_gene_expression_proteincoding.tsv.gz")
file.prot <- file.path(dir.prot, "meso_proteome_proteomediscover_all_log2_dqnorm.tsv.gz")
file.des <- file.path(dir.wrk, "data/annotation/design_table_3p21genes.tsv")

### LOAD DESIGN TABLE: mRNA ---
des <- read.delim(file.des, header=TRUE, stringsAsFactors=FALSE)
grp0 <- des$Sample.ID[which(des$Group.3p21 == 0)]
grp1 <- des$Sample.ID[which(des$Group.3p21 == 1)]

### LOAD mRNA expression ---
expr.mrna <- read.delim(file.mrna, header=T, stringsAsFactors=F, row.names=1)
colnames(expr.mrna) <- str_replace_all(colnames(expr.mrna), "[.]", "-")

### LOAD Protein expression ---
expr.prot <- read.delim(file.prot, header=T, stringsAsFactors=F, row.names=1)[,-c(1:7)]
colnames(expr.prot) <- str_replace_all(colnames(expr.prot), "[.]", "-")
colnames(expr.prot) <- str_replace_all(colnames(expr.prot), "T", "")

### GET COMMON IDS ---
ids <- intersect(colnames(expr.mrna), colnames(expr.prot))
expr.mrna <- subset(expr.mrna, select=ids)
expr.prot <- subset(expr.prot, select=ids)


### SELECT GENES ---
genes <- c("ARID1B","SMARCA2","HIST1H1A","EGFR","ARID1A","SMARCAD1",
				"PSMD8","LATS2","DNAH14")

#### SUBSET DATA ---
expr.mrna <- subset(expr.mrna, rownames(expr.mrna) %in% genes)
expr.prot <- subset(expr.prot, rownames(expr.prot) %in% genes)

expr.mrna <- expr.mrna[match(genes, rownames(expr.mrna)),]
expr.prot <- expr.prot[match(genes, rownames(expr.prot)),]

### MELT DATA ---
df1 <- melt(t(as.matrix(expr.mrna)))
df2 <- melt(t(as.matrix(expr.prot)))
colnames(df1) <- c("SampleID", "Gene", "Expr.mrna")
colnames(df2) <- c("SampleID", "Gene", "Expr.prot")
df2$Expr.prot[which(is.na(df2$Expr.prot))] <- 0

df <- df1
df$Expr.prot <- df2$Expr.prot
df$Gene <- factor(df$Gene, levels=genes)

### ADD GROUP ---
df$SampleID <- as.character(df$SampleID)
df$Group <- ""
df$Group[which(df$SampleID %in% grp0)] <- "BAP1-INTACT"
df$Group[which(df$SampleID %in% grp1)] <- "BAP1-DEL"
df$Group <- factor(df$Group, levels=c("BAP1-INTACT","BAP1-DEL"))

### FUNCTION: getCorrelationPlot ---
getCorrelationPlot <- function(df){
	# LOAD LIBRARIES ---
	require("ggplot2")
	cbPalette <- c("#fffc00","blue")

	#GENERATE PLOT -----
	p <- ggplot(df, aes(x=Expr.mrna, y=Expr.prot)) + 
			scale_fill_manual(values=cbPalette) +
			stat_smooth(method="lm", geom = "smooth", formula = y ~ x, color="red", position = "identity", fullrange = FALSE, se=FALSE) +
			geom_point(aes(fill=Group), color="black", shape=21, size=3, alpha=0.9) +
			#coord_cartesian(ylim=c(-5, 3)) +
			facet_wrap(~Gene, nrow=3, ncol=3, scales="free", drop=TRUE) +
				theme(
					axis.text.x = element_text(size=7, color="black"),
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
			ylab("Protein") + 
			xlab("mRNA") + 
			ggtitle("") 

	return(p)			
}

### GENRATE PLOT ---
file.plot <- file.path(dir.plot, "correlation_mRNA_PROTEIN_selected_genes.pdf")
pdf(file.plot, height=6, width=6)
		grid.arrange(getCorrelationPlot(df), ncol=1, nrow=1)
dev.off()
