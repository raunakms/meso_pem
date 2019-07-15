### Load Libraries -----
library("stringr")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.analysis <- file.path(dir.wrk, "task/07_compare_mRNA_protein")
dir.output <- file.path(dir.analysis, "output")
dir.plot <- file.path(dir.analysis, "plot")
dir.script <- file.path(dir.analysis, "scripts")
dir.cna <- file.path(dir.wrk, "data/cnv/seq_call_refseq_genes_meso")
dir.mrna <- file.path(dir.wrk, "data/expression/analysis")
dir.prot <- file.path(dir.wrk, "data/proteome/processed_data")

### DEFINE FILES ---
file.des.mrna <- file.path(dir.wrk, "data/annotation/design_table_3p21genes.tsv")
file.des.prot <- file.path(dir.wrk, "data/proteome/processed_data/design_table_proteome.tsv")
file.cna <- file.path(dir.cna, "meso_cnv_seg_values_calls_parsed.tsv")
file.mrna <- file.path(dir.mrna, "meso_peritoneal_gene_expression_proteincoding.tsv.gz")
file.prot <- file.path(dir.prot, "meso_proteome_proteomediscover_all_log2_dqnorm.tsv.gz")

### LOAD CNV DATA ---
dat.cna <- read.delim(file.cna, header=T, stringsAsFactors=F, row.names=1)
colnames(dat.cna) <- str_replace_all(colnames(dat.cna), "[.]", "-")

### LOAD DESIGN TABLE: CNV ---
des.mrna <- read.delim(file.des.mrna, header=TRUE, stringsAsFactors=FALSE)
grp0 <- intersect(colnames(dat.cna), des.mrna$Sample.ID[which(des.mrna$Group.3p21 == 0)])
grp1 <- intersect(colnames(dat.cna), des.mrna$Sample.ID[which(des.mrna$Group.3p21 == 1)])

### ARRANGE SAMPLES: CNV ---
dat.cna <- subset(dat.cna, select=c(grp0, grp1))

### LOAD DATA: mRNA ---
expr.mrna <- read.delim(file.mrna, header=TRUE, stringsAsFactors=FALSE, row.names=1)
colnames(expr.mrna) <- str_replace_all(colnames(expr.mrna), "[.]", "-")
expr.mrna[is.na(expr.mrna)] <- 0

### LOAD DESIGN TABLE: mRNA ---
des.mrna <- read.delim(file.des.mrna, header=TRUE, stringsAsFactors=FALSE)
grp0 <- intersect(colnames(expr.mrna), des.mrna$Sample.ID[which(des.mrna$Group.3p21 == 0)])
grp1 <- intersect(colnames(expr.mrna), des.mrna$Sample.ID[which(des.mrna$Group.3p21 == 1)])

### ARRANGE SAMPLES: mRNA ---
expr.mrna <- subset(expr.mrna, select=c(grp0, grp1))


### LOAD DATA: Protein ---
expr.prot <- read.delim(file.prot, header=TRUE, stringsAsFactors=FALSE, row.names=1)
colnames(expr.prot) <- str_replace_all(colnames(expr.prot), "[.]", "-")
expr.prot[is.na(expr.prot)] <- 0

### LOAD DESIGN TABLE: Protein ---
des.prot <- read.delim(file.des.prot, header=T, stringsAsFactors=F)
des.prot <- subset(des.prot, des.prot$SampleType == "Tumor")
grp0 <- des.prot$SequencingID[which(des.prot$Group.3p21 == 0)]
grp1 <- des.prot$SequencingID[which(des.prot$Group.3p21 == 1)]

### ARRANGE EXPRESSION DATA ---
expr.prot <- subset(expr.prot, select=c(grp0, grp1))
colnames(expr.prot) <- str_replace_all(colnames(expr.prot), "T", "")

> colnames(dat.cna)
 [1] "MESO-11"  "MESO-15"  "MESO-18A" "MESO-18E" "MESO-01"  "MESO-08"
 [7] "MESO-04"  "MESO-19"  "MESO-12"  "MESO-13"  "MESO-03"  "MESO-06"
[13] "MESO-17"  "MESO-07"  "MESO-09"  "MESO-14"  "MESO-05"  "MESO-10"
[19] "MESO-02"
> colnames(expr.mrna)
 [1] "MESO-08"  "MESO-11"  "MESO-12"  "MESO-13"  "MESO-19"  "MESO-18A"
 [7] "MESO-18E" "MESO-02"  "MESO-09"  "MESO-10"  "MESO-14"  "MESO-17"
[13] "MESO-05"  "MESO-06"  "MESO-07"
> colnames(expr.prot)
 [1] "MESO-01"  "MESO-03"  "MESO-04"  "MESO-11"  "MESO-12"  "MESO-13"
 [7] "MESO-18A" "MESO-18E" "MESO-02"  "MESO-05"  "MESO-06"  "MESO-07"
[13] "MESO-09"  "MESO-10"  "MESO-14"  "MESO-17"


### LOAD DESIGN TABLE ---
sampleids <- unique(c(colnames(dat.cna), colnames(expr.mrna), colnames(expr.prot)))
des <- read.delim(file.des.mrna, header=TRUE, stringsAsFactors=FALSE)
des <- subset(des, des$Sample.ID %in% sampleids)
des <- des[order(des$Group.3p21, des$Sample.ID),]
#grp0 <- des$Sample.ID[which(des$Group.3p21 == 0)]
#grp1 <- des$Sample.ID[which(des$Group.3p21 == 1)]

### SELECT GENES ---
genes <- c("ARID1A","SMARCC1","ARID1B","PBRM1","PRMT1","SETD2","SETDB1","EP300")

### FUNCTION: ---
prepareData <- function(dat, des, genes, type){
	require("reshape2")

	# SUBSET AND RESHAPE DATA  ----
	dat <- subset(dat, rownames(dat) %in% genes)
	df <- melt(t(as.matrix(dat)))
	colnames(df) <- c("SampleID","Gene","Value")
	df$SampleID <- as.character(df$SampleID)
	df$Gene <- as.character(df$Gene)
	df$Group <- ""
	df$Type <- type

	# GET GROUPS ---
	grp0 <- intersect(colnames(dat), des$Sample.ID[which(des$Group.3p21 == 0)])
	grp1 <- intersect(colnames(dat), des$Sample.ID[which(des$Group.3p21 == 1)])
	df$Group[which(df$SampleID %in% grp0)] <- "BAP1-INTACT"
	df$Group[which(df$SampleID %in% grp1)] <- "BAP1-DEL"

	return(df)
}

### CALL FUNCTION: ---
df.cna <- prepareData(dat=dat.cna, des=des, genes=genes, type="CNA")
df.mrna <- prepareData(dat=expr.mrna, des=des, genes=genes, type="mRNA")
df.prot <- prepareData(dat=expr.prot, des=des, genes=genes, type="Protein")

### COMBINE DATA ---
df <- rbind(df.cna, df.mrna, df.prot)
df$SampleID <- factor(df$SampleID, levels=des$Sample.ID)
df$Gene <- factor(df$Gene, levels=genes)
df$Group <- factor(df$Group, levels=c("BAP1-INTACT","BAP1-DEL"))
df$Type <- factor(df$Type, levels=c("CNA","mRNA","Protein"))

### FUNCTION: ----
get.plot <- function(df){
	require("ggplot2")
	cbPalette <- c("#fffc00","blue")

	# PLOT ---
	p <- ggplot(df, aes(x=SampleID, y=Value)) +
			geom_bar(aes(fill=Group, alpha=0.5), stat = "identity", color="black") +
			scale_fill_manual(values=cbPalette) +
			facet_grid(Type ~ Gene, scales="free_y", drop=FALSE) +
			theme(
				axis.text.x = element_text(size=7, color="black", angle=90, hjust=0, vjust=0),
				axis.text.y = element_text(size=7, color="black"),
				axis.title = element_text(size=7, color="black"),
				plot.title = element_text(size=10, color="black", hjust=0),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),			
				axis.ticks = element_line(size=0.4, color="black"),
				strip.text = element_text(size=10, color="black"),
				strip.background = element_rect(fill="white", color="white"),
				panel.background = element_rect(fill="white", color="black"),
				legend.position="none") + 
			ylab("") +
			xlab("") + ggtitle("") 

	return(p)	
}

### PLOT ---
file.plot <- file.path(dir.plot, "plot_bar_gtp.pdf")
pdf(file.plot, height=6, width=20)
	grid.arrange(get.plot(df), ncol=1, nrow=1)
dev.off()


