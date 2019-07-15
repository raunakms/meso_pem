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
file.des <- file.path(dir.output, "design_table_mrna_protein.tsv")
file.array <- file.path("/Data/Raunak/softwares/bdvtools/array_process/array_preprocess.R")

#### LOAD DESIGN TABLE ----------------
des <- read.delim(file.des, header=T, stringsAsFactors=F)

### LOAD mRNA expression ---
expr.mrna <- read.delim(file.mrna, header=T, stringsAsFactors=F, row.names=1)
colnames(expr.mrna) <- str_replace_all(colnames(expr.mrna), "[.]", "-")

### LOAD Protein expression ---
expr.prot <- read.delim(file.prot, header=T, stringsAsFactors=F, row.names=1)
colnames(expr.prot) <- str_replace_all(colnames(expr.prot), "[.]", "-")
expr.prot[is.na(expr.prot)] <- 0

### MATCH GENES ---
genes <- intersect(rownames(expr.mrna), rownames(expr.prot))
expr.mrna <- subset(expr.mrna, rownames(expr.mrna) %in% genes)
expr.prot <- subset(expr.prot, rownames(expr.prot) %in% genes)

expr.mrna <- expr.mrna[match(genes, rownames(expr.mrna)),]
expr.prot <- expr.prot[match(genes, rownames(expr.prot)),]

#> dim(expr.mrna)
#[1] 8109   15
#> dim(expr.prot)
#[1] 8109   26

y1 <- colnames(expr.mrna)
y2 <- colnames(expr.prot)

### PERFORM QUANTILE NORMALIZATION ---
expr.combo <- cbind(expr.mrna, expr.prot)

source(file.array)
dqnorm <- getQuantile.normalize(expr.combo)

### SEPARATE DATA ---
expr.mrna <- dqnorm[,1:15]
expr.prot <- dqnorm[,16:41]

colnames(expr.mrna) <- y1
colnames(expr.prot) <- y2


# MELT DATA --
df.mrna <- melt(t(as.matrix(expr.mrna)))
df.prot <- melt(t(as.matrix(expr.prot)))

colnames(df.mrna) <- c("SampleID","Gene","Expr")
colnames(df.prot) <- c("SampleID","Gene","Expr")

df.mrna$SampleID <- as.character(df.mrna$SampleID)
df.prot$Gene <- as.character(df.prot$Gene)


### FUNCTION: getCorrelation ---
getCorrelation <- function(df1, df2, s1, s2){
	# GET DATA ---
	d1 <- subset(df1, df1$SampleID == s1)
	d2 <- subset(df2, df2$SampleID == s2)

	# PREPARE DATA ---
	df <- data.frame(SampleID1=d1$SampleID,
					SampleID2=d2$SampleID,
					Gene=d1$Gene,
					Expr1=d1$Expr,
					Expr2=d2$Expr)
	df$Gene <- as.character(df$Gene)
	df$SampleID1 <- as.character(df$SampleID1)
	df$SampleID2 <- as.character(df$SampleID2)

	# FACTORIZE DATA ---
	df$Gene <- factor(df$Gene, levels=df$Gene)

	return(getCorrelationPlot(df))
}


### FUNCTION: getCorrelationPlot ---
getCorrelationPlot <- function(df){
	# LOAD LIBRARIES ---
	require("ggplot2")

	# GET IDS ---
	id1 <- unique(df$SampleID1)
	id2 <- unique(df$SampleID2)

	#GENERATE PLOT -----
	p <- ggplot(df, aes(x=Expr1, y=Expr2)) + 
			geom_point(color="black", alpha=0.5, size=0.1) +
			#scale_color_manual(values=c("#bdbdbd","#bd0026")) +
			#geom_text(aes(x=Expr.rnaseq, y=Expr.marray, label=Label), size=1, color="black", hjust=0, vjust=0)+
			#stat_smooth(method="lm", geom = "smooth", formula = y ~ x, position = "identity", fullrange = FALSE, se=FALSE) +
			geom_density2d(stat = "density2d", lineend = "round", linejoin = "round", alpha=0.9, color="yellow", size=0.25) +
			#coord_cartesian(xlim=c(0, 25), ylim=c(0, 25)) +
			theme(
				axis.text = element_text(size = 5, color="black"),
				axis.title = element_text(size = 8, color="black"),
				#strip.text = element_text(size = 10, color="black"),
				#strip.background = element_rect(fill = "white", colour = "white"),
				plot.title = element_text(size = 8, color="black", hjust=0),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.ticks = element_line(size=0.4, color = "black"),	
				panel.background = element_rect(fill = "white", colour = "black"),
				legend.position="none") +
			ylab("Expression Protein") + 
			xlab("Expression RNAseq") + 
			ggtitle(id1) 

	return(p)			
}



### CORRELATION BY SAMPLEPAIR ----
list.plot <- list()
for(i in 1:nrow(des)){
	s_mrna <- des$SampleID_RNA[i]
	s_prot <- des$SampleID_Protein[i]

	list.plot[[i]] <- getCorrelation(df.mrna, df.prot, s_mrna, s_prot)
}

### PLOT ---
file.plot <- file.path(dir.plot, "correlation_samplewise_mrna_protein.pdf")
pdf(file.plot, height=4, width=4)
	grid.arrange(list.plot[[1]], list.plot[[2]], list.plot[[3]], list.plot[[4]], ncol=2, nrow=2)
	grid.arrange(list.plot[[5]], list.plot[[6]], list.plot[[7]], list.plot[[8]], ncol=2, nrow=2)
	grid.arrange(list.plot[[9]], list.plot[[10]], list.plot[[11]], list.plot[[12]], ncol=2, nrow=2)
	grid.arrange(list.plot[[13]], ncol=2, nrow=2)
dev.off()

