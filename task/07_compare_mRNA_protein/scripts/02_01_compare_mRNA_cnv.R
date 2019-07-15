### LOAD LIBRARIES ---
library("stringr")
library("VennDiagram")
library("reshape2")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.wrk <- file.path("/collinsgroup/Raunak/projects/MESO_peritoneal")
dir.mrna <- file.path(dir.wrk, "data/expression/analysis")
dir.cna <- file.path(dir.wrk, "data/cnv/seq_call_refseq_genes_meso")
dir.analysis <- file.path(dir.wrk, "task/07_compare_mRNA_protein")
dir.output <- file.path(dir.analysis, "output")
dir.plot <- file.path(dir.analysis, "plot")
dir.script <- file.path(dir.analysis, "scripts")
dir.cgc <- file.path("/collinsgroup/Raunak/data_ref/cancer_gene_census")

### DEFINE FILES ---
file.mrna <- file.path(dir.mrna, "meso_peritoneal_gene_expression_proteincoding.tsv")
file.cna <- file.path(dir.cna, "meso_cnv_seg_values_calls_parsed.tsv")
file.ref <- file.path(dir.output, "mart_export.txt")
file.cgc <- file.path(dir.cgc, "Census_allFri_Nov_4_18_17_54_2016.tsv")

### LOAD CGC GENES ---
dat.cgc <- read.delim(file.cgc, header=T, stringsAsFactors=F)
genes.cgc <- unique(dat.cgc$Gene.Symbol)

### LOAD mRNA expression ---
expr.mrna <- read.delim(file.mrna, header=T, stringsAsFactors=F, row.names=1)
colnames(expr.mrna) <- str_replace_all(colnames(expr.mrna), "[.]", "-")

### LOAD CNV ---
expr.cna <- read.delim(file.cna, header=T, stringsAsFactors=F, row.names=1)
colnames(expr.cna) <- str_replace_all(colnames(expr.cna), "[.]", "-")

### GET COMMON IDS ---
ids <- intersect(colnames(expr.mrna), colnames(expr.cna))
expr.mrna <- subset(expr.mrna, select=ids)
expr.cna <- subset(expr.cna, select=ids)

### Compare mRNA & Protein genelist ---
genes.mrna <- rownames(expr.mrna)
genes.cna <- rownames(expr.cna)
list.genes <- list(mRNA=genes.mrna, cnv=genes.cna)


# Venn Diagram ---
pdf(file.path(dir.plot, "venn_genelist_mRNA_cnv.pdf"))
p1 <- venn.diagram(x=list.genes, scaled = TRUE,
		filename = NULL, lwd=2.5, cex=1.5, fontfamily="Helvetica",
		height=2, width=2,  units="inches",
		col="black", fill=c("firebrick1","cornsilk3"),
		main="", main.cex=1.5, main.fontfamily="Helvetica",
		cat.fontfamily="Helvetica", cat.cex=1, cat.default.pos="text",
		cat.pos=c(270, 0),
		#cat.dist=c(0.1,0.05), cat.just = list(c(0.5, 0.5), c(1,1))
		)
grid.draw(p1)
dev.off()


### Get common Genes -----	
genes.common <- Reduce(intersect, list.genes)

### Subset genes ---
expr.mrna <- subset(expr.mrna, rownames(expr.mrna) %in% genes.common)
expr.cna <- subset(expr.cna, rownames(expr.cna) %in% genes.common)

### Match Genes order ---
expr.mrna <- expr.mrna[match(genes.common, rownames(expr.mrna)),]
expr.cna <- expr.cna[match(genes.common, rownames(expr.cna)),]

### GET FUNCTION ---
source(file.path(dir.script, "01_02_function_correlation.R"))

### CORRELATION OF TWO MATRICES ----
#dats <- get.correlation(method="spearman", genes.common, expr.mrna, expr.cna)
datp <- get.correlation(method="pearson", genes.common, expr.mrna, expr.cna, genes.cgc)

write.table(dats, file.path(dir.output, "correlation_summary_mrna_cnv_spearman.tsv"), sep="\t", row.names=F, col.names=T, quote=F)
write.table(datp, file.path(dir.output, "correlation_summary_mrna_cnv_pearson.tsv"), sep="\t", row.names=F, col.names=T, quote=F)






### PLOT HISTOGRAM ---
file.plot <- file.path(dir.plot, "expr_mrna_cnv_correlation_rho_distribution.pdf")
pdf(file.plot, height=6, width=6)
par(mfrow=c(2,1))
	#get.plot(dats, method="Spearman", analysis="mrna_cnv")
	get.plot(datp, method="Pearson", analysis="mrna_cnv")
	get.cgcplot(datp, method="Pearson", analysis="")
dev.off()

write.table(rownames(expr.cna), file.path(dir.output, "genelist_mrna_cna_match.txt"), row.names=F, col.names=F, quote=F)

### STATS ---
get.corr.stat(df=datp, r.threshold=0.5)
get.corr.stat(df=datp, r.threshold=0.1)

#> get.corr.stat(df=datp, r.threshold=0.5)
#                Items    Value
#1         R.Threshold     0.50
#2         Genes.Total 19053.00
#3          Genes.Pass  3158.00
#4 Concordance.Percent    16.57
#> get.corr.stat(df=datp, r.threshold=0.1)
#                Items   Value
#1         R.Threshold     0.1
#2         Genes.Total 19053.0
#3          Genes.Pass 10575.0
#4 Concordance.Percent    55.5


### LOAD MART FILE ----
dat.ref <- read.delim(file.ref, header=TRUE, stringsAsFactors=FALSE)
dat.ref <- dat.ref[-which(str_detect(dat.ref$Chromosome, "CHR_") == TRUE),]
dat.ref <- subset(dat.ref, dat.ref$Gene %in% rownames(expr.cna))

#items <- paste("chr", c(1:22, "X", "Y"), sep="")
#dat.ref <- subset(dat.ref, dat.ref$chrom %in% items)

#### RESHAPE DATA ---
df.mrna <- melt(t(as.matrix(expr.mrna)))
df.cna <- melt(t(as.matrix(expr.cna)))
colnames(df.mrna) <- c("SampleID","Gene","value")
colnames(df.cna) <- c("SampleID","Gene","value")

df.mrna$Type <- "mRNA"
df.cna$Type <- "CNA"

df <- data.frame(SampleID=df.mrna$SampleID,
				Gene=df.mrna$Gene,
				mRNA=df.mrna$value,
				CNA=df.cna$value)
df$SampleID <- as.character(df$SampleID)
df$Gene <- as.character(df$Gene)
df <- subset(df, df$Gene %in% dat.ref$Gene)

### ADD CHROMOSOME INFO ---
df$chrom <- ""
for(i in 1:nrow(dat.ref)){
	index <- which(df$Gene == dat.ref$Gene[i])
	df$chrom[index] <- dat.ref$Chromosome[i]
}

### WRITE OUTPUT ---
file.output <- file.path(dir.output, "meso_pem_mrna_cnv_data_for_plot.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)

df <- read.delim(file.output, header=T, stringsAsFactors=F)

###### GENERATE PLOT ---------------
p <- ggplot(df, aes(x=CNA, y=mRNA)) + 
		#geom_point(shape = 21, fill="white", color="black", size=2, stroke=1, alpha=0.9) +
		geom_point(aes(color=chrom), size=1, alpha=0.7) +
		#scale_color_manual(values=c("green","gray50")) +
		stat_smooth(method="lm", geom = "smooth", formula = y ~ x, position = "identity", fullrange = FALSE, se=TRUE) +
		#geom_density2d(stat = "density2d", lineend = "round", linejoin = "round", alpha=0.9, color="red", size=0.25) +
		#coord_cartesian(xlim=c(-0.5, 18), ylim=c(1, 10)) +
		theme(
			axis.text = element_text(size = 10, color="black"),
			axis.title = element_text(size = 15, color="black"),
			strip.text = element_text(size = 10, color="black"),
			strip.background = element_rect(fill = "white", colour = "white"),
			plot.title = element_text(size = 10, color="black", hjust=0),
			#panel.grid.major = element_line(size=0.6, color = "grey85"),
			#panel.grid.minor = element_line(size=0.4, color = "grey85"),
			#axis.ticks = element_blank(),	
			#panel.background = element_rect(fill = "white", colour = "white"),
			legend.position="none") +
		ylab("mRNA expression (log2)") + 
		xlab("Copy Number Segment Mean (log2)") + 
		ggtitle("") 	

### PLOT ---
file.plot <- file.path(dir.plot, "meso_mrna_cna_correlation_scatterplot.pdf")
pdf(file.plot, height=6, width=6)
	grid.arrange(p, ncol=1, nrow=1)
dev.off()


### PLOT: SCATTER PLOT ---
x.mrna <- as.numeric(as.matrix(expr.mrna))
y.cna <- as.numeric(as.matrix(expr.cna))
y.cna <- log2(10 - y.cna)

pdf(file.path(dir.plot, "expr_mrna_cna_correlation_scatterplot.pdf"), width=6, height=6)
	plot(x=x.mrna, y=y.cna, xlab="mRNA Expression", ylab="Copy Number Segment Mean", 
		tck=-0.03, cex=0.2, pch=19, cex.lab=0.8, cex.axis=0.6, cex.main=0.8,
		main="MESO: Correlation - mRNA Expression vs Copy Number")
	abline(lm(y.cna~x.mrna), col="red", lwd=2) # regression line (y~x) 
	#text(x=12,y=1,labels="R = 0.28", cex=1.2, col="white")
dev.off()

ct.s <- cor.test(x.mrna, y.prot, method="spearman")
ct.p <- cor.test(x.mrna, y.prot, method="pearson")
