### LOAD LIBRARIES ---
library("stringr")
library("VennDiagram")
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

### Compare mRNA & Protein genelist ---
genes.mrna <- rownames(expr.mrna)
genes.prot <- rownames(expr.prot)
list.genes <- list(mRNA=genes.mrna, protein=genes.prot)

# Venn Diagram ---
pdf(file.path(dir.plot, "venn_genelist_mRNA_protein.pdf"))
p1 <- venn.diagram(x=list.genes, scaled = TRUE,
		filename = NULL, lwd=2.5, cex=1.5, fontfamily="Helvetica",
		height=2, width=2,  units="inches",
		col="black", fill=c("firebrick1","deepskyblue"),
		main="", main.cex=1.5, main.fontfamily="Helvetica",
		cat.fontfamily="Helvetica", cat.cex=1, cat.default.pos="text",
		cat.pos=c(0, 270),
		#cat.dist=c(0.1,0.05), cat.just = list(c(0.5, 0.5), c(1,1))
		)
grid.draw(p1)
dev.off()

### Get common Genes -----	
genes.common <- Reduce(intersect, list.genes)

### Subset genes ---
expr.mrna <- subset(expr.mrna, rownames(expr.mrna) %in% genes.common)
expr.prot <- subset(expr.prot, rownames(expr.prot) %in% genes.common)

### Match Genes order ---
expr.mrna <- expr.mrna[match(genes.common, rownames(expr.mrna)),]
expr.prot <- expr.prot[match(genes.common, rownames(expr.prot)),]

### GET FUNCTION ---
source(file.path(dir.script, "01_02_function_correlation.R"))

### CORRELATION OF TWO MATRICES ----
dats <- get.correlation(method="spearman", genes.common, expr.mrna, expr.prot)
datp <- get.correlation(method="pearson", genes.common, expr.mrna, expr.prot)

write.table(dats, file.path(dir.output, "correlation_summary_mrna_protein_spearman.tsv"), sep="\t", row.names=F, col.names=T, quote=F)
write.table(datp, file.path(dir.output, "correlation_summary_mrna_protein_pearson.tsv"), sep="\t", row.names=F, col.names=T, quote=F)


#dat.cor.pass <- subset(dat.cor, dat.cor$R >= 0.75)
#dat.cor.pass <- dat.cor.pass[order(dat.cor.pass$R, decreasing=T),]

file.plot <- file.path(dir.plot, "expr_mrna_prot_correlation_rho_distribution.pdf")
pdf(file.plot, height=6, width=6)
par(mfrow=c(2,1))
	get.plot(dats, method="Spearman", analysis="mrna_protein")
	get.plot(datp, method="Pearson", analysis="mrna_protein")
dev.off()


### STATS ---
get.corr.stat(df=datp, r.threshold=0.5, direction="POSITIVE")
get.corr.stat(df=datp, r.threshold=0.1, direction="POSITIVE")

#> get.corr.stat(df=datp, r.threshold=0.5)
#                Items   Value
#1         R.Threshold    0.50
#2         Genes.Total 8109.00
#3          Genes.Pass 1839.00
#4 Concordance.Percent   22.68
#> get.corr.stat(df=datp, r.threshold=0.1)
#                Items   Value
#1         R.Threshold    0.10
#2         Genes.Total 8109.00
#3          Genes.Pass 4715.00
#4 Concordance.Percent   58.15


get.corr.stat(df=datp, r.threshold=-0.5, direction="NEGATIVE")
#> get.corr.stat(df=datp, r.threshold=-0.5, direction="NEGATIVE")
#                Items    Value
#1         R.Threshold     -0.5
#2           Direction NEGATIVE
#3         Genes.Total     8109
#4          Genes.Pass      194
#5 Concordance.Percent     2.39



df.neg <- subset(datp, datp$R <= -0.5)
write.table(df.neg, file.path(dir.output, "correlation_summary_mrna_protein_pearson_negative.tsv"), sep="\t", row.names=F, col.names=T, quote=F)
write.table(df.neg$Gene, file.path(dir.output, "correlation_summary_mrna_protein_pearson_negative_genelist.txt"), sep="\t", row.names=F, col.names=F, quote=F)


### PLOT: SCATTER PLOT ---
x.mrna <- as.numeric(as.matrix(expr.mrna))
y.prot <- as.numeric(as.matrix(expr.prot))

pdf(file.path(dir.plot, "expr_mrna_prot_correlation_scatterplot.pdf"), width=6, height=6)
	plot(x=x.mrna, y=y.prot, xlab="mRNA Expression", ylab="Protein Expression", 
		tck=-0.03, cex=0.2, pch=19, cex.lab=0.8, cex.axis=0.6, cex.main=0.8,
		main="MESO: Correlation - mRNA (RNAseq) vs Protein (MassSpec) Expression")
	abline(lm(y.prot~x.mrna), col="red", lwd=2) # regression line (y~x) 
	text(x=12,y=1,labels="R = 0.28", cex=1.2, col="white")
dev.off()

ct.s <- cor.test(x.mrna, y.prot, method="spearman")
ct.p <- cor.test(x.mrna, y.prot, method="pearson")
