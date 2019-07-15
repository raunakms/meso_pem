### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.prot <- file.path(dir.wrk, "task/06_proteomics/data")
dir.cna <- file.path(dir.wrk, "data/cnv/seq_call_refseq_genes_meso")
dir.analysis <- file.path(dir.wrk, "task/07_compare_mRNA_protein")
dir.output <- file.path(dir.analysis, "output")
dir.plot <- file.path(dir.analysis, "plot")
dir.script <- file.path(dir.analysis, "scripts")

### DEFINE FILES ---
file.prot <- file.path(dir.prot, "meso_expression_proteome_dqnorm_log2.tsv")
file.cna <- file.path(dir.cna, "meso_cnv_seg_values_calls_parsed.tsv")

### LOAD LIBRARIES ---
library("stringr")
library("VennDiagram")

### LOAD Protein expression ---
expr.prot <- read.delim(file.prot, header=T, stringsAsFactors=F, row.names=1)
colnames(expr.prot) <- str_replace_all(colnames(expr.prot), "[.]", "-")

### LOAD CNV ---
expr.cna <- read.delim(file.cna, header=T, stringsAsFactors=F, row.names=1)
colnames(expr.cna) <- str_replace_all(colnames(expr.cna), "[.]", "-")

### GET COMMON IDS ---
ids <- intersect(colnames(expr.prot), colnames(expr.cna))
expr.prot <- subset(expr.prot, select=ids)
expr.cna <- subset(expr.cna, select=ids)

### Compare mRNA & Protein genelist ---
genes.prot <- rownames(expr.prot)
genes.cna <- rownames(expr.cna)
list.genes <- list(protein=genes.prot, cnv=genes.cna)


# Venn Diagram ---
pdf(file.path(dir.plot, "venn_genelist_protein_cnv.pdf"))
p1 <- venn.diagram(x=list.genes, scaled = TRUE,
		filename = NULL, lwd=2.5, cex=1.5, fontfamily="Helvetica",
		height=2, width=2,  units="inches",
		col="black", fill=c("deepskyblue","cornsilk3"),
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
expr.prot <- subset(expr.prot, rownames(expr.prot) %in% genes.common)
expr.cna <- subset(expr.cna, rownames(expr.cna) %in% genes.common)

### Match Genes order ---
expr.prot <- expr.prot[match(genes.common, rownames(expr.prot)),]
expr.cna <- expr.cna[match(genes.common, rownames(expr.cna)),]

### GET FUNCTION ---
source(file.path(dir.script, "01_02_function_correlation.R"))

### CORRELATION OF TWO MATRICES ----
dats <- get.correlation(method="spearman", genes.common, expr.prot, expr.cna)
datp <- get.correlation(method="pearson", genes.common, expr.prot, expr.cna)

write.table(dats, file.path(dir.output, "correlation_summary_protein_cnv_spearman.tsv"), sep="\t", row.names=F, col.names=T, quote=F)
write.table(datp, file.path(dir.output, "correlation_summary_protein_cnv_pearson.tsv"), sep="\t", row.names=F, col.names=T, quote=F)




file.plot <- file.path(dir.plot, "expr_protein_cnv_correlation_rho_distribution.pdf")
pdf(file.plot, height=6, width=6)
par(mfrow=c(2,1))
	get.plot(dats, method="Spearman", analysis="protein_cnv")
	get.plot(datp, method="Pearson", analysis="protein_cnv")
dev.off()

