### LOAD LIBRAIES ---
library("ggplot2")
library("gridExtra")
library("VennDiagram")


### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/07_compare_mRNA_protein")
dir.output <- file.path(dir.wrk, "output")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")

### DEFINE FILES ---
file.mrna_cnv <- file.path(dir.output, "correlation_summary_mrna_cnv_spearman.tsv")
file.mrna_protein <- file.path(dir.output, "correlation_summary_mrna_protein_spearman.tsv")
file.protein_cnv <- file.path(dir.output, "correlation_summary_protein_cnv_spearman.tsv")

### LOAD DATA ---
dat.mrna_cnv <-  read.delim(file.mrna_cnv, header=T, stringsAsFactors=F)
dat.mrna_protein <-  read.delim(file.mrna_protein, header=T, stringsAsFactors=F)
dat.protein_cnv <-  read.delim(file.protein_cnv, header=T, stringsAsFactors=F)

### Prepare Stats Data ---
df.mrna_cnv <- data.frame(Analysis="mRNA-CNV", Group=names(table(dat.mrna_cnv$Group)), nGenes=as.numeric(table(dat.mrna_cnv$Group)))
df.mrna_protein <- data.frame(Analysis="mRNA-Protein", Group=names(table(dat.mrna_protein$Group)), nGenes=as.numeric(table(dat.mrna_protein$Group)))
df.protein_cnv <- data.frame(Analysis="Protein-CNV", Group=names(table(dat.protein_cnv$Group)), nGenes=as.numeric(table(dat.protein_cnv$Group)))

df.mrna_cnv$Proportion <- df.mrna_cnv$nGenes/nrow(dat.mrna_cnv)
df.mrna_protein$Proportion <- df.mrna_protein$nGenes/nrow(dat.mrna_protein)
df.protein_cnv$Proportion <- df.protein_cnv$nGenes/nrow(dat.protein_cnv)

df <- rbind(df.mrna_cnv, df.mrna_protein, df.protein_cnv)

### Factorize data ---
df$Analysis <- factor(df$Analysis, levels=c("mRNA-CNV","mRNA-Protein","Protein-CNV"))
df$Group <- factor(df$Group, levels=rev(c("Positive","Modrate","Weak","Negative")))

### Generate Barplot ---
p1 <- ggplot(df, aes(y=nGenes, x=Group)) +
		geom_bar(stat="Identity") +
		facet_wrap(~ Analysis, ncol=1) +
		theme(
			axis.text.x = element_text(size = 6, color="black"),
			axis.text.y = element_text(size = 6, color="black"),
			axis.title = element_text(size = 6, color="black"),
			plot.title = element_text(size = 6, color="black"),
			axis.ticks = element_line(size=0.4),
			strip.text = element_text(size = 6, color="black"),
			legend.position = "none") +
		ylab("No. of Genes") +
		xlab("Correlation Type") + 
		ggtitle("")

p2 <- ggplot(df, aes(y=Proportion, x=Group)) +
		geom_bar(stat="Identity") +
		facet_wrap(~ Analysis, ncol=1) +
		theme(
			axis.text.x = element_text(size = 6, color="black"),
			axis.text.y = element_text(size = 6, color="black"),
			axis.title = element_text(size = 6, color="black"),
			plot.title = element_text(size = 6, color="black"),
			axis.ticks = element_line(size=0.4),
			strip.text = element_text(size = 6, color="black"),
			legend.position = "none") +
		ylab("No. of Genes") +
		xlab("Correlation Type") + 
		ggtitle("")

file.plot <- file.path(dir.plot, "correlation_stat.pdf")
pdf(file.plot, width=5, height=4)		
	grid.arrange(p1, p2, nrow=1)
dev.off()

### Get Genes ---
get.genes <- function(dat, grp){
	genes <- dat$Gene[which(dat$Group == grp)]
	return(genes)
}

genes.pos.mrna_cnv <- get.genes(dat=dat.mrna_cnv, grp="Positive")
genes.pos.mrna_protein <- get.genes(dat=dat.mrna_protein, grp="Positive")

write.table(genes.pos.mrna_cnv, file.path(dir.output, "genelist_positive_mRNA-CNV.txt"), row.names=F, col.names=F, quote=F)
write.table(genes.pos.mrna_protein, file.path(dir.output, "genelist_positive_mRNA-Protein.txt"), row.names=F, col.names=F, quote=F)

list.genes <- list(mRNA_CNV=genes.pos.mrna_cnv, mRNA_Protein=genes.pos.mrna_protein)

### Get common Genes -----	
genes.common <- Reduce(intersect, list.genes)
write.table(genes.common, file.path(dir.output, "genelist_positive_common_mRNA-CNV_mRNA-Protein.txt"), row.names=F, col.names=F, quote=F)


# Venn Diagram ---
pdf(file.path(dir.plot, "venn_genelist_positive.pdf"))
p1 <- venn.diagram(x=list.genes, scaled = TRUE,
		filename = NULL, lwd=2.5, cex=1.5, fontfamily="Helvetica",
		height=2, width=2,  units="inches",
		col="black", fill=c("firebrick1","cornsilk3"),
		main="", main.cex=1.5, main.fontfamily="Helvetica",
		cat.fontfamily="Helvetica", cat.cex=1, cat.default.pos="text",
		cat.pos=c(0, 0),
		)
grid.draw(p1)
dev.off()


