### DEFINE PATH ---------------------------------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data")
dir.cna <- file.path(dir.wrk, "cnv/seq_call_refseq_genes_meso")
dir.exp <- file.path(dir.wrk, "expression/analysis")
dir.plot <- file.path(dir.wrk, "cnv/plot")
dir.analysis <- file.path(dir.wrk, "cnv/analysis")

### DEFINE FILES --------------------------------------------------------------------
file.cna <- file.path(dir.cna, "meso_cnv_seg_values_calls_parsed.tsv")
file.exp <- file.path(dir.exp, "meso_peritoneal_gene_expression_proteincoding.tsv.gz")

### Load Libraries ------------------------------------------------------------------
library("stringr")
library("ggplot2")
library("gridExtra")

### LOAD CNV DATA --------------------------------------------------------------------
dat.cna <- read.delim(file.cna, header=T, stringsAsFactors=F, row.names=1)
colnames(dat.cna) <- str_replace_all(colnames(dat.cna), "[.]", "-")

### Remove Genes with CNA = NA  ------
y <- apply(dat.cna, 1, function(x) length(which(is.na(x))))
del.index <- which(y > 8)
dat.cna <- dat.cna[-del.index,]

### LOAD EXPRESSION  DATA ------------------------------------------------------------
dat.exp <- read.delim(file.exp, header=T, stringsAsFactors=F, row.names=1)
colnames(dat.exp) <- str_replace_all(colnames(dat.exp), "[.]", "-")

### Remove Genes with Expr = 0 or NA  ------
y1 <- apply(dat.exp, 1, function(x) length(which(x == 0)))
y2 <- apply(dat.exp, 1, function(x) length(which(is.na(x))))
del.index1 <- which(y1 > 8)
del.index2 <- which(y2 > 8)
del.index <- unique(c(del.index1, del.index2))
dat.exp <- dat.exp[-del.index,]

### MATCH SAMPLES & GENES IN CNA & EXP ------------------------------------------------
dat.cna <- subset(dat.cna, select=colnames(dat.exp))

### GET COMMON GENES ----
genes <- intersect(rownames(dat.cna), rownames(dat.exp))
dat.cna <- subset(dat.cna, rownames(dat.cna) %in% genes)
dat.exp <- subset(dat.exp, rownames(dat.exp) %in% genes)

dat.cna <- dat.cna[match(genes, rownames(dat.cna)),]
dat.exp <- dat.exp[match(genes, rownames(dat.exp)),]

### CORRELATION OF TWO MATRICES ------------------------------------------------------
dat.cor <- data.frame(Gene=rownames(dat.cna))
dat.cor$Gene <- as.character(dat.cor$Gene)
dat.cor$R <- 0

for(i in 1:nrow(dat.cna)){
	x <- as.numeric(dat.cna[i,])
	y <- as.numeric(dat.exp[i,])

	x_na <- length(which(is.na(x)))
	y_na <- length(which(is.na(y)))

	if((x_na > 8) | (x_na > 8)) next
	dat.cor$R[i] <- cor(x, y, method="pearson", use="pairwise.complete.obs")
	#cat("PROCESSED", i, "OF", nrow(dat.cna), "\n", sep="\t")
}

dat.cor <- na.omit(dat.cor)
dat.cor.pass <- subset(dat.cor, dat.cor$R >= 0.75)
dat.cor.pass <- dat.cor.pass[order(dat.cor.pass$R, decreasing=T),]

dat.cor <- na.omit(dat.cor)
dat.cor <- dat.cor[order(dat.cor$R, decreasing=T),]

write.table(dat.cor, file.path(dir.cna, "cnv_expr_correlation.tsv"), sep="\t", row.names=F, col.names=T, quote=F)

### PLOT CORRELATION OF CNA & EXPR ------------------------------------------------------
list.plot <- list()
for(i in 1:nrow(dat.cor.pass)){
	gene <- dat.cor.pass$Gene[i]
	r <- dat.cor.pass$R[i]

	df <- data.frame(Gene=gene, 
					SampleID=unlist(lapply(str_split(colnames(dat.cna), "-"), function(x) x[2])), 
					CNA.value=as.numeric(dat.cna[gene,]),
					Exp.value=as.numeric(dat.exp[gene,]))

	list.plot[[i]] <- ggplot(df, aes(x=CNA.value, y=Exp.value)) + 
						geom_point(size=2, color="grey30", alpha=0.9) +
						stat_smooth(method="lm", geom = "smooth", formula = y ~ x, position = "identity", fullrange = FALSE, se=TRUE) +
						#coord_cartesian(xlim=c(-5, 2), ylim=c(-5, 20)) +
						geom_text(aes(label=SampleID, alpha=0.3, hjust=0, vjust=0, lineheight=1.5), color="black", size=2.5) +
						annotate("text", label = paste("R = ", round(r, 2), sep=""), x = 0, y = 13, size = 3, colour = "red") +
						facet_wrap(~ Gene, ncol = 1) +
						theme(
							axis.text = element_text(size = 8, color="black"),
							axis.title = element_text(size = 8, color="black"),
							plot.title = element_text(size = 8, color="black"),
							legend.position="none") +
						xlab("Copy Number Segment Mean") + 
						ylab("Gene Expression") + 
						ggtitle("") 	

	cat("PROCESSED: ", i, "OF", nrow(dat.cor.pass), "\n", sep="\t")					

}


#### PLOT -----------------------------------------------------------------------------------
file.plot <- file.path(dir.plot, "cnv_expr_correlation.pdf")
pdf(file.plot, height=4, width=4)
	for(i in 1:length(list.plot)){
		grid.arrange(list.plot[[i]] , ncol=1, nrow=1)
		cat("PROCESSED: ", i, "OF", nrow(dat.cor.pass), "\n", sep="\t")
	}
dev.off()



####
dat <- subset(dat.cor, dat.cor$R < -0.5)
write.table(dat, file.path(dir.cna, "cnv_expr_correlation_negative.tsv"), sep="\t", row.names=F, col.names=T, quote=F)


### PLOT CORRELATION OF CNA & EXPR ------------------------------------------------------
list.plot <- list()
for(i in 1:nrow(dat)){
	gene <- dat$Gene[i]
	r <- dat$R[i]

	df <- data.frame(Gene=gene, 
					SampleID=unlist(lapply(str_split(colnames(dat.cna), "-"), function(x) x[2])), 
					CNA.value=as.numeric(dat.cna[gene,]),
					Exp.value=as.numeric(dat.exp[gene,]))

	list.plot[[i]] <- ggplot(df, aes(x=CNA.value, y=Exp.value)) + 
						geom_point(size=2, color="grey30", alpha=0.9) +
						stat_smooth(method="lm", geom = "smooth", formula = y ~ x, position = "identity", fullrange = FALSE, se=TRUE) +
						#coord_cartesian(xlim=c(-5, 2), ylim=c(-5, 20)) +
						#geom_text(aes(label=SampleID, alpha=0.3, hjust=0, vjust=0, lineheight=1.5), color="black", size=2.5) +
						#annotate("text", label = paste("R = ", round(r, 2), sep=""), x = 0, y = 13, size = 3, colour = "red") +
						facet_wrap(~ Gene, ncol = 1) +
						theme(
							axis.text = element_text(size = 8, color="black"),
							axis.title = element_text(size = 8, color="black"),
							plot.title = element_text(size = 8, color="black"),
							legend.position="none") +
						xlab("Copy Number Segment Mean") + 
						ylab("Gene Expression") + 
						ggtitle("") 	

	cat("PROCESSED: ", i, "OF", nrow(dat), "\n", sep="\t")					

}

file.plot <- file.path(dir.plot, "cnv_expr_corr_negcor.pdf")
pdf(file.plot, height=4, width=4)
	for(i in 1:length(list.plot)){
		grid.arrange(list.plot[[i]] , ncol=1, nrow=1)
		cat("PROCESSED: ", i, "OF", nrow(dat), "\n", sep="\t")
	}
dev.off()


