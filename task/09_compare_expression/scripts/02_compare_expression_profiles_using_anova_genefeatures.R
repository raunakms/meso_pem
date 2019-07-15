### SET DIRECTORIES -------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task_comparision/01_expression_comparision")
dir.data <- file.path(dir.wrk, "data")
dir.plot <- file.path(dir.wrk, "plot")

### Load Libraries ----------------------------------------------------------
library("stringr")
library("gplots")
library("RColorBrewer")

### Define Files ------------------------------------------------------------
file.expr1 <- file.path(dir.data, "meso_peritoneal_gene_expression_znorm_proteincoding.tsv")
file.expr2 <- file.path(dir.data, "TCGA-MESO_normalized_expr_znorm.tsv")
file.des <- file.path(dir.data, "design_table.tsv")
file.stat <- file.path("/Data/Raunak/softwares/bdvtools/array_process/diff_expr.R")

### Load Design Table - TCGA --------------------------------------------------
des <- read.delim(file.des, header=T, stringsAsFactors=F)
des <- subset(des, des$SampleType == "Tumor")
des <- des[-which(des$Subtype == ""),]
des <- des[order(des$Subtype, decreasing=F),]

### Load Expression Data -----------------------------------------------------
expr1 <- read.delim(file.expr1, header=T, stringsAsFactors=F)
colnames(expr1) <- str_replace(colnames(expr1), "[.]", "-")

expr2 <- read.delim(file.expr2, header=T, stringsAsFactors=F)
expr2 <- subset(expr2, select=c("X",des$SampleCode))
expr2 <- expr2[,match(c("X",des$SampleCode), colnames(expr2))]

### Merge Expression Data -----------------------------------------------------
expr <- merge(expr1, expr2, by="X")
rownames(expr) <- expr$X
expr <- as.matrix(expr[,-1])
del.index <- as.numeric(which(is.na(rowMeans(expr))))
expr <- expr[-del.index,]

### Remove Genes with Expr = 0 in >4 samples  ---------------------------------
y <- apply(expr, 1, function(x) length(which(x == 0)))
del.index <- which(y > 4)
expr <- expr[-del.index,]

#### Select ANOVA significant genes --------------------------------------------
file.dm <- file.path("/Data/Raunak/HITnDRIVE/datasets/TCGA_MESO/analysis/06_subtype_diffexpr/data", "TCGA-MESO_subtype_anova.tsv")
dat.summary <- read.delim(file.dm, header=T, stringsAsFactors=F)
dat.summary <- subset(dat.summary, dat.summary$anovaFDR <= 0.1)
genes.select <- dat.summary$Gene

write.table(genes.select, file.path(dir.data, "Anova_significant_genelist_for_tcga_meso_subtypes.txt"), row.names=F, col.names=F, quote=F)

#### Select expression matrix with subset gene features ------------------------
expr.var <- subset(expr, rownames(expr) %in% genes.select)

### Visulize Batch Effects -----------------------------------------------------
file.viz <- file.path("/Data/Raunak/softwares/bdvtools/array_process/visualizeBatch.R")
dir.output <- file.path(dir.plot, "viz_batch_anova")

source(file.viz)
visualizeBatch(expr=as.matrix(expr.var), outcome=c(colnames(expr1)[-1], des$Subtype), batch=c(rep(1, 15), rep(2, 22), rep(3, 5), rep(4, 57), rep(5, 2)), outDir=dir.output)

### Visulize HEAT MAP ----------------------------------------------------------
jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))
label.groups <- c(rep("#ffff99", 15), rep("#C44D58", 22), rep("#4ECDC4", 5), rep("#556270", 57), rep("#C7F464", 2))

file.plot <- file.path(dir.plot, "heatmap_expr_anova_genes.pdf")
pdf(file.plot, height=8, width=9)
	heatmap.2(expr.var, 
          col = rev(jColFun(1024)),
          Colv=TRUE, Rowv=TRUE, 
          dendrogram ="both", trace="none",  scale="none",
          cexCol=0.5, cexRow=0.01, symbreaks=TRUE, margin=c(6,20),
          hclustfun = function(x) hclust(x, method = "ward.D2"), 
          distfun = function(x) dist(x, method = "euclidean"),
          #colsep=c(1:500), rowsep=c(1:500),
          #sepcolor="white", sepwidth=c(0.0005,0.0005), 
		  ColSideColors=label.groups,  xlab="Patients", ylab="Gene",
          key="TRUE", keysize=1, density.info="none", symkey=0,
		  key.title = NA, key.xlab = NA, key.ylab = NA,
		  key.par=list(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0)))
		  legend("topright", legend=c("PeritonealMESO","Biphasic","DiffuseMaglignant","Epithelioid","Sarcomatoid"), fill=c("#ffff99","#C44D58", "#4ECDC4", "#556270", "#C7F464"), border=TRUE, bty="n", x.intersp = 1, y.intersp = 1, cex=1)
dev.off()


