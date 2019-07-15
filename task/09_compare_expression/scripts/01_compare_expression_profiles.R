### SET DIRECTORIES -------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task_comparision/01_expression_comparision")
dir.data <- file.path(dir.wrk, "data")
dir.plot <- file.path(dir.wrk, "plot")

### Load Libraries ----------------------------------------------------------
library("stringr")

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

### Principal Component Analysis ----------------------------------------------
dat.pca <- prcomp(t(expr))

dat.eigenvectors <- as.data.frame(dat.pca$rotation[,1:3])
dat.eigenvectors$max.pc <- apply(dat.eigenvectors, 1, max)

cutoff <- sd(dat.eigenvectors$max.pc) * 2
dat.eigenvectors$var.select <- ifelse(dat.eigenvectors$max.pc > cutoff, 1, 0)
genes.select <- rownames(dat.eigenvectors)[which(dat.eigenvectors$var.select == 1)]

#### Select expression matrix with subset gene features ------------------------
expr.var <- subset(expr, rownames(expr) %in% genes.select)

### Visulize Batch Effects -----------------------------------------------------
file.viz <- file.path("/Data/Raunak/softwares/bdvtools/array_process/visualizeBatch.R")
dir.output <- file.path(dir.plot, "viz_batch_1")

source(file.viz)
visualizeBatch(expr=as.matrix(expr.var), outcome=c(colnames(expr1)[-1], des$Subtype), batch=c(rep(1, 15), rep(2, 22), rep(3, 5), rep(4, 57), rep(5, 2)), outDir=dir.output)

#batch=c(rep(1, 15), rep(2, 87))

### T-Test between Peritoneal and Pleural ------------------------------------------
source(file.stat)
ids.peritoneal <- colnames(expr1)[-1]
ids.pleural <- des$SampleCode[which(des$Subtype == "Epithelioid")]
df <- get.diffexpr(expr, class1=ids.peritoneal, class2=ids.pleural)


### Visulize Batch Effects [Peritoneal and TCGA-Eptheleiod only]----------------------
expr.var <- subset(expr.var, select=c(colnames(expr1)[-1], des$SampleCode[which(des$Subtype == "Epithelioid")]))
file.viz <- file.path("/Data/Raunak/softwares/bdvtools/array_process/visualizeBatch.R")
dir.output <- file.path(dir.plot, "viz_batch_2")

source(file.viz)
visualizeBatch(expr=as.matrix(expr.var), outcome=c(colnames(expr1)[-1], rep("Epithelioid", 57)), batch=c(rep(1, 15), rep(4, 57)), outDir=dir.output)
