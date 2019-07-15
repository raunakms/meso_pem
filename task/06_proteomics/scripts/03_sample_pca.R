#### LOAD LIBRARIES -------------------
library("stringr")
library("RColorBrewer")

#### DEFINE PATH ----------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.proteome <- file.path(dir.wrk, "data/proteome/processed_data")
dir.task <- file.path(dir.wrk, "task/06_proteomics")
dir.data <- file.path(dir.task, "data")
dir.output <- file.path(dir.task, "output")
dir.plot <- file.path(dir.task, "plot")

#### DEFINE FILE ---------------------
file.expr <- file.path(dir.proteome, "meso_proteome_proteomediscover_all_log2_dqnorm.tsv.gz")
file.des <- file.path(dir.proteome, "design_table_proteome.tsv")

#### LOAD DESIGN TABLE ----------------
des <- read.delim(file.des, header=T, stringsAsFactors=F)
ids.normal <- des$SequencingID[which(des$SampleType == "Normal")]
ids.tumor <- des$SequencingID[which(des$SampleType == "Tumor")]
ids.celline <- des$SequencingID[which(des$SampleType == "Celline")]
ids <- c(ids.normal, ids.tumor, ids.celline)

#### LOAD PROTEIN EXPRESSION DATA -----
expr <- read.delim(file.expr, header=T, stringsAsFactors=F, row.names=1)
colnames(expr) <- str_replace_all(colnames(expr), "[.]", "-")
expr <- subset(expr, select=ids)
expr[is.na(expr)] <- 0

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
expr.subset <- subset(expr, rownames(expr) %in% genes.select)

dat.pca.new <- prcomp(expr.subset)
dat.eigenvectors <- as.data.frame(dat.pca.new $rotation[,1:2])

### Plot ---
s <- str_replace_all(rownames(dat.eigenvectors), "MESO-", "")
file.plot <- file.path(dir.plot, "PCA_variable_expression_2SD.pdf")
pdf(file.plot, height=4, width=4)
plot(x=dat.eigenvectors$PC1, y=dat.eigenvectors$PC2,
		xlim=c(-0.4,0.3), 
		ylim=c(-0.5,0.3),
		type="p", pch=1, col="black",
		cex=1, cex.axis=0.5, cex.lab=0.5, cex.main=0.5,
		las=1, tck = -.03, xlab="PC1", ylab="PC2",
		main="PCA - 2 StdDev")
text(x=dat.eigenvectors$PC1, y=dat.eigenvectors$PC2, labels=s, cex=0.5, pos=1, offset=0.5)
dev.off()

