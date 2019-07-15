### Load Libraries -----
library("stringr")
library("ggplot2")
library("gridExtra")
library("gplots")
library("RColorBrewer")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.analysis <- file.path(dir.wrk, "task/07_compare_mRNA_protein")
dir.output <- file.path(dir.analysis, "output")
dir.plot <- file.path(dir.analysis, "plot")
dir.script <- file.path(dir.analysis, "scripts")
dir.prot <- file.path(dir.wrk, "data/proteome/processed_data")

### DEFINE FILES ---
file.des <- file.path(dir.wrk, "data/proteome/processed_data/design_table_proteome.tsv")
file.expr <- file.path(dir.prot, "meso_proteome_proteomediscover_all_log2_dqnorm.tsv.gz")

### LOAD DATA: Protein ---
expr <- read.delim(file.expr, header=TRUE, stringsAsFactors=FALSE, row.names=1)
colnames(expr) <- str_replace_all(colnames(expr), "[.]", "-")
expr[is.na(expr)] <- 0

### LOAD DESIGN TABLE: Protein ---
des <- read.delim(file.des, header=T, stringsAsFactors=F)
des <- subset(des, des$SampleType == "Tumor")
grp0 <- des$SequencingID[which(des$Group.3p21 == 0)]
grp1 <- des$SequencingID[which(des$Group.3p21 == 1)]

### ARRANGE EXPRESSION DATA ---
expr <- subset(expr, select=c(grp0, grp1))
colnames(expr) <- str_replace_all(colnames(expr), "T", "")

### GET GENES ---
tag.genes <- c("SETD2","PBRM1","SMARCC1")
all.genes <- row.names(expr)

### PREPARE MATRIX ---
mat <- matrix(0, nrow=length(all.genes), ncol=length(tag.genes), dimnames=list(all.genes, tag.genes))

### COMPUTE PEARSON CORRELATION ---
# FUNCTION --
getCorrelation <- function(gene1, gene2){
	x <- as.numeric(expr[gene1,])
	y <- as.numeric(expr[gene2,])

	correlation.test <- cor.test(x, y, method="pearson")
	r <- as.numeric(correlation.test$estimate)

	return(r)
}

# LOOP --
for(i in 1:nrow(mat)){
	mat[i,1] <- getCorrelation(gene1="SETD2", gene2=rownames(mat)[i])
	mat[i,2] <- getCorrelation(gene1="PBRM1", gene2=rownames(mat)[i])
	mat[i,3] <- getCorrelation(gene1="SMARCC1", gene2=rownames(mat)[i])

	cat("PROCESSED:", i, "OF", nrow(mat), "\n", sep="\t")
}
mat <- na.omit(mat)

### Visulize HEAT MAP ----
jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))

file.plot <- file.path(dir.plot, "heatmap_protein_correlation.pdf")
pdf(file.plot, height=5, width=20)
	heatmap.2(t(as.matrix(mat)), 
          col = rev(jColFun(1024)),
          Colv=TRUE, Rowv=TRUE, 
          dendrogram ="both", trace="none",  scale="none",
          cexCol=0.1, cexRow=1, symbreaks=TRUE, margin=c(5,5),
          hclustfun = function(x) hclust(x, method = "ward.D2"), 
          distfun = function(x) dist(x, method = "euclidean"),
          #colsep=c(1:500), rowsep=c(1:500),
          #sepcolor="white", sepwidth=c(0.0005,0.0005), 
		  #ColSideColors=label.groups,  
		  xlab="", ylab="",
          key="TRUE", keysize=1, density.info="none", symkey=0,
		  key.title = NA, key.xlab = NA, key.ylab = NA,
		  key.par=list(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0)))
dev.off()


#### SELECTED GENES FROM CLIQUE ---
genes.cq <- c("SMARCA1","SMARCA5","SMARCAD1","SMARCAL1","SMARCC1","SMARCD3","SMARCB1","SMARCA2","SMARCC2","ARID1A","SMARCE1","SMARCA4","PBRM1","SMARCD1","SETD2","SMARCD2","ACTL6A")
expr.cq <- subset(expr, rownames(expr) %in% genes.cq)

### PREPARE MATRIX ---
mat <- matrix(0, nrow=nrow(expr.cq), ncol=nrow(expr.cq), dimnames=list(rownames(expr.cq), rownames(expr.cq)))


# LOOP --
for(i in 1:nrow(mat)){
	for(j in 1:ncol(mat)){
		mat[i,j] <- getCorrelation(gene1=rownames(mat)[i], gene2=colnames(mat)[j])
	}
	cat("PROCESSED:", i, "OF", nrow(mat), "\n", sep="\t")
}



### Visulize HEAT MAP ----
jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))

file.plot <- file.path(dir.plot, "heatmap_protein_correlation_clique_SMARCC1.pdf")
pdf(file.plot, height=4, width=4)
	heatmap.2(as.matrix(mat), 
          col = rev(jColFun(1024)),
          Colv=TRUE, Rowv=TRUE, 
          dendrogram ="both", trace="none",  scale="none",
          cexCol=0.8, cexRow=0.8, symbreaks=TRUE, margin=c(8,8),
          hclustfun = function(x) hclust(x, method = "ward.D2"), 
          distfun = function(x) dist(x, method = "euclidean"),
          colsep=c(1:500), rowsep=c(1:500),
          sepcolor="white", sepwidth=c(0.0005,0.0005), 
		  #ColSideColors=label.groups,  
		  xlab="", ylab="",
          key="TRUE", keysize=1, density.info="none", symkey=0,
		  key.title = NA, key.xlab = NA, key.ylab = NA,
		  key.par=list(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0)))
dev.off()
