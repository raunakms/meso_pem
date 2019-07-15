#### LOAD LIBRARIES ------
library("stringr")
library("reshape2")
library("gplots")
library("RColorBrewer")

#### DEFINE PATH ------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/05_consensus_clustering")
dir.data <- file.path(dir.wrk, "data")
dir.plot <- file.path(dir.wrk, "plot")
dir.output <- file.path(dir.wrk, "output")
dir.des3p21 <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/annotation")

#### DEFINE FILE -----
file.dat <- file.path(dir.data, "meso_cnv_seg_values_calls_parsed.tsv")
file.des <- file.path(dir.plot, "02_MESO_consensus_cluster_cnv/02_MESO_consensus_cluster_cnv.tsv")
file.des3p21 <- file.path(dir.des3p21, "design_table_3p21genes.tsv")
file.diffexpr <- file.path("/Data/Raunak/softwares/bdvtools/array_process/diff_expr.R")

#### LOAD DESIGN TABLE ---
des <- read.delim(file.des, header=T, stringsAsFactors=FALSE)
des <- subset(des, des$k == 4)
id_class1 <- des$SampleID[which(des$Class == "Class_1")]
id_class2 <- des$SampleID[which(des$Class == "Class_2")]
id_class3 <- des$SampleID[which(des$Class == "Class_3")]
id_class4 <- des$SampleID[which(des$Class == "Class_4")]
sampleids <- c("MESO-17","MESO-02","MESO-01","MESO-03","MESO-12","MESO-04","MESO-15","MESO-09","MESO-05","MESO-06","MESO-19","MESO-14","MESO-10","MESO-07","MESO-08","MESO-11","MESO-13","MESO-18A","MESO-18E")


#### LOAD DATA ----
dat <- read.delim(file.dat, header=T, row.names=1)
colnames(dat) <- str_replace(colnames(dat), "[.]", "-")


#### Load Expression Data ---------------------------------------------
dat <- read.delim(file.dat, header=T, row.names=1)
colnames(dat) <- str_replace(colnames(dat), "[.]", "-")

### Remove Genes with CNA = NA  ---------------------------------
y <- apply(dat, 1, function(x) length(which(is.na(x))))
del.index <- which(y >= 1)
dat <- dat[-del.index,]

### Principal Component Analysis --------------------------------------
#expr <- na.omit(expr)
#dat.pca <- prcomp(t(dat))

#dat.eigenvectors <- as.data.frame(dat.pca$rotation[,1:3])
#dat.eigenvectors$max.pc <- apply(dat.eigenvectors, 1, max)

#cutoff <- sd(dat.eigenvectors$max.pc) * 2
#dat.eigenvectors$var.select <- ifelse(dat.eigenvectors$max.pc > cutoff, 1, 0)
#genes.select <- rownames(dat.eigenvectors)[which(dat.eigenvectors$var.select == 1)]

#### Select expression matrix with subset gene features ----------------
#dat.var <- subset(dat, rownames(dat) %in% genes.select)
dat.var <- dat

#### MELT DATA ---
df <- melt(as.matrix(t(dat.var)))
colnames(df) <- c("SampleID","Feature","value")
df$SampleID <- as.character(df$SampleID)
df$Feature <- as.character(df$Feature)

#### ADD GROUP INFO ---
df$Group <- ""
df$Group[which(df$SampleID %in% id_class1)] <- "Class_1"
df$Group[which(df$SampleID %in% id_class2)] <- "Class_2"
df$Group[which(df$SampleID %in% id_class3)] <- "Class_3"
df$Group[which(df$SampleID %in% id_class4)] <- "Class_4"

df$Group <- factor(df$Group, levels=c("Class_1","Class_2","Class_3","Class_4"))

### CALL FUNCTION: KruskalWallisTest ----
source(file.diffexpr)
df_kw <- calcKruskalWallisTest_parallel(df)
df_kw_sig <- subset(df_kw, df_kw$pvalue < 0.01)
df_kw_sig$Feature <- as.character(df_kw_sig$Feature)

write.table(df_kw_sig, file.path(dir.output, "meso_cnv_all_KruskalWallisTest_results.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(df_kw_sig$Feature, file.path(dir.output, "meso_cnv_all_KruskalWallisTest_features.txt"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#> dim(df_kw_sig)
#[1] 3383    4


### ADD ANNOTATIONS ----
file.annot <- file.path(dir.output, "meso_cnv_all_KruskalWallisTest_features.tsv")
annot <- read.delim(file.annot, header=T, stringsAsFactors=F)
annot <- annot[-which(str_detect(annot$Chromosome, "CHR_") == TRUE),]
annot$KaryotypeBand <- unlist(lapply(str_split(annot$KaryotypeBand, "[.]"), function(x) x[1]))
annot$Cytoband <- apply(annot, 1, function(x) paste(as.character(x[2]), as.character(x[3]), sep=":"))

# REPORT THIS VALUE ---
#> dim(annot)
#[1] 3028   4


### LOAD KruskalWallisTest SUMMARY ---
df <- read.delim(file.path(dir.output, "meso_cnv_all_KruskalWallisTest_results.tsv"), header=T, stringsAsFactors=F)
df <- subset(df, df$Feature %in% annot$Gene)
a <- annot
a <- a[match(df$Feature, a$Gene),]
df$ChromosomeLocation <- a$Cytoband
df <- subset(df, select=c("Feature","ChromosomeLocation","pvalue","fdr"))
write.table(df, file.path(dir.output, "meso_cnv_all_KruskalWallisTest_results_supplementarytable.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

### SUBSET CNV DATA ---
dat <- subset(dat, rownames(dat) %in% df$Feature)
write.table(dat, file.path(dir.data, "meso_cnv_all_KruskalWallisTest_cnv.tsv"), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)


#### LOAD DATA ----
dat <- read.delim(file.path(dir.data, "meso_cnv_all_KruskalWallisTest_cnv.tsv"), header=T, stringsAsFactors=FALSE)
colnames(dat) <- str_replace_all(colnames(dat), "[.]", "-")
colnames(dat)[1] <- "Chr"
dat <- dat[match(annot$Gene, dat$Chr),]
dat$Chr <- annot$Cytoband


### GET UNIQUE FEATURES ---
file.array <- file.path("/Data/Raunak/softwares/bdvtools/array_process/array_preprocess.R")
source(file.array)
dat <- getunique.gene.expression(dat)
dat <- dat[,match(sampleids, colnames(dat))]

### GET GCT FORMAT ---
#get.gct(dat, file.gct=file.path(dir.data, "meso_cnv_KruskalWallisTest_cnv.gct"))

### LOAD DESIGN TABLE ---
des3p21 <- read.delim(file.des3p21, header=TRUE, stringsAsFactors=FALSE)
grp0 <- des3p21$Sample.ID[which(des3p21$Group.3p21 == 0)]
grp1 <- des3p21$Sample.ID[which(des3p21$Group.3p21 == 1)]

grp.col <- rep(NA, length(sampleids))
for(i in 1:length(sampleids)){
	grp <- des3p21$Group.3p21[which(des3p21$Sample.ID %in% sampleids[i])]

	if(grp == 1){
		grp.col[i] <- "blue"
	} else{
		grp.col[i] <- "#fffc00"
	}
}

### Visulize HEAT MAP ----
jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))
label.groups <- grp.col

file.plot <- file.path(dir.plot, "heatmap_cnv_all_KruskalWallisTest_PeM.pdf")
pdf(file.plot, height=6, width=5)
	heatmap.2(as.matrix(dat), 
          col = rev(jColFun(1024)),
          Colv=FALSE, Rowv=TRUE, 
          dendrogram ="row", trace="none",  scale="row",
          cexCol=1, cexRow=0.5, symbreaks=TRUE, margin=c(10,10),
          hclustfun = function(x) hclust(x, method = "ward.D2"), 
          distfun = function(x) dist(x, method = "euclidean"),
          colsep=c(1:1500), rowsep=c(1:1500),
          sepcolor="white", sepwidth=c(0.0005,0.0005), 
		  ColSideColors=label.groups,  
		  xlab="", ylab="",
          key="TRUE", keysize=1, density.info="none", symkey=0,
		  key.title = NA, key.xlab = NA, key.ylab = NA,
		  key.par=list(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0)))
dev.off()


