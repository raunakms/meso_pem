### LOAD LIBRAIES ---
library("stringr")
library("gplots")
library("RColorBrewer")

### DEFINE PATH ---
dir.wrk <- file.path("/collinsgroup/Raunak/projects/MESO_peritoneal/task/08_compare_cnv")
dir.data <- file.path(dir.wrk, "data")
dir.output <- file.path(dir.wrk, "output")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")
dir.pleural <- file.path("/collinsgroup/Raunak/HITnDRIVE/datasets/TCGA_MESO/data/cnv/seq_call_refseq_genes_meso")
dir.ref <- file.path("/home/Collins/databases/Chr_gff")

### DEFINE FILES ---
file.peritoneal <- file.path(dir.data, "meso_cnv_seg_values_calls_parsed.tsv")
file.pleural <- file.path(dir.pleural, "tcga_meso_cnv_seg_values_calls_parsed.tsv.gz")
#file.sampleid <- file.path(dir.ov, "sample_ids.txt")
file.script <- file.path("/collinsgroup/Raunak/softwares/bdvtools/array_process/wilcox_rank_test.R")
file.ref <- file.path(dir.ref, "Homo_sapiens.GRCh38.87_gene_proteincoding.bed.gz")

## Remove Genes with Expr = NA in > 25% of the samples ------
remove.na.matrix <- function(dat, cutoff){
	y <- apply(dat, 1, function(x) length(which(is.na(x))))
	del.index <- which(y > (ncol(dat) * cutoff))
	if(length(del.index) != 0){
		dat <- dat[-del.index,]
	} else{
		return(dat)
	}
	return(dat)
}

### LOAD REFERENCE GENOME ---
dat.ref <- read.delim(file.ref, header=F, stringsAsFactors=F)

### LOAD MESO-PERITONEAL CNV ----
dat.peritoneal <-  read.delim(file.peritoneal, header=T, stringsAsFactors=F, row.names=1)
colnames(dat.peritoneal) <- str_replace_all(colnames(dat.peritoneal), "[.]", "-")
dat.peritoneal <- remove.na.matrix(dat=dat.peritoneal, cutoff=0.25)

### LOAD MESO-PLEURAL CNV ----
dat.pleural <-  read.delim(file.pleural, header=T, stringsAsFactors=F, row.names=1)
colnames(dat.pleural) <- str_replace_all(colnames(dat.pleural), "[.]", "-")
dat.pleural <- remove.na.matrix(dat=dat.pleural, cutoff=0.25)


#### Subset data by common genes ---
genes <- intersect(intersect(rownames(dat.peritoneal), rownames(dat.pleural)), unique(dat.ref$V4))
dat.peritoneal <- subset(dat.peritoneal, rownames(dat.peritoneal) %in% genes)
dat.pleural <- subset(dat.pleural, rownames(dat.pleural) %in% genes)


### Combine Data ---
sampleids.peritoneal <- colnames(dat.peritoneal)
sampleids.pleural <- colnames(dat.pleural)

sample.class <- c(rep("PeM", length(sampleids.peritoneal)), rep("PM", length(sampleids.pleural)))
dat <- cbind(dat.peritoneal, dat.pleural)


#### WRITE OUTPUT -----------
file.output <- file.path(dir.data, "cnv_meso_pem_pm_seg_mean_combined.tsv")
write.table(dat, file.output, sep="\t", row.names=T, col.names=NA, quote=F)

file.output <- file.path(dir.data, "cnv_meso_pem_pm_combined_class.txt")
write.table(sample.class, file.output, sep="\t", row.names=F, col.names=F, quote=F)

#### RELOAD DATA ---
file.dat <- file.path(dir.data, "cnv_meso_pem_pm_seg_mean_combined.tsv.gz")
dat <- read.delim(file.dat, header=T, row.names=1)
colnames(dat) <- str_replace_all(colnames(dat), "[.]", "-")
#dat <- remove.na.matrix(dat=dat, cutoff=0.25)
dat[is.na(dat)] <- 0

sampleids.peritoneal <- colnames(dat)[1:19]
sampleids.pleural <- colnames(dat)[20:104]

write.table(rownames(dat), file.path(dir.output, "cnv_meso_pem_pm_background_genelist.txt"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


#### DIFFERENTIAL CNA ANALYSIS ----
source(file.script)
df <- get.wilcox.rank.test(dat, class1=sampleids.peritoneal, class2=sampleids.pleural)

df <- subset(df, df$fdr <= 5e-04)
genes.dg <- df$Gene

write.table(df, file.path(dir.output, "cnv_meso_pem_pm_diffexpr_results.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(df$Gene, file.path(dir.output, "cnv_meso_pem_pm_diffexpr_features.txt"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


### ADD ANNOTATIONS ----
file.annot <- file.path(dir.output, "cnv_meso_pem_pm_diffexpr_features.tsv")
annot <- read.delim(file.annot, header=T, stringsAsFactors=F)
annot <- annot[-which(str_detect(annot$Chromosome, "CHR_") == TRUE),]
annot$KaryotypeBand <- unlist(lapply(str_split(annot$KaryotypeBand, "[.]"), function(x) x[1]))
annot$Cytoband <- apply(annot, 1, function(x) paste(as.character(x[2]), as.character(x[3]), sep=":"))

### LOAD DIFF.EXPR SUMMARY ---
df <- read.delim(file.path(dir.output, "cnv_meso_pem_pm_diffexpr_results.tsv"),header=T, stringsAsFactors=F)
df <- subset(df, df$Gene %in% annot$Gene)
a <- annot
a <- a[match(df$Gene, a$Gene),]
df$ChromosomeLocation <- a$Cytoband
df <- subset(df, select=c("Gene","ChromosomeLocation","pvalue","fdr"))
write.table(df, file.path(dir.output, "cnv_meso_pem_pm_diffexpr_results_supplementarytable.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

### SUBSET CNV DATA ---
dat <- subset(dat, rownames(dat) %in% annot$Gene)
write.table(dat, file.path(dir.data, "cnv_meso_pem_pm_diffexpr_seg.tsv"), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

#### LOAD DATA ----
dat <- read.delim(file.path(dir.data, "cnv_meso_pem_pm_diffexpr_seg.tsv"), header=T, stringsAsFactors=FALSE)
colnames(dat) <- str_replace_all(colnames(dat), "[.]", "-")
colnames(dat)[1] <- "Chr"
dat <- dat[match(annot$Gene, dat$Chr),]
dat$Chr <- annot$Cytoband

### REMOVE X Chromosome (due to differences in Male and Female patient composition in Pm and PeM)--- 
dat <- dat[-which(str_detect(dat$Chr, "X") == TRUE),]

### GET UNIQUE FEATURES ---
file.array <- file.path("/collinsgroup/Raunak/softwares/bdvtools/array_process/array_preprocess.R")
source(file.array)
dat <- getunique.gene.expression(dat)


### Visulize HEAT MAP ----
jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))
label.groups <- c(rep("#7A378B", 19), rep("#ffc100", 85))

file.plot <- file.path(dir.plot, "heatmap_cnv_meso_pem_pm_diff_cytoband_scaled.pdf")
pdf(file.plot, height=5, width=6.8)
	heatmap.2(as.matrix(dat), 
          col = rev(jColFun(1024)),
          Colv=FALSE, Rowv=TRUE, 
          dendrogram ="row", trace="none",  scale="row",
          cexCol=0.5, cexRow=0.5, symbreaks=TRUE, margin=c(10,5),
          hclustfun = function(x) hclust(x, method = "ward.D2"), 
          distfun = function(x) dist(x, method = "euclidean"),
          colsep=c(1:500), rowsep=c(1:500),
          sepcolor="#d9d9d9", sepwidth=c(0.0005,0.0005), 
		  ColSideColors=label.groups,  
		  xlab="", ylab="",
          key="TRUE", keysize=1, density.info="none", symkey=0,
		  key.title = NA, key.xlab = NA, key.ylab = NA,
		  key.par=list(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0)))
dev.off()


df$Chromosome <- ""
df$Start <- ""
df$End <- ""
df$Band <- ""
for(i in 1:nrow(df)){
	index <- which(annot$Gene == df$Gene[i])
	if(length(index) == 0) next
	index <- index[1]
	df$Chromosome[i] <- annot$Chromosome[index]
	df$Start[i] <- annot$Start[index]
	df$End[i] <- annot$End[index]
	df$Band[i] <- annot$Band[index]
}

df <- subset(df, df$Chromosome != "")
df <- df[with(df, order(Chromosome, Band, Start)),]



#### WRITE OUTPUT -----------
file.output <- file.path(dir.output, "cnv_meso_pem_pm_diff_summary_table.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)

file.deg <- file.path(dir.output, "cnv_meso_pem_pm_diff_genelist.txt")
write.table(genes.dg , file.deg, sep="\t", row.names=F, col.names=F, quote=F)
