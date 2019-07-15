### DEFINE LIBRARIES ---
library("stringr")
library("gplots")
library("RColorBrewer")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/08_compare_cnv")
dir.data <- file.path(dir.wrk, "data")
dir.output <- file.path(dir.wrk, "output")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")

### DEFINE FILES ---
file.dat <- file.path(dir.output, "cnv_meso_pem_pm_diff_summary_table.tsv")
file.expr <- file.path(dir.data, "cnv_meso_pem_pm_seg_mean_combined.tsv.gz")
file.script <- file.path("/Data/Raunak/softwares/bdvtools/array_process/array_preprocess.R")

### LOAD DATA ---
dat <- read.delim(file.dat, header=TRUE, stringsAsFactors=FALSE)

### LOAD EXPRESSION DATA ---
expr <- read.delim(file.expr, header=TRUE, stringsAsFactors=FALSE, row.names=1)
colnames(expr) <- str_replace_all(colnames(expr), "[.]", "-")

expr <- subset(expr, rownames(expr) %in% dat$Gene)
expr <- expr[match(dat$Gene, rownames(expr)),]

expr[is.na(expr)] <- 0
expr[expr < -1] <- -1


### ADD CHROMOSOME COLOR ---
dat$Color <- ""
dat$Color[which(dat$Chromosome == "1")] <- "#999999" #gray
dat$Color[which(dat$Chromosome == "11")] <- "#f781bf" #orange
dat$Color[which(dat$Chromosome == "12")] <- "#a65628" #brown
dat$Color[which(dat$Chromosome == "17")] <- "#984ea3" #purple
dat$Color[which(dat$Chromosome == "22")] <- "#4daf4a" #green
dat$Color[which(dat$Chromosome == "X")] <- "#377eb8" #blue

### DEFINE GROUPS ---
ids.peritoneal <- colnames(expr)[which(str_detect(colnames(expr), "MESO") == TRUE)]
ids.pleural <- colnames(expr)[which(str_detect(colnames(expr), "TCGA") == TRUE)]

### Visulize HEAT MAP ----
jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))
label.samples <- c(rep("#7A378B", length(ids.peritoneal)), rep("#ffc100", length(ids.pleural)))
label.genes <- dat$Color

file.plot <- file.path(dir.plot, "heatmap_cnv_meso_pem_pm_diff.pdf")
pdf(file.plot, height=4, width=7)
	heatmap.2(as.matrix(expr), 
          col = rev(jColFun(1024)),
          Colv=TRUE, Rowv=TRUE, na.rm=TRUE, revc=TRUE,
          dendrogram ="both", trace="none",  scale="none",
          cexCol=0.01, cexRow=0.01, symbreaks=TRUE, margin=c(5,5),
          hclustfun = function(x) hclust(x, method = "ward.D2"), 
          distfun = function(x) dist(x, method = "euclidean"),
          #colsep=c(1:5000), rowsep=c(1:5000),
          #sepcolor="black", sepwidth=c(0.0005,0.0005), 
		ColSideColors=label.samples,  RowSideColors=label.genes,
          xlab="", ylab="",
          key="TRUE", keysize=1, density.info="none", symkey=0,
		key.title = NA, key.xlab = NA, key.ylab = NA,
		key.par=list(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0)))
	#legend("topright", legend=c("Peritoneal","Pleural"), 
	#	  fill=c("#7A378B","#ffc100"), border=TRUE, bty="n", 
	#	  x.intersp = 1, y.intersp = 1, cex=1)
dev.off()

### CONVERT TO GCT FORMAT ----
source(file.script)
get.gct(dat=expr, file.gct=file.path(dir.output, "cnv_meso_pem_pm_seg_mean_209genes.gct"))

