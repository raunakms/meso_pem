### SET PATHS ----------------------------------------------------------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/03_enrichment_analysis")
dir.data <- file.path(dir.wrk, "data")
dir.plot <- file.path(dir.wrk, "plot")
dir.enrichment <- file.path(dir.wrk, "enrichment")

### LOAD LIBRARIES -----------------------------------------------------------------------------------------
library("gplots")
library("RColorBrewer")
library("stringr")


### DEFINE FILE --------------------------------------------------------------------------------------------
file.name <- "enrichment_h.all.v5.0.symbols.txt"

### DEFINE ENVIRONMENT -------------------------------------------------------------------------------------
batch <- "MESO_OUTLIERS_p01"
dir.batch <- file.path(dir.enrichment, batch)
dir.patients <- list.dirs(path=dir.batch, full.names=T, recursive=F)
ids.patients <- list.dirs(path=dir.batch, full.names=F, recursive=F)

### PROCESS EACH PATIENT DATA ------------------------------------------------------------------------------
list.dat <- list()
for(i in 1:length(ids.patients)){
	file.dat <- file.path(dir.patients[i], file.name)
	dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
	dat$ngenes <- unlist(lapply(str_split(dat$overlap.genes, ","), function(x) length(x)))
	dat <- subset(dat, dat$ngenes >= 2)
	
	if(nrow(dat) != 0){
		list.dat[[i]] <- data.frame(SampleID=ids.patients[i],
									Category=dat$Category,
									pvalue=dat$pvalue,
									fdr=dat$fdr,
									EnrichmentScore=-log10(dat$fdr),
									genes=dat$overlap.genes)
	} else{
		list.dat[[i]] <- data.frame(NULL)
	}
}
df <- do.call(rbind.data.frame, list.dat)
df$SampleID <- as.character(df$SampleID)
df$Category <- as.character(df$Category)
df$genes <- as.character(df$genes)


### CONSTRUCT PATHWAY MATRIX ----------------------------------------------------------------------------------
pathways <- unique(df$Category)
patients <- unique(df$SampleID)
mat <- matrix(0, nrow=length(pathways), ncol=length(patients), dimnames=list(pathways, patients))

for(i in 1:nrow(df)){
	x <- df$Category[i]
	y <- df$SampleID[i]
	mat[x,y] <- df$EnrichmentScore[i]
}


### PLOT PATHWAYS HEATMAP ----------------------------------------------------------------------------------
jColFun <- colorRampPalette(brewer.pal(n = 9, "Blues"))

file.plot <- file.path(dir.plot, "MESO_OUTLIERS_p01_pathway_heatmap_HALLMARKS.pdf")
pdf(file.plot, height=10, width=8)
	heatmap.2(mat, 
          col = jColFun(256),
          Colv=TRUE, Rowv=TRUE, 
          dendrogram ="both", trace="none",  scale="none",
          cexCol=1, cexRow=0.5, margins = c(5,20), 
          hclustfun = function(x) hclust(x, method = "ward.D2"), 
          distfun = function(x) dist(x, method = "euclidean"),
          colsep=c(1:200), rowsep=c(1:200),
          sepcolor="black", sepwidth=c(0.0005,0.0005),
          key="FALSE", keysize=0.8, density.info="none", symkey=0)
dev.off()



