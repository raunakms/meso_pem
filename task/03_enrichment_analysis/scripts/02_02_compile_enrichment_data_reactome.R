### LOAD LIBRARIES -----------------------------------------------------------------------------------------
library("gplots")
library("RColorBrewer")
library("stringr")

### SET PATHS ----------------------------------------------------------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/03_enrichment_analysis")
dir.data <- file.path(dir.wrk, "data")
dir.plot <- file.path(dir.wrk, "plot")
dir.enrichment <- file.path(dir.wrk, "enrichment")


### DEFINE FILE --------------------------------------------------------------------------------------------
file.name <- "enrichment_c2.cp.reactome.v6.0.symbols.txt"

### DEFINE ENVIRONMENT -------------------------------------------------------------------------------------
batch <- "MESO_OUTLIERS_mRNA_p01"
dir.batch <- file.path(dir.enrichment, batch)
dir.patients <- list.dirs(path=dir.batch, full.names=T, recursive=F)
ids.patients <- list.dirs(path=dir.batch, full.names=F, recursive=F)

### PROCESS EACH PATIENT DATA ------------------------------------------------------------------------------
list.dat <- list()
for(i in 1:length(ids.patients)){
	file.dat <- file.path(dir.patients[i], file.name)

	if(!file.exists(file.dat)) next

	dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
	#dat$Category <- unlist(lapply(str_split(dat$Category, "_"), function(x) paste(x[2:length(x)], collapse=" ")))
	dat$ngenes <- unlist(lapply(str_split(dat$overlap.genes, ","), function(x) length(x)))

	dat <- subset(dat, dat$ngenes >= 4)
	
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
mat <- matrix(0, nrow=length(pathways), ncol=length(ids.patients), dimnames=list(pathways, ids.patients))

for(i in 1:nrow(df)){
	x <- df$Category[i]
	y <- df$SampleID[i]
	mat[x,y] <- df$EnrichmentScore[i]
}

write.table(mat, file.path(dir.data, "MESO_OUTLIERS_mRNA_p01_pathway_enrichment_reactome.tsv"), sep="\t", row.names=T, col.names=NA, quote=F)


### RELOAD DATA ----
file.mat <- file.path(dir.data, "MESO_OUTLIERS_mRNA_p01_pathway_enrichment_reactome.tsv")
mat <- read.delim(file.mat, header=T, stringsAsFactors=FALSE, row.names=1)
colnames(mat) <- str_replace_all(colnames(mat), "[.]", "-")

### PLOT PATHWAYS HEATMAP ----------------------------------------------------------------------------------
jColFun <- colorRampPalette(brewer.pal(n = 9, "Greys"))

file.plot <- file.path(dir.plot, "MESO_OUTLIERS_p01_pathway_heatmap_REACTOME.pdf")
pdf(file.plot, height=6, width=6)
	heatmap.2(as.matrix(mat), 
          col = jColFun(256),
          Colv=TRUE, Rowv=TRUE, 
          dendrogram ="both", trace="none",  scale="none",
          cexCol=0.5, cexRow=0.3, margins = c(5,20), 
          hclustfun = function(x) hclust(x, method = "ward.D2"), 
          distfun = function(x) dist(x, method = "euclidean"),
          colsep=c(1:150), rowsep=c(1:150),
          sepcolor="white", sepwidth=c(0.0005,0.0005),
          key="TRUE", keysize=0.8, density.info="none")
dev.off()



