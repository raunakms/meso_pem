### FUNCTON TO COMPUTE PATHWAY Z-SCORE --------
file.array <- file.path("/Data/Raunak/softwares/bdvtools/array_process/array_preprocess.R")
source(file.array)

#### CALLFUNCTIONS -----------
get.pathway.heatmap <- function(expr, file.pathway, database, group.names, file.plot, plot.height, plot.width){
	pathz <- get.pathway.zcore(expr=expr, file.pathway=file.pathway, database=database)
	plot.heatmap.pathway.zscore(dat=pathz, database=database, group.names=group.names, file.plot=file.plot, plot.height=plot.height, plot.width=plot.width)
}

#### GET PATHWAY ZSCORE ----------------------
get.pathway.zcore <- function(expr, file.pathway, database){
	# LOAD LIBRARIES -----
	require("stringr")

	# LOAD Pathway ------
	dat.path <- read.delim(file.pathway, header=T, stringsAsFactors=F)
	dat.path <- subset(dat.path, dat.path$fdr <= 0.01)
	dat.path$ngenes <- unlist(lapply(str_split(dat.path$overlap.genes, ","), function(x) length(x)))
	dat.path <- subset(dat.path, dat.path$ngenes >= 3)

	# FILTER BY DATABASE ----
	db <- c("KEGG","REACTOME","HALLMARK","GOBP")
	#if(database %in% db){
	#	dat.path$Category <- unlist(lapply(lapply(str_split(dat.path$Category, "_"), function(x) x[-1]), function(x) paste(x, collapse=" ")))
	#} else{
	#	dat.path$Category <- unlist(lapply(str_split(dat.path$Category, "_"), function(x) paste(x, collapse=" ")))
	#}

	# Get Zscores -------
	pathz <- matrix(0, nrow=nrow(dat.path), ncol=ncol(expr), dimnames=list(dat.path$Category, colnames(expr)))
	for(i in 1:nrow(dat.path)){
		pathway <- dat.path$Category[i]
		genes <- str_split(dat.path$overlap.genes[i], ",")[[1]]

		expr.path <- subset(expr, rownames(expr) %in% genes)
		exprz <- getZscore(expr.path)
		pathz[i,] <- apply(exprz, 2, mean)
	}
	return(pathz)
}


#### PLOT HEATMAP ----------------------
plot.heatmap.pathway.zscore <- function(dat, database, group.names, file.plot, plot.height, plot.width){
	# LOAD LIBRARIES -----
	require("gplots")
	require("RColorBrewer")

	group1 <- unique(group.names)[1]
	group2 <- unique(group.names)[2]
	group3 <- unique(group.names)[3]

	label.groups <- c( rep("#4ECDC4", length(which(group.names == group1))), rep("#C44D58", length(which(group.names == group2))), rep("#fdae6b", length(which(group.names == group3))) )


	# Generate Color Palette ------
	color_values <- c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#ffffff","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695") #red-blue
	color_palette <- colorRampPalette(color_values)(11)
	jColFun <- colorRampPalette(color_palette)

	# Generate Plot ------
	pdf(file.plot, height=plot.height, width=plot.width)
	heatmap.2(dat, 
		col = rev(jColFun(256)),
		Colv=TRUE, Rowv=TRUE, 
		#main=database, 
		dendrogram ="both", trace="none",  scale="none",
		cexCol=0.8, cexRow=0.5, margins = c(8,30), 
		hclustfun = function(x) hclust(x, method = "ward.D2"), 
		distfun = function(x) dist(x, method = "euclidean"),
		colsep=c(1:10000), rowsep=c(1:10000),
		sepcolor="black", sepwidth=c(0.000005,0.000005),
		ColSideColors=label.groups,
		key="TRUE", keysize=0.5, density.info="none", symkey=0,
		key.title=NA, # no title
		key.xlab=NA,  # no xlab
		key.par=list(mgp=c(1.5, 0.5, 0),
                       mar=c(2.5, 2.5, 1, 0)))
		legend("topright", legend=c("Normal","Tumor","Celline"), fill=c("#4ECDC4","#C44D58","#fdae6b"), border=FALSE, bty="n", x.intersp = 1, y.intersp = 1, cex=0.5)
	dev.off()
}

