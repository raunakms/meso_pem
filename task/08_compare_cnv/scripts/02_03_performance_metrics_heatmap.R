### LOAD LIBRAIES ---
library("stringr")
library("gplots")
library("RColorBrewer")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/08_compare_cnv")
dir.data <- file.path(dir.wrk, "data")
dir.output <- file.path(dir.wrk, "output")
dir.result <- file.path(dir.wrk, "results")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")

### DEFINE FILE ---
file.dat <- file.path(dir.output, "cnv_meso_ov_diff_feature_performance_summary_allgenes.tsv")

### LOAD DATA ---
dat <- read.delim(file.dat, header=T, stringsAsFactors=F, row.names=1)[,1:9]

# COLOR ---
jColFun <- colorRampPalette(brewer.pal(n = 9, "Reds"))

# Generate Plot ------
file.plot <- file.path(dir.plot, "heatmap_performance_metric.pdf")
pdf(file.plot, height=6, width=6)
	heatmap.2(as.matrix(dat), 
		col = jColFun(256),
		Colv=TRUE, Rowv=TRUE, 
		dendrogram ="both", trace="none",  scale="none",
		cexCol=1, cexRow=0.01, margins = c(10,5), 
		hclustfun = function(x) hclust(x, method = "ward.D2"), 
		distfun = function(x) dist(x, method = "euclidean"),
		#colsep=c(1:1000), rowsep=c(1:1000),
		#sepcolor="black", sepwidth=c(0.000005,0.000005),
		key="TRUE", keysize=1, density.info="none", symkey=0,
		key.title=NA, # no title
		key.xlab=NA,  # no xlab
		key.par=list(mgp=c(1.5, 0.5, 0),
                       mar=c(2.5, 2.5, 1, 0)))
dev.off()

### FINAL PASS 30 ---
file.dat <- file.path(dir.output, "cnv_meso_ov_diff_feature_performance_summary_finalpass30.tsv")
dat <- read.delim(file.dat, header=T, stringsAsFactors=F, row.names=1)[,1:9]
jColFun <- colorRampPalette(brewer.pal(n = 9, "Reds"))

# Generate Plot ------
file.plot <- file.path(dir.plot, "heatmap_performance_metric_final_pass_30.pdf")
pdf(file.plot, height=6, width=6)
	heatmap.2(as.matrix(dat), 
		col = jColFun(256),
		Colv=TRUE, Rowv=TRUE, 
		dendrogram ="both", trace="none",  scale="none",
		cexCol=1, cexRow=1, margins = c(5,8), 
		hclustfun = function(x) hclust(x, method = "ward.D2"), 
		distfun = function(x) dist(x, method = "euclidean"),
		colsep=c(1:1000), rowsep=c(1:1000),
		sepcolor="black", sepwidth=c(0.05,0.05),
		key="TRUE", keysize=1, density.info="none", symkey=0,
		key.title=NA, # no title
		key.xlab=NA,  # no xlab
		key.par=list(mgp=c(1.5, 0.5, 0),
                       mar=c(2.5, 2.5, 1, 0)))
dev.off()

### GET PAIRWISE CORRELATION ----
panel.cor <- function(x, y, digits=2, cex.cor=0.6){
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	r <- abs(cor(x, y, use = "na.or.complete"))
	txt <- format(c(r, 0.123456789), digits=digits)[1]
	test <- cor.test(x,y)
	Signif <- ifelse(round(test$p.value,3)<0.001,"p<0.001",paste("p=",round(test$p.value,3)))  
	text(0.5, 0.25, paste("r=",txt))
	text(.5, .75, Signif)
}

panel.hist <- function(x, ...){
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(usr[1:2], 0, 1.5) )
	h <- hist(x, plot = FALSE)
	breaks <- h$breaks; nB <- length(breaks)
	y <- h$counts; y <- y/max(y)
	rect(breaks[-nB], 0, breaks[-1], y, col="gray50", ...)
}

panel.smooth <- function (x, y, col="gray30", bg = NA, pch = 20, cex = 0.8, col.smooth = "red", span = 2/3, iter = 3, ...) {
	points(x, y, pch = pch, col = col, bg = bg, cex = cex)
	ok <- is.finite(x) & is.finite(y)
	if (any(ok))
	{lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), col = col.smooth, ...)}
}


### PLOT CORRELATION ---
file.plot <- file.path(dir.plot, "PlotPairs_performance_metric.pdf")
pdf(file.plot, height=6.8, width=6.8)
pairs(dat, lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist, cex.labels=1)
dev.off()
