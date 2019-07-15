#### LOAD LIBRARIES -------------------
library("stringr")
library("ggplot2")
library("gridExtra")
library("RColorBrewer")

#### DEFINE PATH ----------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.proteome <- file.path(dir.wrk, "data/proteome/processed_data")
dir.task <- file.path(dir.wrk, "task/06_proteomics")
dir.data <- file.path(dir.task, "data")
dir.output <- file.path(dir.task, "output")
dir.plot <- file.path(dir.task, "plot")

#### DEFINE FILE ---------------------
file.expr1 <- file.path(dir.proteome, "meso_proteome_proteomediscover_all_log2.tsv.gz")
file.expr2 <- file.path(dir.proteome, "meso_proteome_ptm_all_log2.tsv.gz")
file.des <- file.path(dir.proteome, "design_table_proteome.tsv")
file.array <- file.path("/Data/Raunak/softwares/bdvtools/array_process/array_preprocess.R")

#### LOAD DESIGN TABLE ----------------
des <- read.delim(file.des, header=T, stringsAsFactors=F)
ids.normal <- des$SequencingID[which(des$SampleType == "Normal")]
ids.tumor <- des$SequencingID[which(des$SampleType == "Tumor")]
ids.celline <- des$SequencingID[which(des$SampleType == "Celline")]

### FUNCTION: getExpr -----
getExpr <- function(file.expr, file.des, ids){
	expr <- read.delim(file.expr, header=T, stringsAsFactors=F, row.names=1)
	colnames(expr) <- str_replace_all(colnames(expr), "[.]", "-")
	expr <- subset(expr, select=ids)
	expr[is.na(expr)] <- 0
	expr[expr == "-Inf"] <- 0

	# LOAD DESIGN TABLE: Protein ---
	des <- read.delim(file.des, header=T, stringsAsFactors=F)
	des <- subset(des, des$SampleType == "Tumor")
	grp0 <- des$SequencingID[which(des$Group.3p21 == 0)]
	grp1 <- des$SequencingID[which(des$Group.3p21 == 1)]

	# ARRANGE EXPRESSION DATA ---
	expr <- subset(expr, select=c(grp0, grp1))
	colnames(expr) <- str_replace_all(colnames(expr), "T", "")

	### REMOVE 0s ---
	y0 <- rep(0, nrow(expr))
	for(i in 1:nrow(expr)){
		y0[i] <- length(which(expr[i,] == 0))
	}

	del.index <- which(y0 == ncol(expr))
	if(length(del.index) != 0){
		expr <- expr[-del.index,]
	}

	return(expr)
}

#### LOAD PROTEIN EXPRESSION DATA -----
expr.wt <- getExpr(file.expr=file.expr1, file.des, ids=ids.tumor)
expr.ptm <- getExpr(file.expr=file.expr2, file.des, ids=ids.tumor)

dim(expr.wt)
dim(expr.ptm)

#> dim(expr.wt)
#[1] 8239   16
#> dim(expr.ptm)
#[1] 8175   16


### GET COMMON GENELIST ---
genes.wt <- rownames(expr.wt)
genes.ptm <- rownames(expr.ptm)
list.genes <- list(wt=genes.wt, ptm=genes.ptm )
genes.common <- Reduce(intersect, list.genes)

#> length(genes.common)
#[1] 7428

### Subset genes ---
expr.wt <- subset(expr.wt, rownames(expr.wt) %in% genes.common)
expr.ptm <- subset(expr.ptm, rownames(expr.ptm) %in% genes.common)

### Match Genes order ---
expr.wt <- expr.wt[match(genes.common, rownames(expr.wt)),]
expr.ptm <- expr.ptm[match(genes.common, rownames(expr.ptm)),]

### MERGE DATA ---
expr <- cbind(expr.wt, expr.ptm)
y <- colnames(expr.wt)

### QUANTILE NORM ---
#source(file.array)
#dqnorm <- getQuantile.normalize(expr)
#expr.wt <- dqnorm[,1:16]
#expr.ptm <- dqnorm[,17:32]
#colnames(expr.wt) <- y
#colnames(expr.ptm) <- y

### FUNCTION: getProteinCorr ----
getProteinCorr <- function(expr1, expr2, method){
	dat.cor <- data.frame(Gene=rownames(expr1))
	dat.cor$Gene <- as.character(dat.cor$Gene)
	dat.cor$R <- NA
	dat.cor$pvalue <- NA

	for(i in 1:nrow(expr1)){
		x <- as.numeric(expr1[i,])
		y <- as.numeric(expr2[i,])

		x0 <- length(which(x == 0))
		y0 <- length(which(y == 0))

		if((x0 > 12) | (y0 > 12)) next

		correlation.test <- cor.test(x, y, method=method)
		dat.cor$R[i] <- as.numeric(correlation.test$estimate)
		pvalue <- as.numeric(correlation.test$p.value)
	
		if(pvalue == 0){
			dat.cor$pvalue[i] <- as.numeric(2.2e-16)
		} else{
			dat.cor$pvalue[i] <- pvalue
		}
	}
	dat.cor$fdr <- p.adjust(dat.cor$pvalue, method="BH", n=length(dat.cor$pvalue))
	dat.cor <- dat.cor[order(dat.cor$R, decreasing=T),]
	dat.cor$Group <- ifelse(dat.cor$R >= 0.5, "Positive", ifelse((dat.cor$R < 0.5) & (dat.cor$R > 0),  "Modrate", ifelse((dat.cor$R <= 0) & (dat.cor$R > -0.5),  "Weak", "Negative")))
	dat.cor$Color <- ifelse(dat.cor$R >= 0.5, "green", ifelse((dat.cor$R < 0.5) & (dat.cor$R > 0),  "yellow", ifelse((dat.cor$R <= 0) & (dat.cor$R > -0.5),  "grey80", "red")))

	return(dat.cor)
}


### FUNCTION: get.plot ----
get.plot <- function(dat, method, analysis){
	del.index <- which(is.na(datp$R) == TRUE)
	if(length(del.index) != 0){
		dat <- dat[-del.index,]
	}

	title.name <- paste("protein WildType vs PTM -", method, "correlation distribution", sep=" ")
	val.ylim <- c(0,250)

	#items <- c("Positive correlation", "Moderate correlation", "Moderately Weak correlation", "Negative correlation")
	h <- hist(dat$R, breaks=100, plot=FALSE) 
	xfit <- seq(min(dat$R), max(dat$R), length=40) 
	yfit <- dnorm(xfit, mean=mean(dat$R), sd=sd(dat$R)) 
	yfit <- yfit*diff(h$mids[1:2])*length(dat$R)
	cuts <- cut(h$breaks, c(-Inf,-0.5,0,0.5,Inf)) 
	plot(h, col=c("red","grey80","yellow","green")[cuts],
		 main=title.name,
		 xlab=paste(method, "correlation coefficient (R)", sep=" "), ylab="No. of Genes",
		 xlim=c(-1,1), ylim=val.ylim,
		 cex.main=0.8, cex.lab=0.7, cex.axis=0.7, las=1, tck=-0.03)
	#legend("topright", items, fill=c("green","yellow","grey80","red"), bty="n", cex=0.7)
	lines(xfit, yfit, col="black", lwd=2)
}

### PROTEIN-PROTEIN CORRELATION ---
datp <- getProteinCorr(expr1=expr.wt, expr2=expr.ptm, method="pearson")
dats <- getProteinCorr(expr1=expr.wt, expr2=expr.ptm, method="spearman")

### PLOT ---
file.plot <- file.path(dir.plot, "protome_wt_ptm_correlation_rho_distribution.pdf")
pdf(file.plot, height=6, width=6)
par(mfrow=c(2,1))
	get.plot(datp, method="Pearson")
	get.plot(dats, method="Spearman")
dev.off()

### FUNCTION: getCorrelationPlot ---
getCorrPlot <- function(df){
	id <- as.character(unique(df$SampleID))

	#GENERATE PLOT -----
	p <- ggplot(df, aes(x=Value.wt, y=Value.ptm)) + 
			geom_point(color="black", alpha=0.5, size=0.1) +
			#scale_color_manual(values=c("#bdbdbd","#bd0026")) +
			#geom_text(aes(x=Expr.rnaseq, y=Expr.marray, label=Label), size=1, color="black", hjust=0, vjust=0)+
			#stat_smooth(method="lm", geom = "smooth", formula = y ~ x, position = "identity", fullrange = FALSE, se=FALSE) +
			geom_density2d(stat = "density2d", lineend = "round", linejoin = "round", alpha=0.9, color="yellow", size=0.25) +
			#coord_cartesian(xlim=c(0, 25), ylim=c(0, 25)) +
			theme(
				axis.text = element_text(size = 5, color="black"),
				axis.title = element_text(size = 8, color="black"),
				#strip.text = element_text(size = 10, color="black"),
				#strip.background = element_rect(fill = "white", colour = "white"),
				plot.title = element_text(size = 8, color="black", hjust=0),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.ticks = element_line(size=0.4, color = "black"),	
				panel.background = element_rect(fill = "white", colour = "black"),
				legend.position="none") +
			ylab("Expression Protein PTM") + 
			xlab("Expression Protein WT") + 
			ggtitle(id) 

	return(p)			
}

### SAMPLEWISE CORRELATION  ---
list.plot <- list()
for(i in 1:ncol(expr.wt)){
	genes <- rownames(expr.wt)
	id <- colnames(expr.wt)[i]
	x <- as.numeric(expr.wt[,i])
	y <- as.numeric(expr.ptm[,i])

	# PREPARE DATA ---
	df <- data.frame(SampleID=id, Gene=genes, Value.wt=x, Value.ptm=y)
	df$SampleID <- as.factor(df$SampleID)
	df$Gene <- as.factor(df$Gene)

	# GET PLOT ---
	list.plot[[i]] <- getCorrPlot(df)
}


### PLOT ---
file.plot <- file.path(dir.plot, "correlation_samplewise_protein_wt_ptm.pdf")
pdf(file.plot, height=4, width=4)
	grid.arrange(list.plot[[1]], list.plot[[2]], list.plot[[3]], list.plot[[4]], ncol=2, nrow=2)
	grid.arrange(list.plot[[5]], list.plot[[6]], list.plot[[7]], list.plot[[8]], ncol=2, nrow=2)
	grid.arrange(list.plot[[9]], list.plot[[10]], list.plot[[11]], list.plot[[12]], ncol=2, nrow=2)
	grid.arrange(list.plot[[13]], list.plot[[14]], list.plot[[15]], list.plot[[16]], ncol=2, nrow=2)
dev.off()

