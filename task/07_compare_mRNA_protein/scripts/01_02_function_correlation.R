get.correlation <- function(method, genes.common, expr1, expr2, genes.cgc=NULL){
	dat.cor <- data.frame(Gene=genes.common)
	dat.cor$Gene <- as.character(dat.cor$Gene)
	dat.cor$R <- 0
	dat.cor$pvalue <- 0

	for(i in 1:length(genes.common)){
		x <- as.numeric(expr1[i,])
		y <- as.numeric(expr2[i,])

		x_na <- length(which(is.na(x)))
		y_na <- length(which(is.na(y)))

		x0 <- length(which(x == 0))
		y0 <- length(which(y == 0))

		if((x_na > 8) | (y_na > 8)) next
		if((x0 > 8) | (y0 > 8)) next

		#correlation.test <- cor.test(x, y, method=method, exact=TRUE, conf.level=0.95, continuity=TRUE)
		correlation.test <- cor.test(x, y, method=method)
		dat.cor$R[i] <- as.numeric(correlation.test$estimate)
		pvalue <- as.numeric(correlation.test$p.value)
	
		if(pvalue == 0){
			dat.cor$pvalue[i] <- as.numeric(2.2e-16)
		} else{
			dat.cor$pvalue[i] <- pvalue
		}

		#cat("PROCESSED", i, "OF", length(genes.common), "\n", sep="\t")
	}
	dat.cor$fdr <- p.adjust(dat.cor$pvalue, method="BH", n=length(dat.cor$pvalue))
	dat.cor <- dat.cor[order(dat.cor$R, decreasing=T),]
	dat.cor$Group <- ifelse(dat.cor$R >= 0.5, "Positive", ifelse((dat.cor$R < 0.5) & (dat.cor$R > 0),  "Modrate", ifelse((dat.cor$R <= 0) & (dat.cor$R > -0.5),  "Weak", "Negative")))
	dat.cor$Color <- ifelse(dat.cor$R >= 0.5, "green", ifelse((dat.cor$R < 0.5) & (dat.cor$R > 0),  "yellow", ifelse((dat.cor$R <= 0) & (dat.cor$R > -0.5),  "grey80", "red")))

	if(!is.null(genes.cgc)){
		dat.cor$Gene.CGC <- 0
		dat.cor$Gene.CGC[which(dat.cor$Gene %in% genes.cgc)] <- 1
	}

	return(dat.cor)
}

get.plot <- function(dat, method, analysis){
	if(analysis == "mrna_protein"){
		title.name <- paste("protein vs mRNA expression -", method, "correlation distribution", sep=" ")
		val.ylim <- c(0,200)
	} else if(analysis == "mrna_cnv"){
		title.name <- paste("Copy Number vs mRNA expression -", method, "correlation distribution", sep=" ")
		val.ylim <- c(0,500)
	} else if(analysis == "protein_cnv"){
		title.name <- paste("Copy Number vs protein expression -", method, "correlation distribution", sep=" ")
		val.ylim <- c(0,500)
	} else{
		title.name <- ""
		val.ylim <- c(0,20)
	}

	items <- c("Positive correlation", "Moderate correlation", "Moderately Weak correlation", "Negative correlation")
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

get.cgcplot <- function(dat, method, analysis){
	dat <- subset(dat, dat$Gene.CGC == 1)
	method <- method
	analysis <- analysis
	return(get.plot(dat, method, analysis))
}

### CALCULATE STATS ---
get.corr.stat <- function(df, r.threshold, direction){
	n <- nrow(df)

	if(direction == "POSITIVE"){
		df <- subset(df, df$R >= r.threshold)
	} else if(direction == "NEGATIVE"){
		df <- subset(df, df$R <= r.threshold)
	}

	x <- nrow(df)
	p <- round(x/n * 100,2)

	d <- data.frame(Items=c("R.Threshold","Direction","Genes.Total","Genes.Pass","Concordance.Percent"), 
					Value=c(r.threshold,direction,n,x,p))
	return(d)
}
