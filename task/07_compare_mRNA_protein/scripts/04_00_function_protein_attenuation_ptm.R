### Load Libraries -----
library("stringr")
library("ggplot2")
library("ggExtra")
library("gridExtra")

### FUNCTION: LOAD CNA ---
loadCNA <- function(file.dat, file.des){
	# LOAD CNV DATA ---
	dat.cna <- read.delim(file.dat, header=T, stringsAsFactors=F, row.names=1)
	colnames(dat.cna) <- str_replace_all(colnames(dat.cna), "[.]", "-")

	# LOAD DESIGN TABLE: CNV ---
	des.mrna <- read.delim(file.des, header=TRUE, stringsAsFactors=FALSE)
	grp0 <- intersect(colnames(dat.cna), des.mrna$Sample.ID[which(des.mrna$Group.3p21 == 0)])
	grp1 <- intersect(colnames(dat.cna), des.mrna$Sample.ID[which(des.mrna$Group.3p21 == 1)])

	# ARRANGE SAMPLES: CNV ---
	dat.cna <- subset(dat.cna, select=c(grp0, grp1))

	return(dat.cna)
}

### FUNCTION: LOAD mRNA EXPRESSION ---
loadmRNAExpr <- function(file.dat, file.des){
	# LOAD DATA: mRNA ---
	expr.mrna <- read.delim(file.dat, header=TRUE, stringsAsFactors=FALSE, row.names=1)
	colnames(expr.mrna) <- str_replace_all(colnames(expr.mrna), "[.]", "-")
	expr.mrna[is.na(expr.mrna)] <- 0

	# LOAD DESIGN TABLE: mRNA ---
	des.mrna <- read.delim(file.des, header=TRUE, stringsAsFactors=FALSE)
	grp0 <- intersect(colnames(expr.mrna), des.mrna$Sample.ID[which(des.mrna$Group.3p21 == 0)])
	grp1 <- intersect(colnames(expr.mrna), des.mrna$Sample.ID[which(des.mrna$Group.3p21 == 1)])

	# ARRANGE SAMPLES: mRNA ---
	expr.mrna <- subset(expr.mrna, select=c(grp0, grp1))

	return(expr.mrna)
}

### FUNCTION: LOAD mRNA EXPRESSION ---
loadProteinExpr <- function(file.dat, file.des){
	# LOAD DATA: Protein ---
	expr.prot <- read.delim(file.dat, header=TRUE, stringsAsFactors=FALSE, row.names=1)
	colnames(expr.prot) <- str_replace_all(colnames(expr.prot), "[.]", "-")
	expr.prot[is.na(expr.prot)] <- 0
	expr.prot[expr.prot == "-Inf"] <- 0

	# LOAD DESIGN TABLE: Protein ---
	des.prot <- read.delim(file.des, header=T, stringsAsFactors=F)
	des.prot <- subset(des.prot, des.prot$SampleType == "Tumor")
	grp0 <- des.prot$SequencingID[which(des.prot$Group.3p21 == 0)]
	grp1 <- des.prot$SequencingID[which(des.prot$Group.3p21 == 1)]

	# ARRANGE EXPRESSION DATA ---
	expr.prot <- subset(expr.prot, select=c(grp0, grp1))
	colnames(expr.prot) <- str_replace_all(colnames(expr.prot), "T", "")

	return(expr.prot)
}


### FUNCTION: getDesignTable ---
getDesignTable <- function(dat.cna, expr.mrna, expr.prot, file.des.mrna){
	ids <- Reduce(intersect, list(colnames(dat.cna),colnames(expr.mrna),colnames(expr.prot))) 
	des <- read.delim(file.des.mrna, header=TRUE, stringsAsFactors=FALSE)
	des <- subset(des, des$Sample.ID %in% ids)
	des <- des[order(des$Group.3p21, decreasing=FALSE),]
	#ids <- des$Sample.ID	
	return(des)
}


### FUNCTION: prepareData ---
prepareData <- function(dat, des, genes, type){
	# SUBSET AND RESHAPE DATA  ----
	dat <- subset(dat, rownames(dat) %in% genes)
	
	# ORDER GENES ---
	dat <- dat[match(genes, rownames(dat)),]	

	# MELT DATA ---
	df <- melt(t(as.matrix(dat)))
	colnames(df) <- c("SampleID","Gene","Value")	

	# ADD ATTRIBUTES ---
	df$SampleID <- as.character(df$SampleID)
	df$Gene <- as.character(df$Gene)
	df$Group <- ""
	df$Type <- type

	# GET GROUPS ---
	grp0 <- des$Sample.ID[which(des$Group.3p21 == 0)]
	grp1 <- des$Sample.ID[which(des$Group.3p21 == 1)]
	df$Group[which(df$SampleID %in% grp0)] <- "BAP1-INTACT"
	df$Group[which(df$SampleID %in% grp1)] <- "BAP1-DEL"

	return(df)
}

### FUNCTION: combineData ---
combineData <- function(dat1, dat2, ids, genes){
	df <- data.frame(SampleID=dat1$SampleID,
					Gene=dat1$Gene,
					Group=dat1$Group,
					Type=paste(dat1$Type, dat2$Type, sep=":"),
					Value1=dat1$Value,
					Value2=dat2$Value)

	df$SampleID <- factor(df$SampleID, levels=ids)
	df$Gene <- factor(df$Gene, levels=genes)
	df$Group <- factor(df$Group, levels=c("BAP1-INTACT","BAP1-DEL"))
	df$Type <- as.factor(df$Type)

	return(df)
}

### FUNCTION: prepareComboData ---
prepareComboData <- function(dat1, dat2, des, type1, type2){
	require("reshape2")

	# GET COMMON GENES ---
	genes <- intersect(rownames(dat1), rownames(dat2))
	ids <- colnames(dat1)

	# GET DATA ---
	df1 <- prepareData(dat=dat1, des=des, genes=genes, type=type1)
	df2 <- prepareData(dat=dat2, des=des, genes=genes, type=type2)

	# COMBINE DATA ----
	df <- combineData(dat1=df1, dat2=df2, ids, genes)

	return(df)
}


### FUNCTION: getCorrelation ---
getCorrelation <- function(dat){
	genes <- unique(dat$Gene)
	type <- unique(dat$Type)

	# PREPARE DATA  ---
	dat.cor <- data.frame(Gene=genes)
	dat.cor$Gene <- as.character(dat.cor$Gene)
	dat.cor$Type <- type
	dat.cor$R <- 0
	#dat.cor$pvalue <- 0

	# GET CORRELATION ---
	for(i in 1:length(genes)){
		df <- subset(dat, dat$Gene == genes[i])

		x_na <- length(which(is.na(df$Value1)))
		y_na <- length(which(is.na(df$Value2)))

		x0 <- length(which(df$Value1 == 0))
		y0 <- length(which(df$Value2 == 0))

		# FOR ORIGINAL: ALL SAMPLES ---
		#if((x_na > 9) | (y_na > 9)) next
		#if((x0 > 9) | (y0 > 9)) next

		# FOR BY SUBTYPE ---
		if((x_na > 3) | (y_na > 3)) next
		if((x0 > 3) | (y0 > 3)) next

		correlation.test <- cor.test(x=df$Value1, y=df$Value2, method="pearson")
		dat.cor$R[i] <- as.numeric(correlation.test$estimate)
		#pvalue <- as.numeric(correlation.test$p.value)
	
		#if(pvalue == 0){
		#	dat.cor$pvalue[i] <- as.numeric(2.2e-16)
		#} else{
		#	dat.cor$pvalue[i] <- pvalue
		#}

		cat("CORRELATION", as.character(type), i, "OF", length(genes), sep="\t", "\n")
	}

	return(dat.cor)
}



### FUNCTION: ----
get.plot <- function(df){
	require("ggplot2")
	require("ggExtra")
	cbPalette1 <- c("#98a3a5","#e85748")
	cbPalette2 <- c("#98a3a5","#e85748","#000000")

	# PLOT ---
	p <- ggplot(df, aes(x=R1, y=R2, label=LabelGenes)) + 
			#geom_point(shape = 21, fill="white", color="black", size=2, stroke=1, alpha=0.9) +
			geom_point(aes(fill=as.factor(Rdiff.group), color=as.factor(TagGenes)), shape = 21, stroke=0.5, size=1, alpha=0.4) +
			scale_fill_manual(values=cbPalette1) +
			scale_color_manual(values=cbPalette2) +
			#stat_smooth(method="lm", geom = "smooth", formula = y ~ x, position = "identity", fullrange = FALSE, se=TRUE) +
			geom_density2d(aes(group=as.factor(Rdiff.group)), stat = "density2d", lineend = "round", linejoin = "round", color="#FFFFFF", alpha=0.9, size=0.25) +
			coord_cartesian(xlim=c(-1, 1), ylim=c(-1, 1)) +
			geom_abline(intercept = 0, linetype = 2, color="#969696") +
			#geom_text(size=1, color="#000000", hjust=0, vjust=0) +
			theme(
				axis.text = element_text(size = 7, color="black"),
				axis.title = element_text(size = 7, color="black"),
				strip.text = element_text(size = 7, color="black"),
				strip.background = element_rect(fill = "white", colour = "white"),
				plot.title = element_text(size = 7, color="black", hjust=0),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.ticks = element_line(size=0.4, color="black"),	
				panel.background = element_rect(fill="white", color="black"),
				legend.position="none") +
			ylab("Copy-number - Proteome") + 
			xlab("Copy-number - Transcriptome") + 
			ggtitle("")

	# add marginal density
	p <- ggExtra::ggMarginal(p, type="density", size = 4, aes(group=as.factor(Rdiff.group)), alpha=0.5) 

	return(p)	
}


### FOR MARGINAL DENSITY PLOT ----
get.densityplot <- function(df, type){
	require("ggplot2")
	cbPalette <- c("#98a3a5","#e85748")

	# SELECT VALUES ---
	colnames(df)[which(colnames(df) == type)] <- "R"

	g <- ggplot(df, aes(R, fill = as.factor(Rdiff.group))) + 
			geom_density(aes(color = as.factor(Rdiff.group)), alpha=0.4, size=0.5) +
			scale_fill_manual(values=cbPalette) +
			scale_color_manual(values=cbPalette) +
			coord_cartesian(xlim=c(-1, 1)) +
			theme(
				axis.text = element_text(size = 7, color="black"),
				axis.title = element_text(size = 7, color="black"),
				strip.text = element_text(size = 7, color="black"),
				strip.background = element_rect(fill = "white", colour = "white"),
				plot.title = element_text(size = 7, color="black", hjust=0),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.ticks = element_line(size=0.4, color="black"),	
				panel.background = element_rect(fill="white", color="white"),
				legend.position="none") +
			ylab("") +
			xlab("") + 
			ggtitle("")

	return(g)
}

