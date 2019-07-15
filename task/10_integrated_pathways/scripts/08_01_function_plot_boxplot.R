### LOAD LIBRARIES  -----
library("stringr")
library("reshape2")
library("ggplot2")

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

	return(list(dat.cna, grp0, grp1))
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

	return(list(expr.mrna, grp0, grp1))
}

### FUNCTION: LOAD PROTEIN EXPRESSION ---
loadProteinExpr <- function(file.dat, file.des){
	# LOAD DATA: Protein ---
	expr.prot <- read.delim(file.dat, header=TRUE, stringsAsFactors=FALSE, row.names=1)
	colnames(expr.prot) <- str_replace_all(colnames(expr.prot), "[.]", "-")
	expr.prot[is.na(expr.prot)] <- 0

	# LOAD DESIGN TABLE: Protein ---
	des.prot <- read.delim(file.des, header=T, stringsAsFactors=F)
	des.prot <- subset(des.prot, des.prot$SampleType == "Tumor")
	grp0 <- des.prot$SequencingID[which(des.prot$Group.3p21 == 0)]
	grp1 <- des.prot$SequencingID[which(des.prot$Group.3p21 == 1)]

	# ARRANGE EXPRESSION DATA ---
	expr.prot <- subset(expr.prot, select=c(grp0, grp1))
	colnames(expr.prot) <- str_replace_all(colnames(expr.prot), "T", "")

	grp0 <- str_replace_all(grp0, "T", "")
	grp1 <- str_replace_all(grp1, "T", "")

	return(list(expr.prot, grp0, grp1))
}

### FUNCTION MAIN: GET PROFILES ---
getProfile <- function(dir.plot, file.cna, file.mrna, file.prot, genes, type){

	# LOAD OMICS DATA---
	if(type == "CNA"){
		list.dat <- loadCNA(file.dat=file.cna, file.des=file.des.mrna)
		dat <- list.dat[[1]]
		grp0 <- list.dat[[2]]
		grp1 <- list.dat[[3]]
	}else if(type == "RNA"){
		list.dat <- loadmRNAExpr(file.dat=file.mrna, file.des=file.des.mrna)
		dat <- list.dat[[1]]
		grp0 <- list.dat[[2]]
		grp1 <- list.dat[[3]]
	}else if(type == "PROTEIN"){
		list.dat <- loadProteinExpr(file.dat=file.prot, file.des=file.des.prot)
		dat <- list.dat[[1]]
		grp0 <- list.dat[[2]]
		grp1 <- list.dat[[3]]		
	}

	# PREPARE DATA ---
	df <- prepareData(dat, genes, grp0, grp1)

	# GET BOXPLOT ---
	p <- get.boxplot(df)

	return(p)
}

### PREPARE DATA ---
prepareData <- function(dat, genes, grp0, grp1){
	# GET COMPLEX ---
	dat <- subset(dat, rownames(dat) %in% genes)
	genes <- subset(genes, genes %in% rownames(dat))

	# MELT DATA ---
	df <- melt(as.matrix(t(dat)))
	colnames(df) <- c("SampleID","Gene","Value")
	df$SampleID <- as.character(df$SampleID)
	df$Gene <- as.character(df$Gene)

	# ADD GROUP ---
	df$Group <- ""
	df$Group[which(df$SampleID %in% grp0)] <- "BAP1-Intact"
	df$Group[which(df$SampleID %in% grp1)] <- "BAP1-Del"

	# FACTORIZE DATA ---
	df$Group <- factor(df$Group, levels=c("BAP1-Intact","BAP1-Del"))
	df$Gene <- factor(df$Gene, levels=genes)

	return(df)
}

### FUNCTION: ----
get.boxplot <- function(df){
	require("ggplot2")
	cbPalette <- c("#fffc00","blue")

	# PLOT ---
	p <- ggplot(df, aes(x=Group, y=Value)) +
			geom_boxplot(aes(fill=Group, alpha=0.3), lwd=0.2, color="black", outlier.size=0.1, outlier.alpha=0.3, notch=FALSE) +
			scale_fill_manual(values=cbPalette) +
			#geom_jitter(width=0.1, color="black", alpha=0.5, size=0.5) +
			facet_wrap(~ Gene, nrow=1) +
			theme(
				axis.text.x = element_text(size=5, color="black", angle=90, hjust=0, vjust=0),
				axis.text.y = element_text(size=7, color="black"),
				axis.title = element_text(size=7, color="black"),
				plot.title = element_text(size=10, color="black", hjust=0),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),			
				axis.ticks = element_line(size=0.1, color="black"),
				strip.text = element_text(size=3, color="black"),
				strip.background = element_rect(fill="white", color="white"),
				panel.background = element_rect(fill="white", color="black", size=0.3),
				legend.position="none") + 
			ylab("Expression") +
			xlab("") + ggtitle("") 

	return(p)	
}

