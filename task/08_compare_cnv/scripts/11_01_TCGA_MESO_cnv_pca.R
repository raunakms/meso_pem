### Load Libraries -----
library("stringr")

### DEFINE PATH ----
dir.tcga <- file.path("/Data/Raunak/HITnDRIVE/datasets/TCGA_MESO/data")
dir.cnv <- file.path(dir.tcga,"cnv/seq_call_refseq_genes_meso")
dir.des <- file.path(dir.tcga,"annotation")
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/08_compare_cnv")
dir.data <- file.path(dir.wrk, "data")
dir.output <- file.path(dir.wrk, "output")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")

### DEFINE FILES -----
file.dat <- file.path(dir.cnv, "tcga_meso_cnv_seg_values_calls_parsed.tsv.gz")
file.des <- file.path(dir.des, "design_table.tsv")

### LOAD DATA ---
dat <- read.delim(file.dat, header=T, stringsAsFactors=F, row.names=1)
colnames(dat) <- str_replace_all(colnames(dat), "[.]", "-")

### LOAD DES ---
des <- read.delim(file.des, header=T, stringsAsFactors=F)
des <- subset(des, des$SampleID %in% colnames(dat))
des <- des[match(colnames(dat), des$SampleID),]


### FUNCTION TO PLOT PCA ---
plot.pca <- function(dat, features, sample.class, file.plot){
	require("scatterplot3d")

	### SUBSET CNV DATA ---
	dat <- subset(dat, rownames(dat) %in% features)

	### PRINCIPAL COMPONENT ANALYSIS ---
	dat.pca <- prcomp(na.omit(dat))
	dat.eigenvectors <- as.data.frame(dat.pca $rotation[,1:3])

	### GET ITEMS ---
	groups <- unique(sample.class)

	### GET DATA GROUPS ---
	class1 <- which(sample.class == unique(sample.class)[1])
	class2 <- which(sample.class == unique(sample.class)[2])
	class3 <- which(sample.class == unique(sample.class)[3])
	class4 <- which(sample.class == unique(sample.class)[4])

	### GET COLOR ---
	color1 <- "#C44D58"
	color2 <- "#4ECDC4"
	color3 <- "#556270"
	color4 <- "#C7F464"

	label.color <- c(rep(color1, length(class1)), 
						rep(color2, length(class2)), 
						rep(color3, length(class3)), 
						rep(color4, length(class4)))

	color1234 <- c(color1,color2,color3,color4)

	### Plot PCA ---
	#pdf(file.plot, height=4, width=4)
	#plot(x=dat.eigenvectors$PC1, y=dat.eigenvectors$PC2,
	#	#xlim=c(-0.3,0.1), ylim=c(-0.2,0.3),
	#	type="p", pch=19, col=label.color,
	#	cex=0.6, cex.axis=0.5, cex.lab=0.5, cex.main=0.5,
	#	las=1, tck = -.03, xlab="PC1", ylab="PC2",
	#	main="PCA")
	#legend("bottom", legend=groups, pch=19, col=color1234, bty='n', cex=.5, horiz=TRUE, inset=c(0,0))
	#dev.off()

	### 3D Scatter Plot  ----
	pdf(file.plot, width=6, height=6)
		with(dat.eigenvectors, {
   			s3d <- scatterplot3d(PC1,   # x axis
						PC2,   # y axis
						PC3,   # z axis
						color=label.color, 
						pch=19, type="p", 
						#lty.hplot=2,
						main="",  
						xlab="PC1",
						ylab="PC2",
						zlab="PC3")
			})
		legend("top", legend=groups, pch=19, col=color1234, bty='n', cex=.5, horiz=TRUE)
	dev.off()
}

### CALL PLOT PCA FUNCTION ---
plot.pca(dat=dat, features=rownames(dat), sample.class=des$Subtype, 
		file.plot=file.path(dir.plot, "PCA_CNA_MESO_PM_allgenes.pdf"))


### FUNCTION: remove NA ---
removeNA <- function(dat){
	ctr <- 1
	del.index <- vector()
	for(i in 1:nrow(dat)){
		y <- as.numeric(dat[i,])
		ind <- which(is.na(y))
		if(length(ind) != 0){
			del.index[ctr] <- i
			ctr <- ctr + 1
		}
	}
	return(del.index)
}

### CALL FUNCTION ---
del.index <- removeNA(dat)
dat <- dat[-del.index,]


### FUNCTION: MELT DATA ---
MeltData <- function(dat, des){
	# LOAD LIBRARIES ---
	require("reshape2")

	# GET SAMPLE CLASS ---
	ids.class1 <- des$SampleID[which(des$Subtype == "Biphasic")] 
	ids.class2 <- des$SampleID[which(des$Subtype == "DiffuseMaglignant")]  
	ids.class3 <- des$SampleID[which(des$Subtype == "Epithelioid")] 
	ids.class4 <- des$SampleID[which(des$Subtype == "Sarcomatoid")] 

	# MELT DATA ---
	dat.cnv <- melt(as.matrix(t(dat)))
	colnames(dat.cnv) <- c("SampleID","Feature","value")

	dat.cnv$SampleID <- as.character(dat.cnv$SampleID)
	dat.cnv$Feature <- as.character(dat.cnv$Feature)

	# ADD GROUP LABEL ---
	dat.cnv$Group <- ""
	dat.cnv$Group[which(dat.cnv$SampleID %in% ids.class1)] <- "B"
	dat.cnv$Group[which(dat.cnv$SampleID %in% ids.class2)] <- "D"
	dat.cnv$Group[which(dat.cnv$SampleID %in% ids.class3)] <- "E"
	dat.cnv$Group[which(dat.cnv$SampleID %in% ids.class4)] <- "S"

	dat.cnv$Group <- factor(dat.cnv$Group, labels=c("B","D","E","S"))
	#dat.cnv$Feature <- as.factor(dat.cnv$Feature)
	#dat.cnv$SampleID <- as.factor(dat.cnv$SampleID)

	return(dat.cnv)
}

### FUNCTION: calcKruskalWallisTest_parallel --------
calcKruskalWallisTest_parallel <- function(df){
	# LOAD LIBRARIES ---
	require("foreach")
	require("doParallel")
	require("PMCMR")

	# GET ATTRIBUTES ---
	features <- as.character(unique(df$Feature))
	nfeatures <- length(features)

	# Declate Cluster ----
	no_cores <- 40
	cl <- makeCluster(no_cores)
	registerDoParallel(cl)

	# LOOP FOR EACH FEATURE ---
	dat.summary <- foreach(i = 1:nfeatures, .combine=rbind) %dopar%{	
						d <- subset(df, df$Feature == features[i])
						kwtest <- kruskal.test(value ~ Group, data = d, na.action=na.omit)
						#ph <-  posthoc.kruskal.nemenyi.test(x=d$value, g=d$Group, method="Tukey", na.action=na.omit)

						dat.summary <- data.frame(Feature=features[i],
											statistic=as.numeric(kwtest$statistic),
											pvalue=kwtest$p.value
											#phtest_D_B=ph$p.value[1,1],
											#phtest_E_B=ph$p.value[2,1],
											#phtest_S_B=ph$p.value[3,1],
											#phtest_E_D=ph$p.value[2,2],
											#phtest_S_D=ph$p.value[3,2],
											#phtest_S_E=ph$p.value[3,3]
											)
					}
	stopCluster(cl)

	dat.summary <- dat.summary[order(dat.summary$pvalue, decreasing=F),]
	dat.summary$fdr <- p.adjust(dat.summary$pvalue, method = "BH", n = length(dat.summary$pvalue))
	
	return(dat.summary)
}


### CALL FUNCTION: ---
df.kw <- calcKruskalWallisTest_parallel(df=MeltData(dat, des))
df.kw.pass <- subset(df.kw, df.kw$pvalue <= 0.01)
df.kw.pass$Feature <- as.character(df.kw.pass$Feature)

### WRITE OUTPUT ---
file.output <- file.path(dir.output, "TCGA-MESO_subtype_cna_KruskalWallisTest.tsv")
write.table(df.kw, file.output, sep="\t", row.names=F, col.names=T, quote=F)


### CALL PLOT PCA FUNCTION ---
plot.pca(dat=dat, features=df.kw.pass$Feature, sample.class=des$Subtype, 
		file.plot=file.path(dir.plot, "PCA_CNA_MESO_PM_kwpassgenes.pdf"))

