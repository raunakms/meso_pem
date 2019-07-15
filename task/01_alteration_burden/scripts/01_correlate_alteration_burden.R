### SetDirectories ----------------------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/01_alteration_burden")
dir.data <- file.path(dir.wrk, "data")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")

### Define Files ------------------------------------------------------------
file.mut <- file.path(dir.data, "somatic_mutation_non_silent_samplewise_stats.tsv")
file.cna <- file.path(dir.data, "cnv_samplewise_stats.tsv")
file.fus <- file.path(dir.data, "fusion_samplewise_stat.tsv")
file.out <- file.path(dir.data, "expr_outlier_genes_frequency.tsv")

### Load Libraries ----------------------------------------------------------
library("stringr")
#library("ggtern")
library("scatterplot3d")
library("reshape2")
library("ggplot2")
library("gridExtra")

### Load Files ----------------------------------------------------------------
dat.mut <- read.delim(file.mut, header=T, stringsAsFactors=F)[,1:2]
dat.cna <- read.delim(file.cna, header=T, stringsAsFactors=F)[,1:4]
dat.fus <- read.delim(file.fus, header=T, stringsAsFactors=F)
dat.out <- read.delim(file.out, header=T, stringsAsFactors=F)
colnames(dat.out)[2] <- "FreqOut"



### Plot Alteration Data ------------------------------------------------------
ids <- unique(c(dat.mut$SampleID, dat.cna$SampleID, dat.fus$SampleID))[-20]

dat <- data.frame(SampleID=ids, FreqMUT=0, FreqCNA=0, FreqFUS=0, FreqOUT=0)
dat$SampleID <- as.character(dat$SampleID)

get.freq <- function(dat, df, item){
	for(i in 1:nrow(df)){
		index <- which(dat$SampleID == df$SampleID[i])
		if(length(index) == 0) next
		dat[index, which(colnames(dat) == item)] <- as.numeric(df[i,2])
	}
	return(dat)
}

dat <- get.freq(dat=dat, df=dat.mut, item="FreqMUT")
dat <- get.freq(dat=dat, df=dat.cna, item="FreqCNA")
dat <- get.freq(dat=dat, df=dat.fus, item="FreqFUS")
dat <- get.freq(dat=dat, df=dat.out, item="FreqOUT")

write.table(dat, file.path(dir.data, "meso_sample_alteration_stats.tsv"), sep="\t", row.names=F, col.names=T, quote=F)

dat <- dat[order(dat$FreqMUT, decreasing=T),]
df <- melt(dat, id="SampleID")
colnames(df) <- c("SampleID","Event","Freq")
df$SampleID <- as.character(df$SampleID)
df$Event <- as.character(df$Event)
df$Event[which(df$Event == "FreqMUT")] <- "MUTATION"
df$Event[which(df$Event == "FreqCNA")] <- "CNA"
df$Event[which(df$Event == "FreqFUS")] <- "GENE-FUSION"
df$Event[which(df$Event == "FreqOUT")] <- "GENE-EXPRESSION OUTLIERS"

df$SampleID <- factor(df$SampleID, levels=dat$SampleID)					 
df$Event <- factor(df$Event, levels=c("MUTATION","CNA","GENE-FUSION","GENE-EXPRESSION OUTLIERS"))					 
				
p1 <- ggplot(df, aes(y=Freq, x=SampleID)) +
		geom_bar(stat="Identity") +
		facet_wrap(~ Event, scales="free_y", nrow=4, ncol = 1) +
		theme(
			axis.text.x = element_text(size = 7, color="black", angle=90, hjust=1, vjust=0.5),
			axis.text.y = element_text(size = 7, color="black"),
			axis.title = element_text(size = 7, color="black"),
			plot.title = element_text(size = 10, color="black"),
			axis.ticks = element_line(size=0.4),
			legend.position = "none") +
		ylab("No. of Events") +
		xlab("") + 
		ggtitle("")

file.plot <- file.path(dir.plot, "alterations_stat.pdf")
pdf(file.plot, width=7, height=6)		
	grid.arrange(p1, ncol=1, nrow=1)
dev.off()


### 3D Scatter Plot  ------------------------------------------------------
file.plot <- file.path(dir.plot, "scatterplot_alterations.pdf")
pdf(file.plot, width=6, height=6)
with(dat, {
   s3d <- scatterplot3d(FreqMUT,   # x axis
						FreqCNA,   # y axis
						FreqFUS,   # z axis
						color="darkred", 
						pch=19, type="h", lty.hplot=2,
						main="",  
						xlab="Mutation",
						ylab="CNA",
						zlab="Gene Fusion")
	s3d.coords <- s3d$xyz.convert(FreqMUT, FreqCNA, FreqFUS)		
	text(s3d.coords$x, s3d.coords$y, labels=unlist(lapply(str_split(dat$SampleID, "-"), function(x) x[2])), cex=.8, pos=4) 	
})
dev.off()

### Compute Correlation ---------------------------------------------------------
get.correlation <- function(dat1, dat2, x.annot, y.annot){
	dat <- merge(dat1, dat2, by="SampleID")
	dat$slabel <- unlist(lapply(str_split(dat$SampleID, "-"), function(x) x[2]))

	if(colnames(dat)[2] == "FreqMut"){
		x.label <- "No. of Mutations"
	} else if(colnames(dat)[2] == "FreqCNV"){
		x.label <- "No. of CNA"
	} else if(colnames(dat)[2] == "FreqFusion"){
		x.label <- "No. of Gene Fusions"
	} else{
		x.label <- "No. of Outlier Genes"
	}

	if(colnames(dat)[3] == "FreqMut"){
		y.label <- "No. of Mutations"
	} else if(colnames(dat)[3] == "FreqCNV"){
		y.label <- "No. of CNA"
	} else if(colnames(dat)[3] == "FreqFusion"){
		y.label <- "No. of Gene Fusions"
	} else{
		y.label <- "No. of Outlier Genes"
	}
	
	colnames(dat) <- c("SampleID","FreqX","FreqY","slabel")
	r <- cor(dat$FreqX, dat$FreqY, use="complete.obs", method="pearson") 
	
	p <- ggplot(dat, aes(x=FreqX, y=FreqY)) + 
			geom_point(size=2, color="grey30") +
			stat_smooth(method="lm", geom = "smooth", formula = y ~ x, position = "identity", fullrange = FALSE, se=TRUE) +
			geom_text(aes(label=slabel, alpha=0.3, hjust=0, vjust=0, lineheight=1), color="black") +
			annotate("text", label = paste("R = ", round(r, 2), sep=""), x = x.annot, y = y.annot, size = 3, colour = "red") +
			theme(
				axis.text = element_text(size = 8, color="black"),
				axis.title = element_text(size = 8, color="black"),
				plot.title = element_text(size = 8, color="black"),
				legend.position="none") +
			xlab(x.label) + 
			ylab(y.label) + 
			ggtitle("") 	

	return(p)
}


p1 <- get.correlation(dat1=dat.mut, dat2=dat.out, x.annot=40, y.annot=800)
p2 <- get.correlation(dat1=dat.cna[,1:2], dat2=dat.out, x.annot=200, y.annot=800)
p3 <- get.correlation(dat1=dat.fus, dat2=dat.out, x.annot=20, y.annot=800)
p4 <- get.correlation(dat1=dat.mut, dat2=dat.cna[,1:2], x.annot=30, y.annot=250)
p5 <- get.correlation(dat1=dat.mut, dat2=dat.fus, x.annot=25, y.annot=20)
p6 <- get.correlation(dat1=dat.cna[,1:2], dat2=dat.fus, x.annot=200, y.annot=20)

file.plot <- file.path(dir.plot, "alteration_outlier_correlation.pdf")
pdf(file.plot, height=4, width=12)
	grid.arrange(p1, p2, p3, ncol=3, nrow=1)
	grid.arrange(p4, p5, p6, ncol=3, nrow=1)
dev.off()


