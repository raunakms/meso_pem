### LOAD LIBRAIES ---
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/08_compare_cnv")
dir.data <- file.path(dir.wrk, "data")
dir.output <- file.path(dir.wrk, "output")
dir.result <- file.path(dir.wrk, "results")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")

### DEFINE FILES ---
file.dat <- file.path(dir.data, "cnv_meso_ov_seg_mean_combined.tsv.gz")
file.des <- file.path(dir.data, "cnv_meso_ov_combined_class.txt")
file.df65 <- file.path(dir.output, "cnv_meso_ov_diff_feature_performance_summary_compilation.tsv")
file.features3899<- file.path(dir.output, "cnv_meso_ov_diff_feature_genelist_wcox_nexus.txt")

### FUNCTION TO PLOT PCA ---
plot.pca <- function(dat, features, sample.class, file.plot){
	### SUBSET CNV DATA ---
	dat <- subset(dat, rownames(dat) %in% features)

	### PRINCIPAL COMPONENT ANALYSIS ---
	dat.pca <- prcomp(na.omit(dat))
	dat.eigenvectors <- as.data.frame(dat.pca $rotation[,1:2])

	### COLOR DATA ---
	sample.class.meso <- which(sample.class == "MESO")
	sample.class.ov <- which(sample.class == "OV")

	label.color <- c(rep("red", length(sample.class.meso)), rep("blue", length(sample.class.ov)))

	### Plot PCA ---
	pdf(file.plot, height=4, width=4)
	plot(x=dat.eigenvectors$PC1, y=dat.eigenvectors$PC2,
		#xlim=c(-0.3,0.1), ylim=c(-0.2,0.3),
		type="p", pch=19, col=label.color,
		cex=0.5, cex.axis=0.5, cex.lab=0.5, cex.main=0.5,
		las=1, tck = -.03, xlab="PC1", ylab="PC2",
		main="PCA: MESO-OV Copy Number Profiles")
	dev.off()

}

### LOAD DESIGN TABLE ---
sample.class <- read.table(file.des, header=F, stringsAsFactors=F)$V1

### LOAD FEATURE TABLE ---
df65 <- read.delim(file.df65, header=T, stringsAsFactors=F)
features.65 <- df65$Gene

features.3899 <- read.table(file.features3899, header=F, stringsAsFactors=F)$V1

### LOAD CNV DATA ---
dat <-  read.delim(file.dat, header=T, stringsAsFactors=F, row.names=1)
colnames(dat) <- str_replace_all(colnames(dat), "[.]", "-")


### CALL PLOT PCA FUNCTION ---
plot.pca(dat=dat, features=features.3000, sample.class, file.plot=file.path(dir.plot, "PCA_CNA_OV_features3899.pdf"))
plot.pca(dat=dat, features=features.65, sample.class, file.plot=file.path(dir.plot, "PCA_CNA_OV_features65.pdf"))
