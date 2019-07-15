### DEFINE LIBRARIES ---
library("maftools")
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/mutation")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.script <- file.path(dir.wrk, "scripts")
dir.maf <- file.path(dir.wrk, "maf_ver2")

### DEFINE FILES ---
file.dat <- file.path(dir.maf, "meso_pem_benign_tumor_noncds.maf.gz")

### LOAD DATA ---
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)

### CHROM ---
dat$Chromosome <- str_replace_all(dat$Chromosome, "chr", "")


### GENERATE KEY ---
dat$key <- apply(dat, 1, function(x) 
			paste(as.character(x[1]), as.character(x[5]),
				as.numeric(x[6]), as.numeric(x[7]) ,
				as.character(x[9]), as.character(x[11]),
				as.character(x[12]), as.character(x[13]),
				sep=":"))


### GROUP SAMPLEIDS ---
ids <- unique(dat$Tumor_Sample_Barcode)
ids.normal <- ids[which(str_detect(ids, "N") == TRUE)]
ids.tumor <- ids[which(str_detect(ids, "T") == TRUE)]

### SEPARATE DATA BY SAMPLE TYPE ---
dat.normal <- subset(dat, dat$Tumor_Sample_Barcode %in% ids.normal)
dat.tumor <- subset(dat, dat$Tumor_Sample_Barcode %in% ids.tumor)

### FLAG POSSIBLE GERMLINE MUTATION ---
dat.tumor$Flag <- 0
for(i in 1:nrow(dat.tumor)){
	index <- which(dat.normal$key == dat.tumor$key[i])
	if(length(index) != 0){
		dat.tumor$Flag[i] <- 1
		cat("FLAG","\n", sep=" ")
	}
	#cat("PROCESSED:", i, "OF", nrow(dat.tumor), "\n", sep="\t")
}

### REMOVE POSSIBLE GERMLINE MUTATION ---
dat.tumor <- subset(dat.tumor, dat.tumor$Flag == 0)
dat.tumor <- dat.tumor[,-c(75,76)]
dat.tumor$Tumor_Sample_Barcode <- str_replace(dat.tumor$Tumor_Sample_Barcode, "T", "")
dat.tumor <- dat.tumor[,1:21]
dat.tumor <- dat.tumor[!duplicated(dat.tumor),]

### WRITE OUTPUT ---
file.output <- file.path(dir.maf, "meso_pem_tumor_noncds.maf")
write.table(dat.tumor, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#cmd <- paste("gzip", file.output, sep=" ")
#system(cmd)



### DEFINE FILE ---
file.maf <- file.path(dir.maf, "meso_pem_benign_tumor_noncds_tumor.maf")

# LOAD DATA ---
dat.maf <- read.maf(file.maf, removeSilent = FALSE, useAll = FALSE,
						gisticAllLesionsFile = NULL, gisticAmpGenesFile = NULL,
						gisticDelGenesFile = NULL, cnTable = NULL,
						removeDuplicatedVariants = FALSE, isTCGA = FALSE, verbose = TRUE)


### DRAW ONCO PLOT ---
file.plot <- file.path(dir.maf, "oncoplot_mutation_noncds.pdf")
pdf(file.plot, width=15, height=10)
oncoplot(maf = dat.maf, writeMatrix = FALSE, top = 50, genes = NULL,
	drawRowBar = TRUE, drawColBar = TRUE, showTumorSampleBarcodes = TRUE,
	annotation = NULL, annotationColor = NULL, genesToIgnore = NULL,
	removeNonMutated = TRUE, colors = c("red","green","blue"), fontSize = 10,
	sortByMutation = FALSE, sortByAnnotation = FALSE)

dev.off()
