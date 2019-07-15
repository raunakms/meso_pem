### DEFINE LIBRARIES --
library("stringr")
#library("ggplot2")
#library("gridExtra")
library("maftools")


### DEFINe PATH ---
dir.wrk <- file.path("//jbrcsrv009/CollinsGroup/Raunak/projects/MESO_peritoneal/data/mutation")
#dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/mutation")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")

### DEFINE FILES ---
file.maf <- file.path(dir.data, "VPCLAGA.MESO_PeM.nonsilent_mutation.maf")

# LOAD DATA ---
#dat.maf <- read.delim(file.maf, header=TRUE, stringsAsFactors=FALSE)
dat.maf <- read.maf(file.maf, removeSilent = TRUE, useAll = TRUE,
						gisticAllLesionsFile = NULL, gisticAmpGenesFile = NULL,
						gisticDelGenesFile = NULL, cnTable = NULL,
						removeDuplicatedVariants = TRUE, isTCGA = FALSE, verbose = TRUE)

### MAF SUMMARY PLOT ---
file.plot <- file.path(dir.plot, "mutation_maf_summary")
plotmafSummary(maf = dat.maf, file = file.plot, rmOutlier = TRUE, dashboard = TRUE,
	width = 6, height = 5, addStat = 'median', showBarcodes = TRUE,
	fs = 8, textSize = 2, color = NULL, statFontSize = 3,
	titvColor = NULL, top = 15)
dev.off()


plotmafSummary(maf = dat.maf, file = NULL, rmOutlier = TRUE, dashboard = FALSE,
	addStat = 'median', showBarcodes = TRUE,
	fs = 5, textSize = 6, color = NULL, statFontSize = 4,
	titvColor = NULL, top = 15)


### DRAW ONCO PLOT ---
oncoplot(maf = dat.maf, writeMatrix = FALSE, top = 20, genes = NULL,
	drawRowBar = TRUE, drawColBar = TRUE, showTumorSampleBarcodes = TRUE,
	annotation = NULL, annotationColor = NULL, genesToIgnore = NULL,
	removeNonMutated = TRUE, colors = NULL, fontSize = 10,
	sortByMutation = FALSE, sortByAnnotation = FALSE)

file.plot <- file.path(dir.plot, "mutation_maf_oncoplot.pdf")
dev.off()


### LOLLIPOP POLT ---------

lpop.bap1 <- lollipopPlot(maf=dat.maf, gene = "BAP1", AACol = "AAChange", labelPos = 'all',
					showMutationRate = FALSE, fn = NULL, showDomainLabel = FALSE,
					cBioPortal = FALSE, refSeqID = NULL, proteinID = NULL, repel = FALSE,
					collapsePosLabel = FALSE, legendTxtSize = 10, labPosSize = 2,
					labPosAngle = 0, domainLabelSize = 2.5, printCount = FALSE,
					colors = NULL, domainColors = NULL, labelOnlyUniqueDoamins = TRUE,
					defaultYaxis = TRUE)











### Classifies SNPs into transitions and transversions ---
# LOAD DATA ---
file.maf <- file.path(dir.data, "VPCLAGA.MESO_PeM.nonsilent_mutation_for_titv.maf")
dat.maf <- read.maf(file.maf, removeSilent = TRUE, useAll = TRUE,
						gisticAllLesionsFile = NULL, gisticAmpGenesFile = NULL,
						gisticDelGenesFile = NULL, cnTable = NULL,
						removeDuplicatedVariants = TRUE, isTCGA = FALSE, verbose = TRUE)


titv(maf = dat.maf, useSyn = FALSE, plot = TRUE, file = NULL)
y <- titv(maf = dat.maf, useSyn = FALSE, plot = TRUE, file = NULL)

plotTiTv(res = y, plotType = "both", file = NULL, width = 4, height = 4, 
			color = NULL, showBarcodes = TRUE, textSize = 2)

d1 <- as.data.frame(y[[1]])
d2 <- as.data.frame(y[[2]])
d3 <- as.data.frame(y[[3]])

write.table(d1, file.path(dir.analysis, "titv_fraction_crontribution_MESO_PeM.tsv"), sep="\t", row.names=F, col.names=T, quote=F)
write.table(d2, file.path(dir.analysis, "titv_raw_counts_MESO_PeM.tsv"), sep="\t", row.names=F, col.names=T, quote=F)
write.table(d3, file.path(dir.analysis, "titv_fractions_MESO_PeM.tsv"), sep="\t", row.names=F, col.names=T, quote=F)



########
dir.maf <- file.path(dir.wrk, "maf_ver2")
file.dat <- file.path(dir.maf, "meso_pem_benign_tumor_cds.maf.gz")
dat <- read.delim(file.dat, header=TRUE, stringsAsFactors=FALSE)
del.index <- which(str_detect(dat$Tumor_Sample_Barcode, "N") == TRUE)
dat <- dat[-del.index,]
write.table(dat, file.path(dir.maf, "meso_pem_tumor_cds.maf"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

file.maf <- file.path(dir.maf, "meso_pem_tumor_cds.maf")
dat.maf <- read.maf(file.maf, removeSilent = TRUE, useAll = TRUE,
						gisticAllLesionsFile = NULL, gisticAmpGenesFile = NULL,
						gisticDelGenesFile = NULL, cnTable = NULL,
						removeDuplicatedVariants = TRUE, isTCGA = FALSE, verbose = TRUE)

file.plot <- file.path(dir.maf, "TiTv.pdf")
pdf(file.plot, width=6, height=6)
y <- titv(maf = dat.maf, useSyn = FALSE, plot = TRUE, file = NULL)
plotTiTv(res = y, plotType = "both", file = NULL, width = 6, height = 6, 
			color = NULL, showBarcodes = TRUE, textSize = 6)
dev.off()



which(str_detect(dat.maf@data$Tumor_Sample_Barcode, "N") == TRUE)
