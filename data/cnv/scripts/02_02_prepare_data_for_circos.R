### SetDirectories ----------------------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/cnv")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")
dir.output <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/02_circos_plot/data")

### Load Libraries ----------------------------------------------------------
library("stringr")

### Define Files ------------------------------------------------------------
file.dat <- file.path(dir.data, "meso_cnv_calls_nexus.tsv")

### Load CNA Data --------------------------------------------------------
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
dat$Chr <- str_replace(dat$Chr, "chr", "hs")
dat <- dat[,c(1,2,3,4,9)]

dat <- subset(dat, dat$Sample != "MESO-20")
ids <- unique(dat$SampleID)

for(i in 1:length(ids)){
	dir.id <- file.path(dir.output, ids[i])
	df <- subset(dat, dat$SampleID == ids[i])[,-1]
	
	file.output <- file.path(dir.id, paste(ids[i], "cna.histo", sep="_"))
	write.table(df, file.output, sep=" ", row.names=F, col.names=F, quote=F)
}
