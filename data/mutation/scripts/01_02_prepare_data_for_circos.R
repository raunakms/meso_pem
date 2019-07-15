### SetDirectories ----------------------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/mutation")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")
dir.output <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/02_circos_plot/data")

### Load Libraries ----------------------------------------------------------
library("stringr")

### Define Files ------------------------------------------------------------
file.dat <- file.path(dir.analysis, "somatic_mutation_non_silent_samplewise.tsv")
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
dat$Chr <- paste("hs", dat$Chr, sep="")

ids <- unique(dat$SampleID)


for(i in 1:length(ids)){
	dir.id <- file.path(dir.output, ids[i])
	dir.create(dir.id, showWarnings=FALSE)

	df <- subset(dat, dat$SampleID == ids[i])[,c(2:4,9)]
	
	file.output <- file.path(dir.id, paste(ids[i], "mutation.tile", sep="_"))
	write.table(df, file.output, sep=" ", row.names=F, col.names=F, quote=F)
}
