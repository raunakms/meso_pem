### DEFINE PATH -------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/02_circos_plot")
dir.data <- file.path(dir.wrk, "data")
dir.output <- file.path(dir.wrk, "data_merge")

dir.samples <- list.dirs(dir.data, full.names=T, recursive=FALSE)

### LOAD DATA ----------------------
list.mut <- list.cna <- list.fush <- list.fusl <- list()
ctr1 <- ctr2 <- ctr3 <- ctr4 <- 1
for(i in 1:length(dir.samples)){
	file.mut <- list.files(path=dir.samples[i], pattern="mutation.tile", full.names=TRUE)
	file.cna <- list.files(path=dir.samples[i], pattern="cna.histo", full.names=TRUE)
	file.fush <- list.files(path=dir.samples[i], pattern="fusion.histo", full.names=TRUE)
	file.fusl <- list.files(path=dir.samples[i], pattern="fusion.links", full.names=TRUE)
	
	if(length(file.mut) != 0){
		list.mut[[ctr1]] <- read.table(file.mut, sep=" ", header=F, stringsAsFactors=F)
		ctr1 <- ctr1 + 1
	}
	
	if(length(file.cna) != 0){
		list.cna[[ctr2]] <- read.table(file.cna, sep=" ", header=F, stringsAsFactors=F)
		ctr2 <- ctr2 + 1
	}
	
	if(length(file.fush) != 0){
		list.fush[[ctr3]] <- read.table(file.fush, sep=" ", header=F, stringsAsFactors=F)
		ctr3 <- ctr3 + 1
	}

	if(length(file.fusl) != 0){
		list.fusl[[ctr4]] <- read.table(file.fusl, sep=" ", header=F, stringsAsFactors=F)
		ctr4 <- ctr4 + 1
	}
}

dat.mut <- do.call(rbind.data.frame, list.mut)
dat.cna <- do.call(rbind.data.frame, list.cna)
dat.fush <- do.call(rbind.data.frame, list.fush)
dat.fusl <- do.call(rbind.data.frame, list.fusl)

dim(dat.mut)
dim(dat.cna)
dim(dat.fush)
dim(dat.fusl)

dat.mut <- dat.mut[!duplicated(dat.mut),]
dat.cna <- dat.cna[!duplicated(dat.cna),]
dat.fush <- dat.fush[!duplicated(dat.fush),]
dat.fusl <- dat.fusl[!duplicated(dat.fusl),]

file.mut <- file.path(dir.output, "MESO_all_mutation.tile")
file.cna <- file.path(dir.output, "MESO_all_cna.histo")
file.fush <- file.path(dir.output, "MESO_all_fusion.histo")
file.fusl <- file.path(dir.output, "MESO_all_fusion.links")

write.table(dat.mut, file.mut, sep=" ", row.names=F, col.names=F, quote=F)
write.table(dat.cna, file.cna, sep=" ", row.names=F, col.names=F, quote=F)
write.table(dat.fush, file.fush, sep=" ", row.names=F, col.names=F, quote=F)
write.table(dat.fusl, file.fusl, sep=" ", row.names=F, col.names=F, quote=F)
