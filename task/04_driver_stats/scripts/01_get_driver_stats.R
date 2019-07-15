### Load LIBRARIES ----
library("stringr")

### DEFINE PATH ----
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/04_driver_stats")
dir.data <- file.path(dir.wrk, "data")
dir.output <- file.path(dir.wrk, "output")
dir.plot <- file.path(dir.wrk, "plot")

### DEFINE FILES -----
file.driver <- file.path(dir.data, "mesoCombo_patientdrivers.txt")

### GET DRIVER DATA ---
dat.driver <- read.delim(file.driver, header=T, stringsAsFactors=F)
dat.driver <- dat.driver[which(str_detect(dat.driver$patient, "M") == TRUE),]

### Create BiPartite Graph ----
gamma <- unique(dat.driver$gamma)

dat.stat <- data.frame(gamma=gamma)
dat.stat$genes <- ""

for(i in 1:length(gamma)){
	dat.stat$genes[i] <- paste(unique(unlist(str_split(dat.driver$driver.genes[which(dat.driver$gamma == gamma[i])], ","))), collapse=":")
}

### CREATE MATRIX ---
gamma <- unique(dat.driver$gamma)
genes <- unique(unlist(str_split(dat.stat$genes, ":")))
bip <- matrix(0, nrow=length(gamma), ncol=length(genes), dimnames=list(gamma, genes))

for(i in 1:nrow(dat.stat)){
	g <- str_split(dat.stat$genes[i], ":")[[1]]
	index <- which(genes %in% g)
	bip[i, index] <- 1
}

bip <- bip[,match(names(sort(colSums(bip), decreasing=TRUE)), colnames(bip))]

write.table(t(bip), file.path(dir.output, "driver_gamma_bip.tsv"), sep="\t", row.names=T, col.names=NA, quote=F)



#### GET UNION OF DRIVERS PER PATIENTS ACROSS ALL GAMMAS ----
ids <- unique(dat.driver$patient)

dat <- data.frame(SampleID=ids)
dat$Genes <- ""

for(i in 1:length(ids)){
	dat.temp <- subset(dat.driver, dat.driver$patient == ids[i])
	dat$Genes[i] <- paste(unique(unlist(str_split(dat.temp$driver.genes, ","))), collapse=":")
}

dat$SampleID <- str_replace(ids, "M", "MESO-")

### WRITE OUTPUT ---
file.output <- file.path(dir.output, "meso_driver_genes.tsv")
write.table(dat, file.output, sep="\t", row.names=F, col.names=T, quote=F)


### GET ALL DRIVERS ---
genes.driver <- unique(unlist(str_split(dat$Genes, ":")))
write.table(genes.driver, file.path(dir.output, "meso_driver_genelist.txt"), row.names=F, col.names=F, quote=F)










