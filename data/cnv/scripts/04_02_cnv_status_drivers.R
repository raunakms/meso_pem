### LOAD LIBRARIES ---
library("stringr")

### DEFINE PATH ----
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/cnv")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")
dir.task <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/04_driver_stats/output")

### DEFINE FILES ---
file.dat <- file.path(dir.data, "meso_cnv_calls_nexus.tsv")
file.driver <- file.path(dir.task, "meso_driver_genelist.txt")

### Load Drivers ----
genes.driver <- read.table(file.driver, header=F, stringsAsFactors=F)$V1

### Load CNA DATA ----
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
dat <- subset(dat, dat$SampleID != "MESO-20")

### PREPARE DATA ---
ids <- unique(dat$SampleID)
ctr <- 1
list.dat <- list()
for(i in 1:length(ids)){
	dat.temp1 <- subset(dat, dat$SampleID == ids[i])
	event <- unique(dat.temp1$Event)

	for(j in 1:length(event)){
		dat.temp2 <- subset(dat.temp1, dat.temp1$Event == event[j])

		g <- paste(unique(unlist(str_split(dat.temp2$Gene.Symbols, ","))), collapse=":")
		list.dat[[ctr]] <- data.frame(SampleID=ids[i], Event=event[j], Genes=g)
		ctr <- ctr + 1
	}
	cat("PROCESSED", i, "OF", length(ids), "\n", sep="\t")
}

dat2 <- do.call(rbind.data.frame, list.dat)
dat2$SampleID <- as.character(dat2$SampleID)
dat2$Event <- as.character(dat2$Event)
dat2$Genes <- as.character(dat2$Genes)
dat2$Status <- ifelse( (dat2$Event == "High Copy Gain") | (dat2$Event == "CN Gain"), "AMP;", "DEL;")

### RETAIN DRIVER GENES ONLY ---
for(i in 1:nrow(dat2)){
	dat2$Genes[i] <- paste(intersect(str_split(dat2$Genes[i], ":")[[1]], genes.driver), collapse=":")
}

dat2 <- subset(dat2, dat2$Genes != "")


### PREPARE DATA ---
list.df <- list()
for(i in 1:nrow(dat2)){
	genes <- str_split(dat2$Genes[i], ":")[[1]]
	list.df[[i]] <- data.frame(SampleID=dat2$SampleID[i],
								Gene=genes,
								Status=dat2$Status[i])
}
df <- do.call(rbind.data.frame, list.df)
df <- df[!duplicated(df),]

### OUTPUT DATA ---
file.output <- file.path(dir.analysis, "MESO-PeM_drivergene_cnv_status.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)
