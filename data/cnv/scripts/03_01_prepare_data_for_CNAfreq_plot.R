### LOAD LIBRARIES ---
library("stringr")

### DEFINE PATH ----
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/cnv")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")

### DEFINE FILES ---
file.dat <- file.path(dir.data, "meso_cnv_nexus_table_2016_07_06.txt")
file.des <- file.path(dir.data, "samplelist_table.tsv")

### LOAD DESIGN TABLE ----
des <- read.delim(file.des, header=T, stringsAsFactors=F)

### Load CNA DATA ----
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)

### ADD SAMPLEIDS ---
dat$SampleID <- ""
for(i in 1:nrow(des)){
	dat$SampleID[which(dat$Sample == des$SampleName[i])] <- des$SampleID[i]
}
dat <- dat[order(dat$SampleID, decreasing=F),]

### ADD ATTRIBUTES ---
dat$Chromosome.Region <- str_replace_all(dat$Chromosome.Region, ",", "")
dat$Chr <- str_replace(unlist(lapply(str_split(dat$Chromosome.Region, ":"), function(x) x[1])), "chr", "")
dat$Start <- unlist(lapply(str_split(unlist(lapply(str_split(dat$Chromosome.Region, ":"), function(x) x[2])), "-"), function(x) as.numeric(x[1])))
dat$End <- unlist(lapply(str_split(unlist(lapply(str_split(dat$Chromosome.Region, ":"), function(x) x[2])), "-"), function(x) as.numeric(x[2])))
dat$Cytoband <- unlist(lapply(str_split(dat$Cytoband, " - "), function(x) paste(x, collapse="-")))
dat$arm <- str_extract(dat$Cytoband, "\\w")

### FILTER OUT SAMPLE ---
dat <- subset(dat, dat$SampleID != "MESO-20")

### TRIM DATA ---
df <- subset(dat, select=c("SampleID", "Chr", "arm", "Start", "End", "Probes", "Probe.Median"))
colnames(df) <- c("sampleID","chrom","arm","start.pos","end.pos","n.probes","mean")

### Chromosome REPLACE WITH NUMERIC VALUES ---
#df$chrom[which(df$chrom == "X")] <- "22"
#df$chrom[which(df$chrom == "Y")] <- "23"
#df$chrom <- as.numeric(df$chrom)

### ORDER DATA ---
df <- df[order(df$sampleID, df$chrom, df$start.pos),]

### OUTPUT DATA ---
file.output <- file.path(dir.data, "MESO-PeM_cnv_seg_calls_for_chromFreq_plot.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)



### DATA FOR CLUSTERED HEATMAP ----
file.dat <- file.path(dir.data, "MESO-PeM_cnv_seg_calls_for_chromFreq_plot.tsv")
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)

file.des <- file.path(dir.data, "sample_table_forcluster.tsv")
des <- read.delim(file.des, header=T, stringsAsFactors=F)
des$ClusterSample <- paste("A", sprintf("%02d", c(19:1)), sep="")

y <- rep("",length=nrow(dat))
for(i in 1:nrow(des)){
	index <- which(dat$sampleID == des$SampleID[i])
	y[index] <- des$ClusterSample[i]
}

dat$sampleID <- y

### OUTPUT DATA ---
file.output <- file.path(dir.data, "MESO-PeM_cnv_seg_calls_forcluster.tsv")
write.table(dat, file.output, sep="\t", row.names=F, col.names=T, quote=F)
