### Load LIBRARIES ----
library("stringr")

### DEFINE PATH ----
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/04_driver_stats")
dir.data <- file.path(dir.wrk, "data")
dir.output <- file.path(dir.wrk, "output")
dir.plot <- file.path(dir.wrk, "plot")

### DEFINE FILES -----
file.driver <- file.path(dir.output, "meso_driver_genes.tsv")
file.alt <- file.path(dir.data, "meso_alterations_samplewise_status_all.tsv")

### GET DRIVER DATA ---
dat.driver <- read.delim(file.driver, header=T, stringsAsFactors=F)
ids <- unique(dat.driver$SampleID)
genes.driver <- unique(unlist(str_split(dat.driver$Genes, ":")))

### GET ALTERATION DATA ---
dat.alt <- read.delim(file.alt, header=T, stringsAsFactors=F)
dat.alt0 <- subset(dat.alt, !(dat.alt$SampleID %in% ids))
dat.alt1 <- subset(dat.alt, dat.alt$SampleID %in% ids)


### GET ALTERATION FOR MATCH ---
list.alt1 <- list()
for(i in 1:nrow(dat.driver)){
	sampleid <- dat.driver$SampleID[i]
	genes <- str_split(dat.driver$Genes[i], ":")[[1]]

	dat.temp <- subset(dat.alt1, dat.alt1$SampleID == sampleid)
	dat.temp <- subset(dat.temp, dat.temp$Gene %in% genes)
	list.alt1[[i]] <- dat.temp
}

df1 <- do.call(rbind.data.frame, list.alt1)

### GET ALTERATION FOR NON-MATCH ---
df2 <- subset(dat.alt0, dat.alt0$Gene %in% genes.driver)

### MERGE DATA ---
df <- rbind(df1, df2)

#df$Status <- str_replace(df$Status, "INDEL;DEL;", "DEL;")
df$Status <- str_replace(df$Status, "SNP;", "MUT;")
df$Status <- str_replace(df$Status, "INDEL;", "IND;")

### WRITE OUTPUT ---
file.output <- file.path(dir.output, "meso_driver_oncoprint.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)


### LOAA DATA ---
file.df <- file.path(dir.data, "MESO-PeM_drivergene_alt_status.tsv")
df <- read.delim(file.df, header=T, stringsAsFactors=F)

df$Status <- str_replace(df$Status, "SNP;", "MUT;")
df$Status <- str_replace(df$Status, "INDEL;", "IND;")

### PLOT ONCOPRINT -----------------------------------------------------------------------
file.script <- file.path("/Data/Raunak/softwares/bdvtools/oncoprint/get_oncoprint_plot.R")
source(file.script)

file.plot <- file.path(dir.plot, "meso_driver_oncoprint_new.pdf")  
pdf(file.plot, width=2, height=4)
	oncoprint(df, orderGenes="memo.sort", 
		file.genes=NA, file.samples=NA, circularize=FALSE, 
		cex.text=4, cex.axis=2, cex.title=2, 
		plot.title.name="")
dev.off()


