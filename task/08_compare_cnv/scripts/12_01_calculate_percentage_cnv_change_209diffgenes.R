### DEFINE LIBRARIES ---
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/08_compare_cnv")
dir.data <- file.path(dir.wrk, "data")
dir.output <- file.path(dir.wrk, "output")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")

### DEFINE FILES ---
file.dat1 <- file.path(dir.data, "cnv_status_209diffgenes_meso_pem.tsv")
file.dat2 <- file.path(dir.data, "cnv_status_209diffgenes_meso_pm.tsv")

### LOAD DATA ---
dat1 <- read.delim(file.dat1, header=TRUE, stringsAsFactors=FALSE, row.names=1)
dat2 <- read.delim(file.dat2, header=TRUE, stringsAsFactors=FALSE, row.names=1)

colnames(dat1) <- str_replace_all(colnames(dat1), "[.]", "-")
colnames(dat2) <- str_replace_all(colnames(dat2), "[.]", "-")

#dat1[is.na(dat1)] <- ""
#dat2[is.na(dat2)] <- ""


### FUNCTION: getCNVchange ----
getCNVchange <- function(dat){
	n <- ncol(dat)
	df <- data.frame(Gene=rownames(dat))
	df$Gene <- as.character(df$Gene)
	df$percentage.loss <- df$percentage.gain <- 0
	for(i in 1:nrow(dat)){
		namp <- length(unique(which((dat[i,] == "CN Gain") | (dat[i,] == "High Copy Gain") )))
		ndel <- length(unique(which((dat[i,] == "CN Loss") | (dat[i,] == "Homozygous Copy Loss") )))

		df$percentage.gain[i] <- (namp / n) * 100
		df$percentage.loss[i] <- (ndel / n) * 100
	}
	return(df)
}

### CALL FUNCTION ----
df.pem <- getCNVchange(dat1)
df.pm <- getCNVchange(dat2)

### COMBINE DATA ---
df <- data.frame(Gene=rownames(dat1))
df$percentage.gain.pem <- df.pem$percentage.gain
df$percentage.loss.pem <- df.pem$percentage.loss
df$percentage.gain.pm <- df.pm$percentage.gain
df$percentage.loss.pm <- df.pm$percentage.loss

### WRITE OUTPUT ---
file.output <- file.path(dir.output, "cnv_percentage_209diffgenes_meso_pem_pm.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)
