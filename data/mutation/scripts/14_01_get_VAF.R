### DEFINE LIBRARIES ---
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/collinsgroup/Raunak/projects/MESO_peritoneal/data/mutation")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.script <- file.path(dir.wrk, "scripts")
dir.annovar <- file.path(dir.wrk, "annovar")

### DEFINE PATH ---
files <- list.files(dir.annovar, pattern="TableAnnovar.hg19_multianno.txt.gz", full.names=FALSE)
file.des <- file.path(dir.analysis, "files_annovar_tumor_normal_pairs.tsv")

### LOAD DESIGN TABLE ---
des <- read.delim(file.des, header=TRUE, stringsAsFactors=FALSE)

### FUNCTION: get.vaf() ---
get.vaf <- function(file.annovar, sampleid){
	# LOAD DATA ---
	dat <- read.delim(file.annovar, header=FALSE, stringsAsFactors=FALSE, skip=1)

	# GET DATA --
	#dat$AO <- as.numeric(unlist(lapply(str_split(unlist(lapply(str_split(dat$V62, ";"), function(x) x[2])), "="), function(x) x[2])))
	#dat$DP <- as.numeric(unlist(lapply(str_split(unlist(lapply(str_split(dat$V62, ";"), function(x) x[3])), "="), function(x) x[2])))

	if(str_detect(dat$V62[1], "FreqT")){
		dat$FreqT <- as.numeric(unlist(lapply(str_split(unlist(lapply(str_split(dat$V62, ";"), function(x) x[38])), "="), function(x) x[2])))
	} else{
		z <- unlist(lapply(str_split(unlist(lapply(str_split(dat$V62, ";"), function(x) x[2])), "="), function(x) x[2]))
		y.AO <- unlist(lapply(lapply(str_split(z,","), function(x) sort(as.numeric(x), decreasing=TRUE)), function(x) x[1]))	
		y.DP <- as.numeric(unlist(lapply(str_split(unlist(lapply(str_split(dat$V62, ";"), function(x) x[3])), "="), function(x) x[2])))
		dat$FreqT <- y.AO/y.DP
	}		

	# TRIM DATA ---
	df <- dat[,c(1:5,7,65)]
	df$V1 <- str_replace(df$V1, "chr", "")
	colnames(df) <- c("Chromosome","Start_Position","End_Position","Ref","Alt","Gene","FreqT")
	df$SampleID <- sampleid

	return(df)
}


### PROCESS DATA ---
list.df <- list()
for(i in 1:nrow(des)){
	sampleid <- des$SampleID[i]
	fileid <- des$Tumor[i]
	file.annovar <- file.path(dir.annovar, fileid)

	# CALL FUNCTION ---
	list.df[[i]] <- get.vaf(file.annovar, sampleid)

	cat("PROCESSED:", sampleid, "\n", sep="\t")
}

### AGGREGATE DATA ---
df <- do.call(rbind.data.frame, list.df)

### WRITE OUTPUT ---
file.output <- file.path(dir.wrk, "main_calls/annovar_all_vaf.tsv")
write.table(df, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

