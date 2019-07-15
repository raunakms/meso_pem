### DEFINE LIBRARIES ----
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/proteome")
dir.data <- file.path(dir.wrk, "data_ptm_categories")
dir.analysis <- file.path(dir.wrk, "processed_data_ptm")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")

### DEFINE FILES ---
file.dat <- file.path(dir.data, "step6.txt.gz")
file.des <- file.path(dir.data, "design_table_proteome_ptm.tsv")
file.array <- file.path("/Data/Raunak/softwares/bdvtools/array_process/array_preprocess.R")

### LOAD DESIGN TABLE ---
des <- read.delim(file.des, header=T, stringsAsFactors=F)
des <- des[which(str_detect(des$Column.Names.Tbl, "Abundance") == TRUE),]

### LOAD DATA ---
dat <- read.delim(file.dat, header=TRUE, stringsAsFactors=FALSE)

# COLUMNS ---
# 1:3 	: DESCRIPTION
# 4:30	: CONFIDENCE SCORE
# 34:60	: ABUNDANCE VALUE

### FUNCTION: getData ---
getData <- function(dat, des, analysis){
	# GET DESIGN TABLE IDS ---
	ids.normal <- sort(des$SequencingID[which(des$SampleType == "Normal")], decreasing=FALSE)
	ids.tumor <- sort(des$SequencingID[which(des$SampleType == "Tumor")], decreasing=FALSE)
	ids.celline <- sort(des$SequencingID[which(des$SampleType == "Celline")], decreasing=FALSE)

	# TRIM DATA ---
	if(analysis == "CONFIDENCE"){
		df <- dat[,c(2:3,4:30)]
	}else if(analysis == "ABUNDANCE"){
		df <- dat[,c(2:3,34:60)]
	}

	# REPLACE COLNAMES ---
	colnames(df) <- c("Gene","PTM",des$SequencingID[1:27])
	df <- subset(df, select=c("Gene","PTM",ids.normal, ids.tumor))

	# MODIFY PTM NAME ---
	df$PTM[which(df$PTM == "Acetyl")] <- "Acetylation"
	df$PTM[which(df$PTM == "Methyl")] <- "Methylation"
	df$PTM[which(df$PTM == "Phospho")] <- "Phosphorylation"
	df$PTM[which(df$PTM == "Label:2H(4)+GG")] <- "Ubiquitination"
	df$PTM[which(df$PTM == "Carbamidomethyl")] <- "Carboxyamidomethylation"
	df$PTM[which(df$PTM == "Oxidation")] <- "Oxidation"
	df$PTM[which(df$PTM == "Gln->pyro-Glu")] <- "PyroglutamicAcid"

	# ORDER DATA ---
	df <- df[order(df$Gene, df$PTM),]

	return(df)
}

### FUNCTION: getPTMAbundance ---
getDataByType <- function(dat, ptm){
	# GET PTM CATEGORIES ---
	df <- subset(dat, dat$PTM == ptm)

	# CREATE MATRIX ---
	mat <- as.matrix(df[,3:ncol(df)])
	rownames(mat) <- df$Gene

	return(mat)
}

### FUNCTION: generateData ---
generateData <- function(dir.output, dat, ptm, analysis){
	# GET DATA ---
	mat <- getDataByType(dat, ptm)

	# WRITE OUPUT ---
	file.output <- file.path(dir.output, paste("protein_", analysis, "_", ptm, ".tsv", sep=""))
	write.table(mat, file.output, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
}

### GET PTM DATA ---
dat.ptm <- getData(dat, des, analysis="ABUNDANCE")
generateData(dir.output=dir.analysis, dat=dat.ptm, ptm="Acetylation", analysis="abundance")
generateData(dir.output=dir.analysis, dat=dat.ptm, ptm="Methylation", analysis="abundance")
generateData(dir.output=dir.analysis, dat=dat.ptm, ptm="Phosphorylation", analysis="abundance")
generateData(dir.output=dir.analysis, dat=dat.ptm, ptm="Ubiquitination", analysis="abundance")
generateData(dir.output=dir.analysis, dat=dat.ptm, ptm="Carboxyamidomethylation", analysis="abundance")
generateData(dir.output=dir.analysis, dat=dat.ptm, ptm="Oxidation", analysis="abundance")
generateData(dir.output=dir.analysis, dat=dat.ptm, ptm="PyroglutamicAcid", analysis="abundance")

### GET PTM CONFIDENCE DATA ---
dat.conf <- getData(dat, des, analysis="CONFIDENCE")
generateData(dir.output=dir.analysis, dat=dat.conf, ptm="Acetylation", analysis="confidence")
generateData(dir.output=dir.analysis, dat=dat.conf, ptm="Methylation", analysis="confidence")
generateData(dir.output=dir.analysis, dat=dat.conf, ptm="Phosphorylation", analysis="confidence")
generateData(dir.output=dir.analysis, dat=dat.conf, ptm="Ubiquitination", analysis="confidence")
generateData(dir.output=dir.analysis, dat=dat.conf, ptm="Carboxyamidomethylation", analysis="confidence")
generateData(dir.output=dir.analysis, dat=dat.conf, ptm="Oxidation", analysis="confidence")
generateData(dir.output=dir.analysis, dat=dat.conf, ptm="PyroglutamicAcid", analysis="confidence")

