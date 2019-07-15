#### LOAD LIBRARIES -------------------
library("stringr")

#### DEFINE PATH ----------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.proteome <- file.path(dir.wrk, "data/proteome/processed_data_ptm")
dir.des <- file.path(dir.wrk, "data/proteome/data_ptm_categories")
dir.task <- file.path(dir.wrk, "task/06_proteomics")
dir.data <- file.path(dir.task, "data")
dir.output <- file.path(dir.task, "output/ptm_by_categories_wilcoxn_test")
dir.plot <- file.path(dir.task, "plot")

### LOAD FUNCTION ---
file.script <- file.path("/Data/Raunak/softwares/bdvtools/array_process/diff_expr.R")
source(file.script)

### DEFINE FILES ---
files <- list.files(dir.proteome, pattern="protein_abundance", full.names=FALSE)
file.des <- file.path(dir.des, "design_table_proteome_ptm.tsv")

# LOAD DESIGN TABLE ---
des <- read.delim(file.des, header=TRUE, stringsAsFactors=FALSE)
des <- des[which(str_detect(des$Column.Names.Tbl, "Abundance") == TRUE),]

ids.normal <- des$SequencingID[which(des$SampleType == "Normal")]
ids.tumor <- des$SequencingID[which(des$SampleType == "Tumor")]
ids.celline <- des$SequencingID[which(des$SampleType == "Celline")]
ids.bap1del <- des$SequencingID[which((des$SampleType == "Tumor") & (des$PeM.Subtype == "BAP1-DEL"))]
ids.bap1int <- des$SequencingID[which((des$SampleType == "Tumor") & (des$PeM.Subtype == "BAP1-INTACT"))]

### FUNCTION: ---
testDifference <- function(file.dat, analysis, group1, group2){
	# LOAD DATA ---
	dat <- read.delim(file.dat, header=TRUE, stringsAsFactors=FALSE, row.names=1)
	colnames(dat) <- str_replace_all(colnames(dat), "[.]", "-")

	# SELECT DATA ---
	dat <- subset(dat, select=c(group1, group2))

	# PERFORM WILXON RANK SUM TEST ---
	df <- get.wilcox.rank.test(dat, class1=group1, class2=group2)

	# COMPILE DATA ---
	df <- df[order(df$pvalue, decreasing=FALSE),]
	#df <- subset(df, df$pvalue <= 0.1)

	return(df)
}


### FUNCTION: getAnalysis ---
getAnalysis <- function(files, type, group1, group2){
	for(i in 1:length(files)){

		file <- files[i]
		file.dat <- file.path(dir.proteome, file)
		analysis <- str_replace(str_split(file, "_")[[1]][3], ".tsv", "")

		cat("START:", analysis, type, "\n", sep="\t")

		# PERFORM TEST ---
		df <- testDifference(file.dat, analysis, group1, group2)

		# WRITE OUTPUT ---
		file.output <- file.path(dir.output, paste("ptm_abundance_", analysis, "_", type, ".tsv", sep=""))
		write.table(df, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

		cat("PROCESSED:", analysis, type, "\n", sep="\t")
	}
}


### CALL FUNCTION ---
getAnalysis(files, type="normal-tumor", group1=ids.normal, group2=ids.tumor)
getAnalysis(files, type="pem-subtypes", group1=ids.bap1int, group2=ids.bap1del)
