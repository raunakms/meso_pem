### DEFINE LIBRARIES ---
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/collinsgroup/Raunak/projects/MESO_peritoneal/data/mutation")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.script <- file.path(dir.wrk, "scripts")
dir.annovar <- file.path( "/home/collins/IonProton/Mesothelioma/annovar")

### DEFINE PATH ---
files <- list.files(dir.annovar, pattern="TableAnnovar.hg19_multianno.txt.gz", full.names=FALSE)
file.des <- file.path(dir.analysis, "files_annovar_tumor_normal_pairs.tsv")

### LOAD DESIGN TABLE ---
des <- read.delim(file.des, header=TRUE, stringsAsFactors=FALSE)

### ITEMS COLNAMES ---
items <- c("Chr","Start","End","Ref","Alt",
			"Func.refGene","Gene.refGene","GeneDetail.refGene",
			"ExonicFunc.refGene","AAChange.refGene","cytoBand",
			"ExAC_ALL","ExAC_AFR","ExAC_AMR","ExAC_EAS",
			"ExAC_FIN","ExAC_NFE","ExAC_OTH","ExAC_SAS",
			"avsnp147","SIFT_score","SIFT_pred",
			"Polyphen2_HDIV_score","Polyphen2_HDIV_pred",
			"Polyphen2_HVAR_score","Polyphen2_HVAR_pred",
			"LRT_score","LRT_pred",
			"MutationTaster_score","MutationTaster_pred",
			"MutationAssessor_score","MutationAssessor_pred",
			"FATHMM_score","FATHMM_pred","PROVEAN_score","PROVEAN_pred",
			"VEST3_score","CADD_raw","CADD_phred","DANN_score",
			"fathmm-MKL_coding_score","fathmm-MKL_coding_pred",
			"MetaSVM_score","MetaSVM_pred","MetaLR_score","MetaLR_pred",
			"integrated_fitCons_score","integrated_confidence_value",
			"GERP++_RS","phyloP7way_vertebrate","phyloP20way_mammalian",
			"phastCons7way_vertebrate","phastCons20way_mammalian",
			"SiPhy_29way_logOdds","Otherinfo",
			"V1","V2","V3","V4","V5","V6","V7","V8","V9")


### NORMALIZE VARIANTS BY NORMAL ---
for(i in 1:nrow(des)){
	if(des$Normal != ""){
		# DEFINE FILE ---
		file.normal <- file.path(dir.annovar, des$Normal[i])
		file.tumor <- file.path(dir.annovar, des$Tumor[i])

		# READ DATA --
		dat.normal <- read.delim(file.normal, header=F, stringsAsFactors=F, skip=1)
		dat.tumor <- read.delim(file.tumor, header=F, stringsAsFactors=F, skip=1)

		# DEFINE PRIMARY KEY ---
		dat.normal$Key <- apply(dat.normal, 1, function(x) paste(as.character(x[1]), as.numeric(x[2]), as.numeric(x[3]), as.character(x[4]), as.character(x[5]), sep=":"))
		dat.tumor$Key <- apply(dat.tumor, 1, function(x) paste(as.character(x[1]), as.numeric(x[2]), as.numeric(x[3]), as.character(x[4]), as.character(x[5]), sep=":"))

		# MATCH VARIANT IN TUMOR WITH THAT OF NORMAL ---
		dat.tumor$Status <- 1
		dat.tumor$Status[which(dat.tumor$Key %in% dat.normal$Key)] <- 0

		# RETAIN VARIAN WITH STATUS == 1 ---
		dat.tumor <- subset(dat.tumor, dat.tumor$Status == 1)
		dat <- dat.tumor[,-which(colnames(dat.tumor) == "Status")]

	} else{
		# DEFINE FILE ---
		file.tumor <- file.path(dir.annovar, des$Tumor[i])

		# READ DATA --
		dat.tumor <- read.delim(file.tumor, header=F, stringsAsFactors=F, skip=1)
		colnames(dat.tumor) <- items

	}
}
