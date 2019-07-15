### DEFINE LIBRARIES ---
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/mutation")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.script <- file.path(dir.wrk, "scripts")
dir.annovar <- file.path( "/home/collins/IonProton/Mesothelioma/annovar")

### DEFINE PATH ---
files <- list.files(dir.annovar, pattern="TableAnnovar.hg19_multianno.txt.gz", full.names=FALSE)
ids <- unlist(lapply(str_split(files, "[.]"), function(x) x[1]))

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

### 
list.df <- list()
ctr <- 1
for(i in 1:length(files)){
	# LOAD DATA --
	file.dat <- file.path(dir.annovar, files[i])
	dat <- read.delim(file.dat, header=F, stringsAsFactors=F, skip=1)
	colnames(dat) <- items

	# GET BAP1 mutation ---
	index <- which(dat$Gene.refGene == "BAP1")

	if(length(index) != 0){
		d <- dat[index,]
		d$SampleID <- ids[i]
		list.df[[ctr]] <- d
		ctr <- ctr + 1
	}

	cat(ids[i], "BAP1 :", length(index), "\n", sep="\t")
}

df <- do.call(rbind.data.frame, list.df)

write.table(df, file.path(dir.analysis, "meso_germline_TableAnnovar.hg19_multianno.tsv"), sep="\t", row.names=F, col.names=T, quote=F)
