### SetDirectories ----------------------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/mutation")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")

### Define Files ------------------------------------------------------------
file.dat <- file.path(dir.data, "All_Combined.variant.vcf2annovar.TableAnnovar.hg19_multianno.txt.AAA_Splicing.S2.non-dbsnp.Germline.tsv")

### Load Libraries ----------------------------------------------------------
library("stringr")

### Load Mutation Data --------------------------------------------------------
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
colnames(dat)[136:171] <- str_replace_all(colnames(dat)[136:171], "[.]", "-")
#dat <- dat[,c(1:10, 136:171)]

#### GT:AF:AO:DP:FAO:FDP:FRO:FSAF:FSAR:FSRF:FSRR:GQ:RO:SAF:SAR:SRF:SRR
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency based on Flow Evaluator observation counts">
##INFO=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observations">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
##INFO=<ID=FAO,Number=A,Type=Integer,Description="Flow Evaluator Alternate allele observations">
##INFO=<ID=FDP,Number=1,Type=Integer,Description="Flow Evaluator read depth at the locus">
##INFO=<ID=FRO,Number=1,Type=Integer,Description="Flow Evaluator Reference allele observations">
##INFO=<ID=FSAF,Number=A,Type=Integer,Description="Flow Evaluator Alternate allele observations on the forward strand">
##INFO=<ID=FSAR,Number=A,Type=Integer,Description="Flow Evaluator Alternate allele observations on the reverse strand">
##INFO=<ID=FSRF,Number=1,Type=Integer,Description="Flow Evaluator Reference observations on the forward strand">
##INFO=<ID=FSRR,Number=1,Type=Integer,Description="Flow Evaluator Reference observations on the reverse strand">
##INFO=<ID=RO,Number=1,Type=Integer,Description="Reference allele observations">
##INFO=<ID=SAF,Number=A,Type=Integer,Description="Alternate allele observations on the forward strand">
##INFO=<ID=SAR,Number=A,Type=Integer,Description="Alternate allele observations on the reverse strand">
##INFO=<ID=SRF,Number=1,Type=Integer,Description="Number of reference observations on the forward strand">
##INFO=<ID=SRR,Number=1,Type=Integer,Description="Number of reference observations on the reverse strand">

### Parse Files ---------------------------------------------------------------
dat.list <- list()
ctr <- 1
for(i in 136:171){
	id <- colnames(dat)[i]
	d1 <- dat[which(dat[,i] != "."),1:26]
	d2 <- dat[which(dat[,i] != "."),i]

	dat.list[[ctr]] <- cbind(d1, AF_AO=d2, SampleID=id)
	ctr <- ctr + 1
}
df <- do.call(rbind.data.frame, dat.list)
df$SampleID <- as.character(df$SampleID)

file.output <- file.path(dir.analysis, "germline_mutation_samplewise.vcf")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)

### Group Data ------------------------------------------------------------------
ids <- unique(df$SampleID)
ids.normal <- ids[which(str_detect(ids, "-N") == TRUE)]
ids.tumor <- ids[which(str_detect(ids, "-N") == FALSE)]
ids.target <- c("MESO-18A","MESO-18E","MESO-19")

### GET POOLED NORMAL ---
df.normal <- subset(df, df$SampleID %in% ids.normal)
df.normal$SampleID <- "MESO-POOLED-NORMAL"

### GET TUMOR SAMPLES ---
df.tumor <- subset(df, df$SampleID %in% ids.tumor)
df.tumor <- subset(df.tumor, !(df.tumor$SampleID %in% ids.target))

### GET TARGET SAMPLES ---
df.target <- subset(df, df$SampleID %in% ids.target)

### FILTER BY DBSNP138 ---
df.target <- df.target[which(df.target$snp138 == "."),]

### FUNCTION : GET KEY ---
get.key <- function(dat){
	key <- apply(dat, 1, function(x) paste( as.character(x[1]), as.numeric(x[2]), as.numeric(x[3]), 
											as.character(x[4]), as.character(x[5]), sep=":"))	
	return(key)
}

### CALL FUNCTION ---
df.normal$key <- get.key(df.normal)
df.target$key <- get.key(df.target)

### MATCH GERMLINE VARIANTS ---
common.key <- intersect(df.target$key, df.normal$key)

#yn <- subset(df.normal, df.normal$key %in% common.key)
#yt <- subset(df.target, df.target$key %in% common.key)
#yn <- yn[match(common.key, yn$key),]
#yt <- yt[match(common.key, yt$key),]
#index.germline.mut <- which(df.target$key %in% df.normal$key)

### FILTER OUT ENTRIES IN DBSNP ---
df.target <- subset(df.target, !(df.target$key %in% common.key))

### REPLACE EXONIC FUNCTION ---
df.target$ExonicFunc.ensGene[which(df.target$ExonicFunc.ensGene == "frameshift deletion;frameshift deletion")] <- "frameshift deletion"
df.target$ExonicFunc.ensGene[which(df.target$ExonicFunc.ensGene == "frameshift insertion;frameshift insertion")] <- "frameshift insertion"
df.target$Func.ensGene[which(df.target$Func.ensGene == "exonic;exonic")] <- "exonic"

### SPLIT AF and AO ---
df.target$AF <- unlist(lapply(str_split(df.target$AF_AO, "_"), function(x) as.numeric(x[1])))
df.target$AO <- unlist(lapply(str_split(df.target$AF_AO, "_"), function(x) as.numeric(x[2])))

### FILTER OUT VARIANTS WITH > 80% ALLELE FREQ ---
# VARIANT WITH ALLELE FREQ 100% ARE USUALLY GERMLINE 
# UNLESS TUMOR PURITY IS 100%
df.target <- subset(df.target, df.target$AF < 80)

### FILTER BASED ON GENE IN TUMOR SAMPLES ---
genes.tumor <- unique(df.tumor$Gene.ensGene)
genes.target <- unique(df.target$Gene.ensGene)

genes.rescue1 <- intersect(genes.tumor, genes.target)

### GENES REPORTED TO BE MUTATED IN --------------------------
# MESO-PeM: Alakus et al. 2015, AACR PROJECT GENIE,
# MESO-PM: TCGA-MESO, Bueno et al. 2016,
# CANCER GENE CENSUS ---
# RESULTS: genes rescued = 28 (~32 ENSG.IDS)
# 	MPHOSPH9, ABCB11, MOCOS, APMAP, TSC22D1, NUP98, 
# 	ZAP70, APOA1, NCAPH, USP9X, KRT34, DNAJB1, 
# 	MATN3, RNF17, MYO1F, FMN2, ERG, TMPRSS3, NOS3, 
# 	FAT3, SAAL1, TLN2, CSPG4, GPR149, TTLL11, ARID2, 
#	SPTAN1, CUX1

genes.rescue2 <- c("ENSG00000276582","ENSG00000073734","ENSG00000101474","ENSG00000118137","ENSG00000189079",
					"ENSG00000173546","ENSG00000257923","ENSG00000132002","ENSG00000157554","ENSG00000165323",
					"ENSG00000282908","ENSG00000155816","ENSG00000174948","ENSG00000262045","ENSG00000131737",
					"ENSG00000132031","ENSG00000075643","ENSG00000051825","ENSG00000142347","ENSG00000121152",
					"ENSG00000164867","ENSG00000110713","ENSG00000132972","ENSG00000166788","ENSG00000197694",
					"ENSG00000171914","ENSG00000160183","ENSG00000102804","ENSG00000175764","ENSG00000124486",
					"ENSG00000115085")

genes.rescue <- unique(c(genes.rescue1, genes.rescue2))

### FILTER GENES BY GENES RESCUE ---
df.target <- subset(df.target, df.target$Gene.ensGene %in% genes.rescue)


### TRIM DATA ---
items <- c("Chr","Start","End","Ref","Alt",
			"Func.ensGene","Gene.ensGene","GeneDetail.ensGene",
			"ExonicFunc.ensGene","AAChange.ensGene",
			"esp6500si_all","cg69","X1000g2012apr_all",
			"snp138","snp138NonFlagged","cosmic68","SampleID")

dat <- subset(df.target, select=items)

### WRITE OUTPUT ---
file.output <- file.path(dir.data, "somatic_mutation_filtered_by_germline_annovar2maf.tsv")
write.table(dat, file.output, sep="\t", row.names=F, col.names=T, quote=F)


