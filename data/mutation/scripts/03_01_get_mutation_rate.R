### DEFINE LIBRARIES --
library("stringr")

### DEFINe PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/mutation")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")
dir.hg <- file.path("/Data/Raunak/data_ref/hg_chromosome/hg38")
dir.script <- file.path("/Data/Raunak/softwares/bdvtools/NGS")

### DEFINE FILES ---
file.maf <- file.path(dir.data, "VPCLAGA.MESO_PeM.nonsilent_mutation.maf")
file.hg <- file.path(dir.hg, "Homo_sapiens_GRCh38_gene_chromInfo.tsv")
file.ngs <- file.path(dir.script, "get_mutation_Rate.R")

### FUNCTION: GET MUTATION RATE ---
get.mutation.rate <- function(file.maf, file.hg){
	# DEFINE LIBRARIES ---
	require("stringr")

	# LOAD CHROMOSOME DATA ---
	dat.hg <- read.delim(file.hg, header=TRUE, stringsAsFactors=FALSE)

	# LOAD DATA ---
	dat.maf <- read.delim(file.maf, header=TRUE, stringsAsFactors=FALSE)

	# GET SAMPLEID ---
	#dat.maf$SampleID <- unlist(lapply(lapply(str_split(dat.maf$Tumor_Sample_Barcode, "-"), function(x){paste(x[1],x[2],x[3],x[4], sep="-")}), function(x){str_sub(x, 1, 15)}))

	# SUBSET ONLY PROTEIN-CODING MUTATIONS ---
	#dat.maf <- subset(dat.maf, dat.maf$BIOTYPE == "protein_coding")
	#dat.maf <- subset(dat.maf, dat.maf$Variant_Classification != "Silent")

	# SUBSET FIELDS ---
	items <- c("Tumor_Sample_Barcode","Hugo_Symbol","Chromosome","Start_Position","End_Position","Strand","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2")
	dat.maf <- subset(dat.maf, select=items)
	dat.maf$Chromosome <- paste("chr", dat.maf$Chromosome, sep="")

	# GET SAMPLEIDS ---
	ids <- sort(unique(dat.maf$Tumor_Sample_Barcode), decreasing=FALSE)
	df <- data.frame(SampleID=ids)
	df$SampleID <- as.character(df$SampleID)
	df$Rate.SNV <- 0

	# GET MUTATION RATE---
	for(ctr in 1:length(ids)){
		df.hg <- dat.hg
		dat.temp <- subset(dat.maf, dat.maf$Tumor_Sample_Barcode %in% ids[ctr])

		# GET MUTATION BY CHR ---
		df.freq <- data.frame(chrom=names(table(dat.temp$Chromosome)), Freq=as.numeric(table(dat.temp$Chromosome)))
		df.freq$chrom <- as.character(df.freq$chrom)

		# ADD Freq TO CHR DATA ---
		df.hg$Freq <- 0
		for(j in 1:nrow(df.freq)){
			df.hg$Freq[which(df.hg$Chrom == df.freq$chrom[j])] <- df.freq$Freq[j]
		}

		# MUTATION RATE ---
		df$Rate.SNV[ctr] <- sum(df.hg$Freq)/sum(df.hg$Size.mb)
	}

	return(df)
}


### CALL FUNCTION ---
df <- get.mutation.rate(file.maf, file.hg)
df <- cbind(CancerType="MESO-PeM", df)

# WRITE OUTPUT ---
file.output <- file.path(dir.analysis, paste("snv_rate_", "MESO-PeM", ".tsv", sep=""))
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)
