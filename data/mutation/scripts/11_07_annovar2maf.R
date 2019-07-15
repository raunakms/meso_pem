### DEFINE LIBRARIES ---
library("maftools")
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/mutation")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.script <- file.path(dir.wrk, "scripts")
dir.maf <- file.path(dir.wrk, "maf_ver2")
dir.annovar <- file.path( "/home/collins/IonProton/Mesothelioma/annovar")

### DEFINE PATH ---
files <- list.files(dir.annovar, pattern="input_for_annovar_to_maf_", full.names=FALSE)


### FUNCTION: get.annovar2maf ----
get.annovar2maf <- function(file.dat){
	# GET MAF ---
	df <- annovarToMaf(annovar = file.dat , 
						Center = 'LAGA', 
						refBuild = 'hg19',
						tsbCol = 'SampleID', 
						table = 'refGene')

	# REMOVE GENES WITH NA ----
	index <- which(is.na(df$Hugo_Symbol))

	if(length(index) != 0){
		df <- df[-index,]	
	}
	
	# SORT BY GENE SYMBOL ----
	df <- df[order(df$Hugo_Symbol, decreasing=FALSE),]

	return(df)
}


### GENERATE INPUT FILE FOR ANNOVAR2MAF ---
list.df <- list()
for(i in 1:length(files)){
	file.dat <- file.path(dir.annovar, files[i])

	# CALL FUNCTION ---
	list.df[[i]] <- get.annovar2maf(file.dat)

 	cat("PROCESSED:", i, "OF", length(files), "\n", sep=" ")
}


### AGGREGRATE DATA ---
df <- do.call(rbind.data.frame, list.df)

# WRITE OUTPUT ---
file.output <- file.path(dir.maf, "meo_pem_benign_tumor_all.maf")
write.table(df, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# COMPRESS ---
cmd <- paste("gzip", file.output, sep=" ")
system(cmd)



### VARIANT ANNOTATION
#> unique(df$Variant_Classification)
# [1] "Intron"            "Silent"            "Missense_Mutation"
# [4] "5'UTR"             "3'UTR"             "Frame_Shift_Ins"
# [7] "Nonsense_Mutation" "Frame_Shift_Del"   "5'Flank"
#[10] "In_Frame_Ins"      "RNA"               "Splice_Site"
#[13] "IGR"               "UNKNOWN"           "Nonstop_Mutation"
#[16] "3'Flank"           "In_Frame_Del"      NA


items.cds <- c("Missense_Mutation", "Nonsense_Mutation", 
				"Splice_Site","Nonstop_Mutation",
				"Frame_Shift_Ins","Frame_Shift_Del",
				"In_Frame_Ins","In_Frame_Del")

items.noncds <- c("5'UTR","3'UTR","5'Flank")

### SEPARATE CDS & NON-CDS ---
df.cds <- subset(df, df$Variant_Classification %in% items.cds)
df.noncds <- subset(df, df$Variant_Classification %in% items.noncds)

### WRITE OUTPUT ---
file.output <- file.path(dir.maf, "meo_pem_benign_tumor_cds.maf")
write.table(df.cds, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
cmd <- paste("gzip", file.output, sep=" ")
system(cmd)

### WRITE OUTPUT ---
file.output <- file.path(dir.maf, "meo_pem_benign_tumor_noncds.maf")
write.table(df.noncds, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
cmd <- paste("gzip", file.output, sep=" ")
system(cmd)
