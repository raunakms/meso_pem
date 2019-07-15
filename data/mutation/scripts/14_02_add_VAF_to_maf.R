### DEFINE LIBRARIES ---
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/collinsgroup/Raunak/projects/MESO_peritoneal/data/mutation")
dir.data <- file.path(dir.wrk, "main_calls")

### DEFINE PATH ---
file.maf <- file.path(dir.data, "VPCLAGA.MESO_PeM.nonsilent_mutation.maf")
file.vaf <- file.path(dir.data, "annovar_all_vaf.tsv")

### LOAD DATA ---
dat.maf <- read.delim(file.maf, header=TRUE, stringsAsFactors=FALSE)
dat.vaf <- read.delim(file.vaf, header=TRUE, stringsAsFactors=FALSE)


### ADD KEY ---
dat.maf$key <- apply(dat.maf, 1, function(x) paste(as.character(x[5]), as.numeric(x[6]), as.numeric(x[7]), as.character(x[15]), sep=":") )
dat.vaf$key <- apply(dat.vaf, 1, function(x) paste(as.character(x[1]), as.numeric(x[2]), as.numeric(x[3]), as.character(x[8]), sep=":") )

### ADD VAF ---
dat.maf$VAF <- NA
for(i in 1:nrow(dat.maf)){
	index <- which(dat.vaf$key == dat.maf$key[i])
	dat.maf$VAF[i] <- dat.vaf$FreqT[index[1]]
}

dat.maf <- dat.maf[,-which(colnames(dat.maf) == "key")]

### WRITE OUTPUT ---
file.output <- file.path(dir.data, "VPCLAGA.MESO_PeM.nonsilent_mutation_withVAF.maf")
write.table(dat.maf, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
