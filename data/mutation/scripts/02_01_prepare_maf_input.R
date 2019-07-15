### DEFINE LIBRARIES --
library("stringr")

### DEFINe PATH ---
dir.wrk <- file.path("/collinsgroup/Raunak/projects/MESO_peritoneal/data/mutation")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")

### DEFINE FILES ---
file.dat <- file.path(dir.data, "somatic_mutation_hg19_non_silent_non_dbsnp_trimmed.tsv")

### LOAD DATA ---
dat <- read.delim(file.dat, header=TRUE, stringsAsFactors=FALSE)
y <- str_replace_all(colnames(dat)[28:43], "[.]", "-")
colnames(dat)[28:43] <- y

### SEPARATE DATA BY SAMPLES ---
list.df <- list()
ctr <- 1
for(i in 28:43){
	sampleid <- colnames(dat)[i]
	index <- which(dat[,i] != ".")
	df <- dat[index, 1:16]
	df$SampleID <- sampleid
	list.df[[ctr]] <- df
	ctr <- ctr + 1
}

df <- do.call(rbind.data.frame, list.df)

### REPLACE EXONIC FUNCTION ---
df$ExonicFunc.ensGene[which(df$ExonicFunc.ensGene == "frameshift insertion;frameshift insertion")] <- "frameshift insertion"
df$Func.ensGene[which(df$Func.ensGene == "exonic;exonic")] <- "exonic"

### FOR SAMPLES WITHOUT NORMALS ---
file.dat2 <- file.path(dir.data, "somatic_mutation_filtered_by_germline_annovar2maf.tsv")
dat2 <- read.delim(file.dat2, header=TRUE, stringsAsFactors=FALSE)

### COMBINE DATA ---
df.combined <- rbind(df, dat2)

### WRITE OUTPUT ---
file.output <- file.path(dir.data, "somatic_mutation_non_silent_samplewise_input_for_annovar2maf.tsv")
write.table(df.combined, file.output, sep="\t", row.names=F, col.names=T, quote=F)
