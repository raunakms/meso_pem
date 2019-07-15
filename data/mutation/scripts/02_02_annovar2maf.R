### DEFINE LIBRARIES ---
library("maftools")

### DEFINE PATH ---
dir.wrk <- file.path("/collinsgroup/Raunak/projects/MESO_peritoneal/data/mutation")
dir.data <- file.path(dir.wrk, "main_calls")

### DEFINE FILES ---
file.dat <- file.path(dir.data, "somatic_mutation_non_silent_samplewise_input_for_annovar2maf.tsv")

### FUNCTION: get.annovar2maf ----
df <- annovarToMaf(annovar = file.dat , 
					Center = 'VPC-LAGA', 
					refBuild = 'hg19',
					tsbCol = 'SampleID', 
					table = 'ensGene')

df <- df[,c(1:26,29)]
df <- df[order(df$Hugo_Symbol, decreasing=FALSE),]


### WRITE OUTPUT ---
file.output <- file.path(dir.data, "VPCLAGA.MESO_PeM.nonsilent_mutation.maf")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)

# GENE --
y <- unique(df$Hugo_Symbol)
write.table(y, file.path(dir.data, "genelist_meso_pem_VPC.txt"), row.names=FALSE, col.names=FALSE, quote=F)

### CREATE MAF FILE ---
file.maf <- file.path(dir.data, "VPCLAGA.MESO_PeM.nonsilent_mutation.maf")
dat.maf <- read.delim(file.maf, header=TRUE, stringsAsFactors=FALSE)
dat.maf <- subset(dat.maf, dat.maf$AAChange != ".")

# WRITE OUTPUT
file.maf <- file.path(dir.data, "VPCLAGA.MESO_PeM.nonsilent_mutation_for_titv.maf")
write.table(dat.maf, file.maf, sep="\t", row.names=F, col.names=T, quote=F)



### FUNCTION ---
get.key <- function(dat){
	key <- apply(dat, 1, function(x) paste( as.character(x[5]), as.numeric(x[6]), as.numeric(x[7]), 
											as.character(x[12]), as.character(x[13]), sep=":"))	
	return(key)
}

### CALL FUNCTION ---
df$key <- get.key(df)

u.mut <- unique(df$key)
d <- df[,c(28,26)]
d <- d[!duplicated(d),]

#> length(which(d$cosmic68 == "."))
#[1] 313
