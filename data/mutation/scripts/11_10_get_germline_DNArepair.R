### DEFINE LIBRARIES ---
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/mutation")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.script <- file.path(dir.wrk, "scripts")
dir.maf <- file.path(dir.wrk, "maf_ver2")

### DEFINE FILES ---
file.dat <- file.path(dir.maf, "meso_pem_benign_tumor_cds.maf.gz")
file.genelist <- file.path(dir.maf, "genelist_dna_repair.txt")

### LOAD GENELIST ---
genelist <- read.delim(file.genelist, header=F, stringsAsFactors=F)$V1

### LOAD DATA ---
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)

### CHROM ---
dat$Chromosome <- str_replace_all(dat$Chromosome, "chr", "")


### GROUP SAMPLEIDS ---
ids <- unique(dat$Tumor_Sample_Barcode)
ids.normal <- ids[which(str_detect(ids, "N") == TRUE)]
ids.tumor <- ids[which(str_detect(ids, "T") == TRUE)]

### SEPARATE DATA BY SAMPLE TYPE ---
dat.normal <- subset(dat, dat$Tumor_Sample_Barcode %in% ids.normal)
dat.tumor <- subset(dat, dat$Tumor_Sample_Barcode %in% ids.tumor)

### SUBSET BY GENELIST ---
dat.germline <- subset(dat.normal, dat.normal$Hugo_Symbol %in% genelist)

### WRITE OUTPUT ---
file.output <- file.path(dir.maf, "meso_pem_germline_dnarepair.maf")
write.table(dat.germline, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
