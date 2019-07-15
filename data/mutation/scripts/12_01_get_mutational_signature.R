### DEFINE LIBRARIES ---
library("deconstructSigs")
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/mutation")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.script <- file.path(dir.wrk, "scripts")
dir.maf <- file.path(dir.wrk, "maf_ver2")
dir.annovar <- file.path( "/home/collins/IonProton/Mesothelioma/annovar")

### DEFINE FILE ----
file.dat <- file.path(dir.maf, "meso_pem_benign_tumor_all.maf.gz")

### LOAD DATA ---
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
ids <- unique(dat$Tumor_Sample_Barcode)
ids.tumor <- ids[which(str_detect(ids, "T") == TRUE)]
dat <- subset(dat, dat$Tumor_Sample_Barcode %in% ids.tumor)
dat$Tumor_Sample_Barcode <- str_replace(dat$Tumor_Sample_Barcode, "T", "")
#dat <- subset(dat, dat$Variant_Type == "SNP")

### SUBSET DATA ---
dat.mut <- subset(dat, select=c("Tumor_Sample_Barcode","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2"))
#dat.mut$Chromosome <- paste("chr", dat.mut$Chromosome, sep="")
colnames(dat.mut) <- c("Sample","chr","pos","ref","alt")

#> unique(dat.mut$Sample)
#[1] "MESO-04"  "MESO-10"  "MESO-12"  "MESO-17"  "MESO-18A" "MESO-18E" "MESO-19"

# MESO-04  MESO-10  MESO-12  MESO-17 MESO-18A MESO-18E  MESO-19
#       1        1        1        1     8424     8374     8823

#dat.mut <- subset(dat.mut, dat.mut$Sample %in% c("MESO-18A","MESO-18E","MESO-19"))

###
sigs.input <- mut.to.sigs.input(mut.ref = dat.mut, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt")

### Determine the signatures contributing to the samples
sample_1 <- whichSignatures(tumor.ref = sigs.input, 
                           signatures.ref = signatures.cosmic, 
                           sample.id = "MESO-19", 
                           contexts.needed = TRUE,
                           tri.counts.method = "default")

