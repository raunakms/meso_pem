### DEFINE LIBRARIES --
library("stringr")
#library("ggplot2")
#library("gridExtra")
library("maftools")

### DEFINe PATH ---
dir.fa <- file.path("//jbrcsrv009/CollinsGroup/Raunak/data_ref/hg_chromosome/GRCh37")
dir.wrk <- file.path("//jbrcsrv009/CollinsGroup/Raunak/projects/MESO_peritoneal/data/mutation")
#dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/mutation")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")


### DEFINE FILES ---
file.maf <- file.path(dir.data, "VPCLAGA.MESO_PeM.nonsilent_mutation.maf")
file.ref <- file.path(dir.fa, "Homo_sapiens.GRCh37.dna.alt.fa")

# LOAD DATA ---
#dat.maf <- read.delim(file.maf, header=TRUE, stringsAsFactors=FALSE)
dat.maf <- read.maf(file.maf, removeSilent = TRUE, useAll = TRUE,
						gisticAllLesionsFile = NULL, gisticAmpGenesFile = NULL,
						gisticDelGenesFile = NULL, cnTable = NULL,
						removeDuplicatedVariants = TRUE, isTCGA = FALSE, verbose = TRUE)

### ADD TRINUCLEOTIDE INFO ----
meso.tnm <- trinucleotideMatrix(maf = dat.maf, ref_genome = file.ref, 
                               prefix = 'chr', add = TRUE, ignoreChr = 'chr23', useSyn = TRUE)
