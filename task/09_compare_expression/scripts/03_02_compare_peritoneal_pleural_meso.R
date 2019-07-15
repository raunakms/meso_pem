### LOAD LIBRARIES ----
library("stringr")

#### DEFINE PATH ----
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/09_compare_expression")
dir.data <- file.path(dir.wrk, "data")
dir.plot <- file.path(dir.wrk, "plot")

### DEFINE FILE ---
file.dat <- file.path(dir.data, "expr_combined_meso_pleural_peritoneal.tsv.gz")
file.des <- file.path(dir.data, "expr_combined_meso_pleural_peritoneal.txt")
file.diffexpr <- file.path("/Data/Raunak/softwares/bdvtools/array_process/wilcox_rank_test.R")

## Remove Genes with Expr = 0 in > 50% of the samples ------
remove.0.matrix <- function(dat, cutoff){
	y <- apply(dat, 1, function(x) length(which(x == 0)))
	del.index <- which(y > (ncol(dat) * cutoff))
	if(length(del.index) != 0){
		dat <- dat[-del.index,]
	} else{
		return(dat)
	}
	return(dat)
}

### LOAD DATA ---
dat <- read.delim(file.dat, header=T, stringsAsFactors=F, row.names=1)
colnames(dat) <- str_replace_all(colnames(dat), "[.]", "-")
dat <- remove.0.matrix(dat=dat, cutoff=0.25)

### SAMPLE CLASS ---
sample.class <- read.delim(file.des, header=F, stringsAsFactors=F)$V1
des <- data.frame(SampleID=colnames(dat), Class=sample.class)
des$SampleID <- as.character(des$SampleID)
des$Class <- as.character(des$Class)

sampleids.pleu <- des$SampleID[which(des$Class == "MESO-PLEURAL")]
sampleids.peri <- des$SampleID[which(des$Class == "MESO-PERITONEAL")]



#### DIFFERENTIAL ANALYSIS ----
source(file.diffexpr)
df <- get.wilcox.rank.test(dat, class1=sampleids.pleu, class2=sampleids.peri)
dfx <- subset(df, df$fdr <= 1e-04)
genes.dg <- df$Gene

