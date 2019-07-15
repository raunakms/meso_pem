### LOAD LIBRAIES ---
library("stringr")
library("bigmemory")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/08_compare_cnv")
dir.data <- file.path(dir.wrk, "data")
dir.output <- file.path(dir.wrk, "output")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")
dir.ov <- file.path("/Data/Raunak/HITnDRIVE/datasets/TCGA_OV/data/cnv/seq_call_refseq_genes_ov")
dir.ref <- file.path("/home/Collins/databases/Chr_gff")

### DEFINE FILES ---
file.meso <- file.path(dir.data, "meso_cnv_seg_values_calls_parsed.tsv")
file.ov <- file.path(dir.ov, "tcga_ov_seg_values_calls_with_ids_unique.tsv")
file.sampleid <- file.path(dir.ov, "sample_ids.txt")
file.script <- file.path("/Data/Raunak/softwares/bdvtools/array_process/wilcox_rank_test.R")
file.ref <- file.path(dir.ref, "Homo_sapiens.GRCh38.87_gene_proteincoding.bed.gz")


## Remove Genes with Expr = NA in > 25% of the samples ------
remove.na.matrix <- function(dat, cutoff){
	y <- apply(dat, 1, function(x) length(which(is.na(x))))
	del.index <- which(y > (ncol(dat) * cutoff))
	if(length(del.index) != 0){
		dat <- dat[-del.index,]
	} else{
		return(dat)
	}
	return(dat)
}

### LOAD REFERENCE GENOME ---
dat.ref <- read.delim(file.ref, header=F, stringsAsFactors=F)

### LOAD MESO CNV ----
dat.meso <-  read.delim(file.meso, header=T, stringsAsFactors=F, row.names=1)
colnames(dat.meso) <- str_replace_all(colnames(dat.meso), "[.]", "-")
dat.meso <- remove.na.matrix(dat=dat.meso, cutoff=0.25)

### LOAD OV CNV ----
sampleid <-  read.delim(file.sampleid, header=F, stringsAsFactors=F)$V1
dat.ov <- read.big.matrix(filename=file.ov, sep="\t", header=TRUE, 
				col.names=sampleid, has.row.names=TRUE, ignore.row.names=FALSE,
				type="double", skip=0)

dat.ov <- as.matrix(dat.ov)
dat.ov <- remove.na.matrix(dat=dat.ov, cutoff=0.25)

#### Subset data by common genes ---
genes <- intersect(rownames(dat.meso), rownames(dat.ov))
dat.meso <- subset(dat.meso, rownames(dat.meso) %in% genes)
dat.ov <- subset(dat.ov, rownames(dat.ov) %in% genes)

### Combine Data ---
sampleids.meso <- colnames(dat.meso)
sampleids.ov <- colnames(dat.ov)

sample.class <- c(rep("MESO", length(sampleids.meso)), rep("OV", length(sampleids.ov)))

dat <- cbind(dat.meso, dat.ov)

#### WRITE OUTPUT -----------
file.output <- file.path(dir.data, "cnv_meso_ov_seg_mean_combined.tsv")
write.table(dat, file.output, sep="\t", row.names=T, col.names=NA, quote=F)

file.output <- file.path(dir.data, "cnv_meso_ov_combined_class.txt")
write.table(sample.class, file.output, sep="\t", row.names=F, col.names=F, quote=F)


#### DIFFERENTIAL CNA ANALYSIS ----
source(file.script)
df <- get.wilcox.rank.test(dat, class1=sampleids.meso, class2=sampleids.ov)

#### SUBSET ONLY PROTEINCODING GENES ----
df.proteincoding <- subset(df, df$Gene %in% unique(dat.ref$V4))
df.proteincoding <- subset(df.proteincoding, df.proteincoding$fdr <= 1e-04)

genes.dg <- df.proteincoding$Gene

#### WRITE OUTPUT -----------
file.output <- file.path(dir.output, "cnv_meso_ov_diff_summary_table.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)

file.deg <- file.path(dir.output, "cnv_meso_ov_diff_genelist.txt")
write.table(genes.dg , file.deg, sep="\t", row.names=F, col.names=F, quote=F)


####
file.annot <- file.path(dir.output, "biomart_gene_annotation.tsv")
file.dat <- file.path(dir.output, "cnv_meso_ov_diff_summary_table.tsv")

### LOAD BIOMART ANNOTATION FILE ----
annot <- read.delim(file.annot, header=T, stringsAsFactors=F)
annot <- annot[-which(str_detect(annot$Chromosome, "CHR_") == TRUE),]

### LOAD FILES ---
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)


### ADD ANNOTATION ---
dat$Chromosome <- ""
dat$Start <- ""
dat$End <- ""
dat$Band <- ""
for(i in 1:nrow(dat)){
	index <- which(annot$Gene == dat$Gene[i])
	if(length(index) == 0) next
	index <- index[1]
	dat$Chromosome[i] <- annot$Chromosome[index]
	dat$Start[i] <- annot$Start[index]
	dat$End[i] <- annot$End[index]
	dat$Band[i] <- annot$Band[index]

	#cat("i= ", i, "index= ", index, "\n", sep="\t")
}

### TRIM DATA ---
dat <- subset(dat, dat$Chromosome != "")
dat <- dat[with(dat, order(Chromosome, Band, Start)),]

items <- c("Gene","Chromosome","Start","End","Band",       
			"MedianA","StdevA","MedianB","StdevB",
			"pvalue","fdr")
dat <- subset(dat, select=items)

#### WRITE OUTPUT -----------
file.output <- file.path(dir.output, "cnv_meso_ov_diff_summary_table_by_chromosome.tsv")
write.table(dat, file.output, sep="\t", row.names=F, col.names=T, quote=F)
