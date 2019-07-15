### DEFINE LIBRARIES ----
library("stringr")

### SET PATHS ---------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/03_enrichment_analysis")
dir.data <- file.path(dir.wrk, "data")
dir.enrichment <- file.path(dir.wrk, "enrichment/cnv_meso_pem_pm")
dir.plot <- file.path(dir.wrk, "plot")
dir.output <- file.path(dir.wrk, "output")
dir.cnv <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/08_compare_cnv/output")

### DEFINE FILES ---
file.cnv <- file.path(dir.cnv, "cnv_meso_pem_pm_diffexpr_results_supplementarytable.tsv")
files.enrich <- list.files(dir.enrichment, pattern=".txt", full.names=TRUE)


### FUNCTION: ---
prepareData <- function(file.enrich, file.cnv){
	# LOAD CNV ---
	df <- read.delim(file.cnv, header=TRUE, stringsAsFactors=FALSE)

	# LOAD ENRICHMENT DATA ---
	dat.enrich <- read.delim(file.enrich, header=TRUE, stringsAsFactors=FALSE)

	# PREPARE DATA ---
	list.dat <- list()
	for(i in 1:nrow(dat.enrich)){
		list.dat[[i]] <- data.frame(Category=dat.enrich$Category[i], Gene=unlist(str_split(dat.enrich$overlap.genes[i], ",")))
	}

	### COMPILE DATA ---
	dat <- do.call(rbind.data.frame, list.dat)
	dat$Category <- as.character(dat$Category)
	dat$Gene <- as.character(dat$Gene)

	### ADD CHROM LOC. ---
	dat$ChromLoc <- ""
	for(i in 1:nrow(dat)){
		dat$ChromLoc[i] <- df$ChromosomeLocation[which(df$Gene == dat$Gene[i])]
	}

	return(dat)
}

### CALL FUNCTION : REACTOME---
file.enrich <- files.enrich[6]
dat <- prepareData(file.enrich, file.cnv)
file.output <- file.path(dir.output, "cnv_meso_pem_pm_diffexpr_reactome.tsv")
write.table(dat, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
