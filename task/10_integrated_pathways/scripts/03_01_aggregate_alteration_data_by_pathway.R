### Load LIBRARIES -------
library("stringr")
library("reshape2")

### SetDirectories -----
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.task <- file.path(dir.wrk, "task/10_integrated_pathways")
dir.data <- file.path(dir.task, "data")
dir.analysis <- file.path(dir.task, "output")
dir.plot <- file.path(dir.task, "plots")

### DEFINE FILES ---
file.snv <- file.path(dir.wrk, "data/mutation/analysis/snv_alterations_status.tsv")
file.cnv <- file.path(dir.wrk, "data/cnv/analysis/alterations_samplewise_cna.tsv")
file.fus <- file.path(dir.wrk, "data/genefusion/analysis/alterations_samplewise_fusions.tsv")
file.pathway <- file.path(dir.data, "pathway_gene_relationship.tsv")

### LOAD DATA ---
dat.snv <- read.delim(file.snv, header=T, stringsAsFactors=F)[,1:3]
dat.cnv <- read.delim(file.cnv, header=T, stringsAsFactors=F)[,1:3]
dat.fus <- read.delim(file.fus, header=T, stringsAsFactors=F)[,1:3]

### FILTER ---
dat.snv <- subset(dat.snv, dat.snv$Gene != "")
dat.cnv <- subset(dat.cnv, dat.cnv$Gene != "")
dat.fus <- subset(dat.fus, dat.fus$Gene != "")

### LOAD PATHWAY DATA ---
dat.pathway <- read.delim(file.pathway, header=T, stringsAsFactors=F)

### GET UNIQUE PATHWAY ---
pathways <- unique(dat.pathway$Pathway)

### GET UNIQUE LIST OF SAMPLEIDS ---
ids <- unique(c(unique(dat.snv$SampleID), unique(dat.cnv$SampleID), unique(dat.fus$SampleID)))


### FUNCTION: getAlteration ---
getAlteration <- function(ids, dat.snv, dat.cnv, dat.fus, dat.pathway, p){
	# SUBSET PATHWAY ---
	d <- subset(dat.pathway, dat.pathway == p)

	# SUBSET ALTERATION ---
	alt1 <- subset(dat.snv, dat.snv$Gene %in% d$Gene)
	alt2 <- subset(dat.cnv, dat.cnv$Gene %in% d$Gene)
	alt3 <- subset(dat.fus, dat.fus$Gene %in% d$Gene)

	# MERGE DATA ---
	dat.alt <- rbind(alt1, alt2, alt3)

	# GET NULL IDS ---
	ids.null <- setdiff(ids, unique(dat.alt$SampleID))

	# CREATE NULL DATAFRAME ---
	if((length(ids.null) != 0) & (nrow(dat.alt) != 0)){
		dat.null <- data.frame(SampleID=ids.null, Gene=dat.alt$Gene[1], Status="NULL;")
		df.alt <- rbind(dat.alt, dat.null)
		df.alt <- cbind(Pathway=p, df.alt)
	} else if((length(ids.null) == 0) & (nrow(dat.alt) != 0)){
		df.alt <- cbind(Pathway=p, dat.alt)
	} else if((length(ids.null) == 0) & (nrow(dat.alt) == 0)){
		df.alt <- data.frame()
	}

	return(df.alt)
}

### CALL FUNCTION ---
list.pathway <- list()
for(i in 1:length(pathways)){
	p <- pathways[i]
	list.pathway[[i]] <- getAlteration(ids, dat.snv, dat.cnv, dat.fus, dat.pathway, p)
}

### COMBINE DATA ---
df <- do.call(rbind.data.frame, list.pathway)

### REPLACE STATUS VALUES ---
df$Status[which(df$Status == "INDEL;")] <- "IND;"
df$Status[which(df$Status == "SNP;")] <- "MUT;"
df$Status[which(df$Status == "SNP;INDEL;")] <- "MUT;IND;"
df$Status[which(df$Status == "NULL;")] <- "EMPTY;"

### WRITE OUTPUT ---
file.output <- file.path(dir.data, "aggregate_alteration_oncoprint_by_pathway.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)

### ONCOPRINT INPUT ---
gnames <- apply(df, 1, function(x) paste(as.character(x[1]), as.character(x[3]), sep=":"))
df$Gene <- gnames
df <- df[,-1]

### WRITE OUTPUT ---
file.output <- file.path(dir.data, "aggregate_alteration_oncoprint_by_pathway_input.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)

