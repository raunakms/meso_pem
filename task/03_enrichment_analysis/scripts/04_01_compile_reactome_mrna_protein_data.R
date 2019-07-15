### LOAD LIBRARIES -----
library("stringr")

### SET PATHS --------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/03_enrichment_analysis")
dir.data <- file.path(dir.wrk, "data")
dir.plot <- file.path(dir.wrk, "plot")
dir.enrichment <- file.path(dir.wrk, "enrichment")

### DEFINE FILE ----------
file.mrna <- file.path(dir.data, "MESO_OUTLIERS_mRNA_p01_pathway_enrichment_reactome.tsv")
file.prot <- file.path(dir.data, "MESO_OUTLIERS_protein_p03_pathway_enrichment_reactome.tsv")

### LOAD FILES ---
dat.mrna <- read.delim(file.mrna, header=TRUE, stringsAsFactors=FALSE, row.names=1)
colnames(dat.mrna)  <- str_replace_all(colnames(dat.mrna), "[.]", "-")

dat.prot <- read.delim(file.prot, header=TRUE, stringsAsFactors=FALSE, row.names=1)
colnames(dat.prot) <- str_replace_all(colnames(dat.prot), "[.]", "-")

### GET UNIQUE PATHWAYS ---
p <- unique(c(rownames(dat.mrna), rownames(dat.prot)))
write.table(p, file.path(dir.data, "MESO_OUTLIERS_reactome_pathwaylist.txt"), sep="\t", row.names=F, col.names=F, quote=F)

### LOAD PATHWAY DESIGN TABLE ---
file.des <- file.path(dir.data, "design_table_pathway.tsv")
des <- read.delim(file.des, header=TRUE, stringsAsFactors=FALSE)

### SUBSET PATHWAYS ---
dat.mrna <- subset(dat.mrna, rownames(dat.mrna) %in% des$Pathway)
dat.prot <- subset(dat.prot, rownames(dat.prot) %in% des$Pathway)

### CREATE EMPTY MATRIX ----
mat.mrna <- matrix(0, nrow=nrow(des), ncol=ncol(dat.mrna), dimnames=list(des$Pathway, colnames(dat.mrna)))
mat.prot <- matrix(0, nrow=nrow(des), ncol=ncol(dat.prot), dimnames=list(des$Pathway, colnames(dat.prot)))

### FUNCTION: FILL MATRIX ---
fillmatrix <- function(des, mat, dat){
	for(i in 1:nrow(dat)){
		index <- which(rownames(mat) == rownames(dat)[i])
		mat[index,] <- as.numeric(dat[i,])
	}
	return(mat)
}

### CALL FUNCTION ---
mat.mrna <- fillmatrix(des=des, mat=mat.mrna, dat=dat.mrna)
mat.prot <- fillmatrix(des=des, mat=mat.prot, dat=dat.prot)

### WRITE OTUPTU ---
file.output <- file.path(dir.data, "MESO_OUTLIERS_reactome_pathway_filtered_mrna.tsv")
write.table(mat.mrna, file.output, sep="\t", row.names=T, col.names=NA, quote=F)

file.output <- file.path(dir.data, "MESO_OUTLIERS_reactome_pathway_filtered_protein.tsv")
write.table(mat.prot, file.output, sep="\t", row.names=T, col.names=NA, quote=F)

