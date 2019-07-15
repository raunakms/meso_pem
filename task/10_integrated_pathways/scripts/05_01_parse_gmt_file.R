### Parse GMT Genesets ----------------
parseGMT <- function(file.geneset){
	# Load GeneSet Data
	no_col <- max(count.fields(file.geneset,  sep="\t"))
	dat.geneset <- read.delim(file.geneset, header=FALSE, fill=TRUE, col.names=1:no_col, as.is=paste("X",1:no_col, sep=""))		
	
	dat <- data.frame(Category=dat.geneset[,1])
	dat$Genesets <- ""
	#dat$Description <- ""
	for(i in 1:nrow(dat)){
		elements <- as.character(dat.geneset[i,3:ncol(dat.geneset)])
		elements <- paste(elements[which(elements != "")], collapse=":")
  
		dat$Genesets[i] <- elements
		#dat$Description[i] <- dat.geneset[i,2]
	} 

	dat$Category <- as.character(dat$Category)
	dat$Genesets <- as.character(dat$Genesets)
	
	return(dat)
}

### DEFINE LIBRARIES ---
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.task <- file.path(dir.wrk, "task/10_integrated_pathways")
dir.data <- file.path(dir.task, "data")
dir.analysis <- file.path(dir.task, "output")
dir.plot <- file.path(dir.task, "plots")

### DNA REPAIR ---
file.geneset <- file.path(dir.data, "genesets_DNA_Repair.gmt")
dat <- parseGMT(file.geneset)
write.table(dat, file.path(dir.data, "genesets_DNA_Repair.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

### CHROMATIN MODIFICATION ---
file.geneset <- file.path(dir.data, "genesets_Chromatin_Modification.gmt")
dat <- parseGMT(file.geneset)
write.table(dat, file.path(dir.data, "genesets_Chromatin_Modification.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
