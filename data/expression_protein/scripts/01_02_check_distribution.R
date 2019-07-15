### DEFINE LIBRARIES ----
library("stringr")

### DEFINE PATH ---
#dir.wrk <- file.path("//jbrcsrv009/CollinsGroup/Raunak/projects/MESO_peritoneal/data/proteome")
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/proteome")
dir.data <- file.path(dir.wrk, "data_proteome_discover")
dir.analysis <- file.path(dir.wrk, "processed_data")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")

### DEFINE FILES ---
file.expr <- file.path(dir.analysis, "meso_proteome_proteomediscover_all_human_log2_dqnorm_znorm.tsv")

### LOAD DATA ---
expr <- read.delim(file.expr, header=T, stringsAsFactors=F, row.names=1)
colnames(expr) <- str_replace_all(colnames(expr), "[.]", "-")

### COUNT NA ---
y.na <- rep(0, nrow(expr))
for(i in 1:nrow(expr)){
	y.na[i] <- length(which(is.na(expr[i,]) == TRUE))
}

genes0 <- rownames(expr)[which(y.na == 0)]
genes <- c(sample(genes0, 50))

### FUNCTION: GET DISTRIBUTION ---
get.distribution <- function(expr, genes, dir.plot){
	file.plot <- file.path(dir.plot, paste("protein_expr_distribution_meso.pdf", sep=""))
	pdf(file.plot, height=4, width=4)
	for(i in 1:length(genes)){
		gene <- genes[i]
		h <- hist(as.numeric(expr[gene,]), plot=FALSE) 
		plot(h, 
			main=gene, 
		 	xlab="protein expression", ylab="No. of Samples",
		 	cex.main=0.8, cex.lab=0.7, cex.axis=0.7, las=1, tck=-0.03)
	}
	dev.off()
}

### CALL FUNCTION: ---
get.distribution(expr, genes, dir.plot)


