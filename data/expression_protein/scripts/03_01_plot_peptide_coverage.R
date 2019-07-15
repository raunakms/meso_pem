### DEFINE LIBRARIES ----
library("stringr")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/proteome")
dir.data <- file.path(dir.wrk, "data_proteome_discover")
dir.analysis <- file.path(dir.wrk, "processed_data")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")
dir.task <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/07_compare_mRNA_protein/output")

### DEFINE FILES ---
file.dat <- file.path(dir.data, "MesoData_YY-DB-prot_only.converted.tsv")
file.des <- file.path(dir.data, "design_table_proteome.tsv")
file.attn <- file.path(dir.task, "meso_protein_attenuated_genes.tsv")

### LOAD DESIGN TABLE ---
des <- read.delim(file.des, header=T, stringsAsFactors=F)
ids.normal <- sort(des$SequencingID[which(des$SampleType == "Normal")], decreasing=FALSE)
ids.tumor <- sort(des$SequencingID[which(des$SampleType == "Tumor")], decreasing=FALSE)
ids.celline <- sort(des$SequencingID[which(des$SampleType == "Celline")], decreasing=FALSE)

### LOAD ATTENUATED PROTEIN LIST ---
proteins.attn <- read.table(file.attn, header=F, stringsAsFactors=F)$V1

### LOAD DATA ---
dat <- read.delim(file.dat, header=T, stringsAsFactors=F, skip=1)

#> str(dat)
#'data.frame':   8901 obs. of  133 variables:
# $ Checked               : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
# $ Master                : chr  "Master Protein" "Master Protein" "Master Protein" "Master Protein" ...
# $ Accession             : chr  "ENSP00000346879.5" "ENSP00000268057.4" "ENSP00000427271.1" "ENSP00000323036.8" ...
# $ Description           : chr  "NKX2-1" "BBS4" "SLC23A1" "ACPP" ...
# $ Coverage              : num  2 7.32 3.55 2.39 6.08 ...
# $ X..Peptides           : int  1 1 1 1 1 1 1 1 1 1 ...
# $ X..PSMs               : int  2 3 3 1 1 7 4 1 1 1 ...
# $ X..Unique.Peptides    : int  1 1 1 1 1 1 1 1 1 1 ...
# $ X..Protein.Groups     : int  1 1 1 1 1 1 1 1 1 1 ...
# $ X..AAs                : int  401 519 169 418 263 682 101 468 388 887 ...
# $ MW..kDa.              : num  41.7 58.2 18.9 48.3 29.6 74 11.3 53.9 42.8 99.2 ...
# $ calc..pI              : num  10.05 7.31 8.97 7.02 10.7 ...

### RETAIN CALLS WITH MASTER == "Master Protein"
dat <- subset(dat, dat$Master == "Master Protein")

### PREPARE DATA ---
df <- subset(dat, select=c("Accession","Description","Coverage"))

### FUNCTION: get.densityplot1 ----
get.densityplot1 <- function(df){
	require("ggplot2")

	g <- ggplot(df, aes(x=Coverage, fill="gray30")) + 
			geom_density(color=NA, fill="gray30", alpha=0.6, size=0.5) +
			geom_vline(xintercept=median(df$Coverage), colour="gray70", linetype=4, alpha=1) +
			theme(
				axis.text = element_text(size = 7, color="black"),
				axis.title = element_text(size = 7, color="black"),
				strip.text = element_text(size = 7, color="black"),
				strip.background = element_rect(fill = "white", colour = "white"),
				plot.title = element_text(size = 7, color="black", hjust=0),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.ticks = element_line(size=0.4, color="black"),	
				panel.background = element_rect(fill="white", color="black"),
				legend.position="none") +
			ylab("") +
			xlab("Protein Coverage") + 
			ggtitle("")

	return(g)
}

### PLOT ---
file.plot <- file.path(dir.plot, "meso_protein_coverage.pdf")
pdf(file.plot, height=2, width=4)
	grid.arrange(get.densityplot1(df), ncol=1, nrow=1)
dev.off()



### FUNCTION: get.densityplot1 ----
get.densityplot2 <- function(df, proteins.attn){
	require("ggplot2")
	cbPalette <- c("#98a3a5","#e85748")

	# LABEL ATTN DATA ---
	df$Group <- "NON-ATTN"
	df$Group[which(df$Description %in% proteins.attn)] <- "ATTN"
	df$Group <- factor(df$Group, levels=c("NON-ATTN","ATTN"))

	# PLOT 
	g <- ggplot(df, aes(x=Coverage, fill=Group)) + 
			geom_density(aes(fill=Group), color=NA, alpha=0.4, size=0.5) +
			scale_fill_manual(values=cbPalette) +
			#scale_color_manual(values=cbPalette) +
			#geom_vline(xintercept=median(df$Coverage), colour="gray70", linetype=4, alpha=1) +
			theme(
				axis.text = element_text(size = 7, color="black"),
				axis.title = element_text(size = 7, color="black"),
				strip.text = element_text(size = 7, color="black"),
				strip.background = element_rect(fill = "white", colour = "white"),
				plot.title = element_text(size = 7, color="black", hjust=0),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.ticks = element_line(size=0.4, color="black"),	
				panel.background = element_rect(fill="white", color="black"),
				legend.position="none") +
			ylab("") +
			xlab("Protein Coverage") + 
			ggtitle("")

	return(g)
}

### PLOT ---
file.plot <- file.path(dir.plot, "meso_protein_coverage_bygroup.pdf")
pdf(file.plot, height=2, width=4)
	grid.arrange(get.densityplot2(df, proteins.attn), ncol=1, nrow=1)
dev.off()

