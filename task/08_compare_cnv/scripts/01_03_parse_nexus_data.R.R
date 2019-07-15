### LOAD LIBRAIES ---
library("stringr")
library("VennDiagram")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/08_compare_cnv")
dir.data <- file.path(dir.wrk, "data")
dir.output <- file.path(dir.wrk, "output")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")
dir.ref <- file.path("/home/Collins/databases/Chr_gff")

### DEFINE FILES ---
file.nexus <- file.path(dir.data, "nexus_comparision_meso_ov.tsv")
file.features.wcox <- file.path(dir.output, "cnv_meso_ov_diff_summary_table_by_chromosome.tsv")
file.ref <- file.path(dir.ref, "Homo_sapiens.GRCh38.87_gene_proteincoding.bed.gz")

### LOAD REFERENCE GENOME ---
dat.ref <- read.delim(file.ref, header=F, stringsAsFactors=F)

### LOAD FEATURES ---
dat.features.wcox <- read.table(file.features.wcox, header=T, stringsAsFactors=F)
genes.wcox <- dat.features.wcox$Gene
genelist.background <- unique(dat.ref$V4)

### LOAD FILE ---
dat <- read.delim(file.nexus, header=T, stringsAsFactors=F)

### FORMAT DATA ---
dat$Region <- str_replace_all(dat$Region, ",", "")
dat$Chr <- unlist(lapply(str_split(dat$Region, ":"), function(x) x[1]))
dat$Start <- unlist(lapply(str_split(unlist(lapply(str_split(dat$Region, ":"), function(x) x[2])), "-"), function(x) as.numeric(x[1])))
dat$End <- unlist(lapply(str_split(unlist(lapply(str_split(dat$Region, ":"), function(x) x[2])), "-"), function(x) as.numeric(x[2])))
dat$Cytoband <- unlist(lapply(str_split(dat$Cytoband.Location, " - "), function(x) paste(x, collapse="-")))
dat$Gene.Symbols <- str_replace_all(dat$Gene.Symbols, ", ",",")

### FILTER DATA ---
dat <- dat[,c(15,16,17,2,13,3,4,5,6,7,9,10,11)]

### WRITE OUTPUT ---
file.output <- file.path(dir.data, "nexus_comparision_meso_ov_parsed.tsv")
write.table(dat, file.output, sep="\t", row.names=F, col.names=T, quote=F) 


### RELOAD DATA ---
file.dat <- file.path(dir.data, "nexus_comparision_meso_ov_parsed.tsv")
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)

dat <- subset(dat, dat$q.bound <= 1e-04)

### GET UNIQUE SET OF GENES ---
genes.nexus <- unique(unlist(str_split(dat$Gene.Symbols, ",")))
genes.nexus <- subset(genes.nexus, genes.nexus %in% genelist.background)

### INTERSECT FEATURES ----
genes.common.features <- intersect(genes.wcox, genes.nexus)
file.output <- file.path(dir.output, "cnv_meso_ov_diff_feature_genelist_wcox_nexus.txt")
write.table(genes.common.features, file.output, row.names=F, col.names=F, quote=F) 

### GET SUBSET DATA ----
df <- subset(dat.features.wcox, dat.features.wcox$Gene %in% genes.common.features)

#### WRITE OUTPUT -----------
file.output <- file.path(dir.output, "cnv_meso_ov_diff_summary_table_by_chromosome_wcox_nexus.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)

# Venn Diagram ---
list.genes <- list(WilcoxRank=genes.wcox, NexusCN=genes.nexus)
pdf(file.path(dir.plot, "venn_genelist_CNA_features.pdf"))
p1 <- venn.diagram(x=list.genes, scaled = TRUE,
		filename = NULL, lwd=2.5, cex=1.5, fontfamily="Helvetica",
		height=2, width=2,  units="inches",
		col="black", fill=c("firebrick1","deepskyblue"),
		main="", main.cex=1.5, main.fontfamily="Helvetica",
		cat.fontfamily="Helvetica", cat.cex=1, cat.default.pos="text",
		cat.pos=c(0, 0),
		#cat.dist=c(0.1,0.05), cat.just = list(c(0.5, 0.5), c(1,1))
		)
grid.draw(p1)
dev.off()
