### Load Libraries -----
library("stringr")
library("ggplot2")
library("ggExtra")
library("gridExtra")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.analysis <- file.path(dir.wrk, "task/07_compare_mRNA_protein")
dir.output <- file.path(dir.analysis, "output")
dir.plot <- file.path(dir.analysis, "plot")
dir.script <- file.path(dir.analysis, "scripts")
dir.cna <- file.path(dir.wrk, "data/cnv/seq_call_refseq_genes_meso")
dir.mrna <- file.path(dir.wrk, "data/expression/analysis")
dir.prot <- file.path(dir.wrk, "data/proteome/processed_data")

### DEFINE FILES ---
file.des.mrna <- file.path(dir.wrk, "data/annotation/design_table_3p21genes.tsv")
file.des.prot <- file.path(dir.wrk, "data/proteome/processed_data/design_table_proteome.tsv")
file.cna <- file.path(dir.cna, "meso_cnv_seg_values_calls_parsed.tsv")
file.mrna <- file.path(dir.mrna, "meso_peritoneal_gene_expression_proteincoding.tsv.gz")
file.prot <- file.path(dir.prot, "meso_proteome_proteomediscover_all_log2_dqnorm.tsv.gz")

# ----------- CNV -----------
### LOAD CNV DATA ---
dat.cna <- read.delim(file.cna, header=T, stringsAsFactors=F, row.names=1)
colnames(dat.cna) <- str_replace_all(colnames(dat.cna), "[.]", "-")

### LOAD DESIGN TABLE: CNV ---
des.mrna <- read.delim(file.des.mrna, header=TRUE, stringsAsFactors=FALSE)
grp0 <- intersect(colnames(dat.cna), des.mrna$Sample.ID[which(des.mrna$Group.3p21 == 0)])
grp1 <- intersect(colnames(dat.cna), des.mrna$Sample.ID[which(des.mrna$Group.3p21 == 1)])

### ARRANGE SAMPLES: CNV ---
dat.cna <- subset(dat.cna, select=c(grp0, grp1))

# ----------- mRNA -----------
### LOAD DATA: mRNA ---
expr.mrna <- read.delim(file.mrna, header=TRUE, stringsAsFactors=FALSE, row.names=1)
colnames(expr.mrna) <- str_replace_all(colnames(expr.mrna), "[.]", "-")
expr.mrna[is.na(expr.mrna)] <- 0

### LOAD DESIGN TABLE: mRNA ---
des.mrna <- read.delim(file.des.mrna, header=TRUE, stringsAsFactors=FALSE)
grp0 <- intersect(colnames(expr.mrna), des.mrna$Sample.ID[which(des.mrna$Group.3p21 == 0)])
grp1 <- intersect(colnames(expr.mrna), des.mrna$Sample.ID[which(des.mrna$Group.3p21 == 1)])

### ARRANGE SAMPLES: mRNA ---
expr.mrna <- subset(expr.mrna, select=c(grp0, grp1))

# ----------- Protein -----------
### LOAD DATA: Protein ---
expr.prot <- read.delim(file.prot, header=TRUE, stringsAsFactors=FALSE, row.names=1)
colnames(expr.prot) <- str_replace_all(colnames(expr.prot), "[.]", "-")
expr.prot[is.na(expr.prot)] <- 0

### LOAD DESIGN TABLE: Protein ---
des.prot <- read.delim(file.des.prot, header=T, stringsAsFactors=F)
des.prot <- subset(des.prot, des.prot$SampleType == "Tumor")
grp0 <- des.prot$SequencingID[which(des.prot$Group.3p21 == 0)]
grp1 <- des.prot$SequencingID[which(des.prot$Group.3p21 == 1)]

### ARRANGE EXPRESSION DATA ---
expr.prot <- subset(expr.prot, select=c(grp0, grp1))
colnames(expr.prot) <- str_replace_all(colnames(expr.prot), "T", "")

### SAMPLEIDS ------
> colnames(dat.cna)
 [1] "MESO-11"  "MESO-15"  "MESO-18A" "MESO-18E" "MESO-01"  "MESO-08"
 [7] "MESO-04"  "MESO-19"  "MESO-12"  "MESO-13"  "MESO-03"  "MESO-06"
[13] "MESO-17"  "MESO-07"  "MESO-09"  "MESO-14"  "MESO-05"  "MESO-10"
[19] "MESO-02"
> colnames(expr.mrna)
 [1] "MESO-08"  "MESO-11"  "MESO-12"  "MESO-13"  "MESO-19"  "MESO-18A"
 [7] "MESO-18E" "MESO-02"  "MESO-09"  "MESO-10"  "MESO-14"  "MESO-17"
[13] "MESO-05"  "MESO-06"  "MESO-07"
> colnames(expr.prot)
 [1] "MESO-01"  "MESO-03"  "MESO-04"  "MESO-11"  "MESO-12"  "MESO-13"
 [7] "MESO-18A" "MESO-18E" "MESO-02"  "MESO-05"  "MESO-06"  "MESO-07"
[13] "MESO-09"  "MESO-10"  "MESO-14"  "MESO-17"

### GET COMMON SAMPLEIDS ---
ids <- Reduce(intersect, list(colnames(dat.cna),colnames(expr.mrna),colnames(expr.prot))) 
des <- subset(des.mrna, des.mrna$Sample.ID %in% ids)
des <- des[order(des$Group.3p21, decreasing=FALSE),]
ids <- des$Sample.ID

### SUBSET AND ARRANGE DATA ---
dat.cna <- subset(dat.cna, select=ids)
expr.mrna <- subset(expr.mrna, select=ids)
expr.prot <- subset(expr.prot, select=ids)


### FUNCTION: prepareData ---
prepareData <- function(dat, des, genes, type){
	# SUBSET AND RESHAPE DATA  ----
	dat <- subset(dat, rownames(dat) %in% genes)
	
	# ORDER GENES ---
	dat <- dat[match(genes, rownames(dat)),]	

	# MELT DATA ---
	df <- melt(t(as.matrix(dat)))
	colnames(df) <- c("SampleID","Gene","Value")	

	# ADD ATTRIBUTES ---
	df$SampleID <- as.character(df$SampleID)
	df$Gene <- as.character(df$Gene)
	df$Group <- ""
	df$Type <- type

	# GET GROUPS ---
	grp0 <- des$Sample.ID[which(des$Group.3p21 == 0)]
	grp1 <- des$Sample.ID[which(des$Group.3p21 == 1)]
	df$Group[which(df$SampleID %in% grp0)] <- "BAP1-INTACT"
	df$Group[which(df$SampleID %in% grp1)] <- "BAP1-DEL"

	return(df)
}

### FUNCTION: combineData ---
combineData <- function(dat1, dat2, ids, genes){
	df <- data.frame(SampleID=dat1$SampleID,
					Gene=dat1$Gene,
					Group=dat1$Group,
					Type=paste(dat1$Type, dat2$Type, sep=":"),
					Value1=dat1$Value,
					Value2=dat2$Value)

	df$SampleID <- factor(df$SampleID, levels=ids)
	df$Gene <- factor(df$Gene, levels=genes)
	df$Group <- factor(df$Group, levels=c("BAP1-INTACT","BAP1-DEL"))
	df$Type <- as.factor(df$Type)

	return(df)
}

### FUNCTION: prepareComboData ---
prepareComboData <- function(dat1, dat2, des, type1, type2){
	require("reshape2")

	# GET COMMON GENES ---
	genes <- intersect(rownames(dat1), rownames(dat2))
	ids <- colnames(dat1)

	# GET DATA ---
	df1 <- prepareData(dat=dat1, des=des, genes=genes, type=type1)
	df2 <- prepareData(dat=dat2, des=des, genes=genes, type=type2)

	# COMBINE DATA ----
	df <- combineData(dat1=df1, dat2=df2, ids, genes)

	return(df)
}


### COMBINE DATA: CNA-mRNA --- 
df.cna_mrna <- prepareComboData(dat1=dat.cna, dat2=expr.mrna, des=des, type1="CNA", type2="mRNA")
df.cna_prot <- prepareComboData(dat1=dat.cna, dat2=expr.prot, des=des, type1="CNA", type2="PROTEIN")


### FUNCTION: getCorrelation ---
getCorrelation <- function(dat){
	genes <- unique(dat$Gene)
	type <- unique(dat$Type)

	# PREPARE DATA  ---
	dat.cor <- data.frame(Gene=genes)
	dat.cor$Gene <- as.character(dat.cor$Gene)
	dat.cor$Type <- type
	dat.cor$R <- 0
	dat.cor$pvalue <- 0

	# GET CORRELATION ---
	for(i in 1:length(genes)){
		df <- subset(dat, dat$Gene == genes[i])

		x_na <- length(which(is.na(df$Value1)))
		y_na <- length(which(is.na(df$Value2)))

		x0 <- length(which(df$Value1 == 0))
		y0 <- length(which(df$Value2 == 0))

		if((x_na > 9) | (y_na > 9)) next
		if((x0 > 9) | (y0 > 9)) next

		correlation.test <- cor.test(x=df$Value1, y=df$Value2, method="pearson")
		dat.cor$R[i] <- as.numeric(correlation.test$estimate)
		pvalue <- as.numeric(correlation.test$p.value)
	
		if(pvalue == 0){
			dat.cor$pvalue[i] <- as.numeric(2.2e-16)
		} else{
			dat.cor$pvalue[i] <- pvalue
		}
	}

	return(dat.cor)
}

### CALL FUNCTION: getCorrelation -----
dfcor.cna_mrna <- getCorrelation(dat=df.cna_mrna)
dfcor.cna_prot <- getCorrelation(dat=df.cna_prot)

### WRITE OUTPUT ---
file.output1 <- file.path(dir.output, "meso_correlation_cna_mrna.tsv")
file.output2 <- file.path(dir.output, "meso_correlation_cna_protein.tsv")

write.table(dfcor.cna_mrna, file.output1, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(dfcor.cna_prot, file.output2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)



### RE-LOAD FILES ----
file.dat1 <- file.path(dir.output, "meso_correlation_cna_mrna.tsv")
file.dat2 <- file.path(dir.output, "meso_correlation_cna_protein.tsv")

dat1 <- read.delim(file.dat1, header=TRUE, stringsAsFactors=FALSE)
dat2 <- read.delim(file.dat2, header=TRUE, stringsAsFactors=FALSE)

### GET COMMON GENES ---
genes <- intersect(dat1$Gene, dat2$Gene)

### SUBSET BY GENES ---
dat1 <- subset(dat1, dat1$Gene %in% genes)
dat2 <- subset(dat2, dat2$Gene %in% genes)

### COMBINE DATA ---
df <- data.frame(Gene=dat1$Gene,
				Type1=dat1$Type,
				R1=dat1$R,
				pvalue1=dat1$pvalue,
				Type2=dat2$Type,
				R2=dat2$R,
				pvalue2=dat2$pvalue)

df$Gene <- as.character(df$Gene)
df$Type1 <- as.character(df$Type1)
df$Type2 <- as.character(df$Type2)

### REMOVE GENES WITH NO DATA ---
ind.del <- which((df$R1 == 0) | (df$R2 == 0))
df <- df[-ind.del,]

#> dim(df)
#[1] 7462    7

### COMPUTE ABSOLUTE DIFFERENCE BETWEEN R1-R2
df$Rdiff <- apply(df, 1, function(x) as.numeric(x[3]) - as.numeric(x[6]))
df$Rdiff.group <- ifelse(df$Rdiff >= 0.45, 1, 0)

#> table(df$Rdiff.group)
#
#   0    1
#5591 1871

### TAG GENES ---
tag.genes <- c("APC","CHEK1","EGFR","HDAC7","MAP3K4","MTAP","NF2","PBRM1","PIK3CA","RAD50","SETD2","SMARCC1")
df$TagGenes <- 0
df$TagGenes[which(df$Rdiff.group == 1)] <- 1
df$TagGenes[which(df$Gene %in% tag.genes)] <- 2

### ORDER DATA ---
df <- df[order(df$Rdiff.group, df$TagGenes),]

### ADD GENE LABEL ---
df$LabelGenes <- ""
df$LabelGenes[which(df$TagGenes == 2)] <- as.character(df$Gene[which(df$TagGenes == 2)])

### WRITE OUTPUT ---
file.output <- file.path(dir.output, "meso_correlation_cna-mrna_cna-protein.tsv")
write.table(df, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


### SUBSET PROTEIN ATTENUATED GENES ---
d <- subset(df, df$Rdiff.group == 1)
file.output <- file.path(dir.output, "meso_protein_attenuated_genes.tsv")
write.table(d$Gene, file.output, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


### FUNCTION: ----
get.plot <- function(df){
	require("ggplot2")
	require("ggExtra")
	cbPalette1 <- c("#98a3a5","#e85748")
	cbPalette2 <- c("#98a3a5","#e85748","#000000")

	# PLOT ---
	p <- ggplot(df, aes(x=R1, y=R2, label=LabelGenes)) + 
			#geom_point(shape = 21, fill="white", color="black", size=2, stroke=1, alpha=0.9) +
			geom_point(aes(fill=as.factor(Rdiff.group), color=as.factor(TagGenes)), shape = 21, stroke=0.5, size=1, alpha=0.4) +
			scale_fill_manual(values=cbPalette1) +
			scale_color_manual(values=cbPalette2) +
			#stat_smooth(method="lm", geom = "smooth", formula = y ~ x, position = "identity", fullrange = FALSE, se=TRUE) +
			geom_density2d(aes(group=as.factor(Rdiff.group)), stat = "density2d", lineend = "round", linejoin = "round", color="#FFFFFF", alpha=0.9, size=0.25) +
			coord_cartesian(xlim=c(-1, 1), ylim=c(-1, 1)) +
			geom_abline(intercept = 0, linetype = 2, color="#969696") +
			#geom_text(size=1, color="#000000", hjust=0, vjust=0) +
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
			ylab("Copy-number - Proteome") + 
			xlab("Copy-number - Transcriptome") + 
			ggtitle("")

	# add marginal density
	p <- ggExtra::ggMarginal(p, type="density", size = 4, aes(group=as.factor(Rdiff.group)), alpha=0.5) 

	return(p)	
}


### FOR MARGINAL DENSITY PLOT ----
get.densityplot <- function(df, type){
	require("ggplot2")
	cbPalette <- c("#98a3a5","#e85748")

	# SELECT VALUES ---
	colnames(df)[which(colnames(df) == type)] <- "R"

	g <- ggplot(df, aes(R, fill = as.factor(Rdiff.group))) + 
			geom_density(aes(color = as.factor(Rdiff.group)), alpha=0.4, size=0.5) +
			scale_fill_manual(values=cbPalette) +
			scale_color_manual(values=cbPalette) +
			coord_cartesian(xlim=c(-1, 1)) +
			theme(
				axis.text = element_text(size = 7, color="black"),
				axis.title = element_text(size = 7, color="black"),
				strip.text = element_text(size = 7, color="black"),
				strip.background = element_rect(fill = "white", colour = "white"),
				plot.title = element_text(size = 7, color="black", hjust=0),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.ticks = element_line(size=0.4, color="black"),	
				panel.background = element_rect(fill="white", color="white"),
				legend.position="none") +
			ylab("") +
			xlab("") + 
			ggtitle("")

	return(g)
}


### PLOT ---
file.plot <- file.path(dir.plot, "meso_correlation_cna-mrna_cna-protein.pdf")
pdf(file.plot, height=4, width=4)
	grid.arrange(get.plot(df), ncol=1, nrow=1)
	grid.arrange(get.densityplot(df, type="R1"), 
				get.densityplot(df, type="R2"),
				ncol=1, nrow=3)
dev.off()


