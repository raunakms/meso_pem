### Load Libraries -----
library("stringr")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ----
dir.wrk <- file.path("/collinsgroup/Raunak/projects/MESO_peritoneal/data")
dir.cna <- file.path(dir.wrk, "cnv/seq_call_refseq_genes_meso")
dir.exp <- file.path(dir.wrk, "expression/analysis")
dir.plot <- file.path(dir.wrk, "cnv/plot")
dir.analysis <- file.path(dir.wrk, "cnv/analysis")

### DEFINE FILES ---
file.cna <- file.path(dir.cna, "meso_cnv_seg_values_calls_parsed.tsv")
file.exp <- file.path(dir.exp, "meso_peritoneal_gene_expression_proteincoding.tsv.gz")
file.des <- file.path("/collinsgroup/Raunak/projects/MESO_peritoneal/data/annotation/design_table_3p21genes.tsv")

### LOAD DESIGN TABLE ---
des <- read.delim(file.des, header=TRUE, stringsAsFactors=FALSE)
grp0 <- des$Sample.ID[which(des$Group.3p21 == 0)]
grp1 <- des$Sample.ID[which(des$Group.3p21 == 1)]

### LOAD CNV DATA ---
dat.cna <- read.delim(file.cna, header=T, stringsAsFactors=F, row.names=1)
colnames(dat.cna) <- str_replace_all(colnames(dat.cna), "[.]", "-")

### LOAD EXPRESSION  DATA ----
dat.exp <- read.delim(file.exp, header=T, stringsAsFactors=F, row.names=1)
colnames(dat.exp) <- str_replace_all(colnames(dat.exp), "[.]", "-")


### FUNCTION: GET CORRELATION --
get.corr <- function(des, dat.cna, dat.exp, genes){
	# GET SAMPLE GROUPS ---
	grp0 <- des$Sample.ID[which(des$Group.3p21 == 0)]
	grp1 <- des$Sample.ID[which(des$Group.3p21 == 1)]

	### MATCH SAMPLES IN CNA & EXP ----
	sampleids <- intersect(colnames(dat.cna), colnames(dat.exp))
	dat.cna <- subset(dat.cna, select=sampleids)
	dat.exp <- subset(dat.exp, select=sampleids)

	### SUBSET BY ASSIGNED GENES ---
	dat.cna <- subset(dat.cna, rownames(dat.cna) %in% genes)
	dat.exp <- subset(dat.exp, rownames(dat.exp) %in% genes)

	### CORRELATION OF TWO MATRICES ----
	dat.cor <- data.frame(Gene=genes)
	dat.cor$Gene <- as.character(dat.cor$Gene)
	dat.cor$R <- 0

	list.plot <- list()
	for(i in 1:nrow(dat.cna)){
		gene <- genes[i]
		x <- as.numeric(dat.cna[gene,])
		y <- as.numeric(dat.exp[gene,])

		r <- cor(x, y, method="pearson", use="pairwise.complete.obs")
		dat.cor$R[i] <- r

		# PREPARE DATA ----
		df <- data.frame(Gene=gene, 
					SampleID=unlist(lapply(str_split(colnames(dat.cna), "-"), function(x) x[2])), 
					SampleIDfull=colnames(dat.cna),
					CNA.value=x,
					Exp.value=y)
		df$Group <- ""
		df$Group[which(df$SampleIDfull %in% grp0)] <- "BAP1-intact"
		df$Group[which(df$SampleIDfull %in% grp1)] <- "BAP1-del"
		df$Group <- factor(df$Group, levels=c("BAP1-del","BAP1-intact"))

		list.plot[[i]] <- ggplot(df, aes(x=CNA.value, y=Exp.value)) + 
							stat_smooth(method="lm", geom="smooth", formula=y~x, color="red", position="identity", fullrange=FALSE, se=TRUE) +
							scale_fill_manual(values=c("blue","#fffc00")) +
							geom_point(aes(fill=Group), color="black", shape = 21, size=3, alpha=0.9) +
							#coord_cartesian(xlim=c(-5, 2), ylim=c(-5, 20)) +
							geom_text(aes(label=SampleID,  hjust=1, vjust=-1, lineheight=1.5), color="black", size=3) +
							#annotate("text", label = paste("R = ", round(r, 2), sep=""), x = 0, y = 13, size = 3, colour = "red") +
							facet_wrap(~ Gene, ncol = 1) +
							theme(
								axis.text = element_text(size = 5, color="black"),
								axis.title = element_text(size = 5, color="black"),
								plot.title = element_text(size = 5, color="black", hjust=0),
								panel.grid.major = element_blank(),
								panel.grid.minor = element_blank(),
								axis.ticks = element_line(size=0.4, color="black"),
								strip.text = element_text(size=5, color="black"),
								strip.background = element_rect(fill="white", color="white"),
								panel.background = element_rect(fill="white", color="black"),
								legend.position="none") +
							xlab("Copy Number Segment Mean") + 
							ylab("Gene Expression") + 
							ggtitle("") 

		cat("PLOT GENERATED", i, "OF", nrow(dat.cna), "\n", sep="\t")							
	}

	#return(dat.cor)
	return(list.plot)
}

#### CALL FUNCTION ---
genes <- c("BAP1","PBRM1","SETD2","SMARCC1")
p <- get.corr(des, dat.cna, dat.exp, genes)


### PLOT ---
file.plot <- file.path(dir.plot, "cnv_expr_corr_3p21genes.pdf")
pdf(file.plot, height=6, width=6)
	grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], ncol=2, nrow=2)
dev.off()

file.plot <- file.path(dir.plot, "cnv_expr_corr_3p21genes_refined.pdf")
pdf(file.plot, height=6, width=6)
	grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], ncol=2, nrow=2)
dev.off()


#> dat.cor
#     Gene         R
#1   SETD2 0.7161982
#2 SMARCC1 0.9011789
#3    BAP1 0.8967636
#4   PBRM1 0.7948309

### FOR CN GAIN GENES ---
#genes <- c("BRD9","SETDB1","MDM4","RICTOR","ETF1","FOXL1","EGFR","HRAS","CSPG4","STK11","NF2")
genes <- c("BRD9","EGFR","ETF1","NF2")
#genes <- sort(genes, decreasing=FALSE)
p <- get.corr(des, dat.cna, dat.exp, genes)


### PLOT ---
file.plot <- file.path(dir.plot, "cnv_expr_corr_cngain_refined.pdf")
pdf(file.plot, height=2, width=8)
	grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], ncol=4, nrow=1)
dev.off()



genes <- c("FOXL1","PDGFA")
p <- get.corr(des, dat.cna, dat.exp, genes)

### PLOT ---
file.plot <- file.path(dir.plot, "cnv_expr_corr_cngain_refined2.pdf")
pdf(file.plot, height=2, width=4)
	grid.arrange(p[[1]], p[[2]], ncol=2, nrow=1)
dev.off()


#     Gene            R
#1    BRD9  0.906220185
#2   CSPG4 -0.140371307
#3    EGFR  0.556834534
#4    ETF1  0.587899108
#5   FOXL1  0.620301130
#6    HRAS  0.134207439
#7    MDM4 -0.009422827
#8     NF2  0.792163546
#9  RICTOR  0.135549693
#10 SETDB1 -0.581994306
#11  STK11  0.316987285


### FOR CN GAIN GENES ---
genes <- c("BRD9","ETF1","FOXL1","EGFR","NF2")
genes <- sort(genes, decreasing=FALSE)
p <- get.corr(des, dat.cna, dat.exp, genes)

### PLOT ---
file.plot <- file.path(dir.plot, "cnv_expr_corr_cngain.pdf")
pdf(file.plot, height=6, width=6.8)
	grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], ncol=2, nrow=3)
dev.off()


### FOR CNV LOSS GENES ----
genes <- c("LATS1","CDKN2A","ERCC2","MAP2K2")
genes <- sort(genes, decreasing=FALSE)
p <- get.corrdes, (dat.cna, dat.exp, genes)


### PLOT ---
file.plot <- file.path(dir.plot, "cnv_expr_corr_cnloss.pdf")
pdf(file.plot, height=6, width=6)
	grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], 
				ncol=2, nrow=2)
dev.off()

#    Gene         R
#1 CDKN2A 0.6997592
#2  ERCC2 0.1316195
#3  LATS1 0.8692454
#4 MAP2K2 0.3418527
#5   TP53 0.5528620


