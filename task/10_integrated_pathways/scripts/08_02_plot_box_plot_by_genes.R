### DEFINE LIBRARIES ---
library("gridExtra")


### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.task <- file.path(dir.wrk, "task/10_integrated_pathways")
dir.data <- file.path(dir.task, "data")
dir.analysis <- file.path(dir.task, "output")
dir.plot <- file.path(dir.task, "plots")
dir.script <- file.path(dir.task, "scripts")
dir.cna <- file.path(dir.wrk, "data/cnv/seq_call_refseq_genes_meso")
dir.mrna <- file.path(dir.wrk, "data/expression/analysis")
dir.prot <- file.path(dir.wrk, "data/proteome/processed_data")

### DEFINE FILES ---
file.des.mrna <- file.path(dir.wrk, "data/annotation/design_table_3p21genes.tsv")
file.des.prot <- file.path(dir.wrk, "data/proteome/processed_data/design_table_proteome.tsv")
file.cna <- file.path(dir.cna, "meso_cnv_seg_values_calls_parsed.tsv")
file.mrna <- file.path(dir.mrna, "meso_peritoneal_gene_expression_invnorm_proteincoding.tsv.gz")
file.prot <- file.path(dir.prot, "meso_proteome_proteomediscover_all_log2_dqnorm_znorm.tsv.gz")

### LOAD FUNCTIONS ---
file.script <- file.path(dir.script, "08_01_function_plot_boxplot.R")
source(file.script)

#### HDAC + HDAC Complex ---
#genes <- c("HDAC1","HDAC2","HDAC3","HDAC4","HDAC5","HDAC6",
#			"SIRT4","SIRT2","KDM1A","MTA2","GTF2I","RBBP4","CHD3","MTA1",
#			"ZMYM3","GSE1","CHD4","RBBP7","PHF21A","SIN3A","HMG20B","ZMYM2","RCOR1")

genes <- c("HDAC1","HDAC2","HDAC4","HDAC6","CHD4","ZMYM2","ZMYM3")
			
file.plot <- file.path(dir.plot, "boxplot_HDAC.pdf")
pdf(file.plot, height=3, width=3)
	#p1 <- getProfile(dir.plot, file.cna, file.mrna, file.prot, genes, type="CNA")
	p2 <- getProfile(dir.plot, file.cna, file.mrna, file.prot, genes, type="RNA")
	p3 <- getProfile(dir.plot, file.cna, file.mrna, file.prot, genes, type="PROTEIN")
	grid.arrange(p2, p3, ncol=1, nrow=2)
dev.off()




### SWI/SNF
#genes <- c("ARID1A","ARID1B","ARID2","ACTB","ACTG1","ACTL6A","ACTL6B","PBRM1",
#			"SMARCA4","SMARCB1","SMARCC1","SMARCC2","SMARCD1","SMARCE1","SMARCA2")

genes <- c("ARID2","ACTL6A","PBRM1",
			"SMARCA4","SMARCB1","SMARCC1","SMARCC2")

file.plot <- file.path(dir.plot, "boxplot_SWISNF.pdf")
pdf(file.plot, height=3, width=3)
	#p1 <- getProfile(dir.plot, file.cna, file.mrna, file.prot, genes, type="CNA")
	p2 <- getProfile(dir.plot, file.cna, file.mrna, file.prot, genes, type="RNA")
	p3 <- getProfile(dir.plot, file.cna, file.mrna, file.prot, genes, type="PROTEIN")
	grid.arrange(p2, p3, ncol=1, nrow=2)
dev.off()


