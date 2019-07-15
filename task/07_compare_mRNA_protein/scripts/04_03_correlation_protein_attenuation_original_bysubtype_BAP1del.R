### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.analysis <- file.path(dir.wrk, "task/07_compare_mRNA_protein")
dir.output <- file.path(dir.analysis, "output")
dir.plot <- file.path(dir.analysis, "plot/protein_attenuation_correlation")
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
file.script <- file.path(dir.script, "04_00_function_protein_attenuation_ptm.R")

### CALL FUNCTION ---
source(file.script)

### GET DATA ---
dat.cna <- loadCNA(file.dat=file.cna, file.des=file.des.mrna)
expr.mrna <- loadmRNAExpr(file.dat=file.mrna, file.des=file.des.mrna)
expr.prot <- loadProteinExpr(file.dat=file.prot, file.des=file.des.prot)

### GET COMMON SAMPLEIDS ---
des <- getDesignTable(dat.cna, expr.mrna, expr.prot, file.des.mrna)
des <- subset(des, des$Group.3p21 == 1)
ids <- des$Sample.ID

### SUBSET AND ARRANGE DATA ---
dat.cna <- subset(dat.cna, select=ids)
expr.mrna <- subset(expr.mrna, select=ids)
expr.prot <- subset(expr.prot, select=ids)


### COMBINE DATA: CNA-mRNA --- 
df.cna_mrna <- prepareComboData(dat1=dat.cna, dat2=expr.mrna, des=des, type1="CNA", type2="mRNA")
df.cna_prot <- prepareComboData(dat1=dat.cna, dat2=expr.prot, des=des, type1="CNA", type2="PROTEIN")


### CALL FUNCTION: getCorrelation -----
dfcor.cna_mrna <- getCorrelation(dat=df.cna_mrna)
dfcor.cna_prot <- getCorrelation(dat=df.cna_prot)

### WRITE OUTPUT ---
file.output1 <- file.path(dir.output, "meso_correlation_original_BAP1DEL_cna_mrna.tsv")
file.output2 <- file.path(dir.output, "meso_correlation_original_BAP1DEL_cna_protein.tsv")

write.table(dfcor.cna_mrna, file.output1, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(dfcor.cna_prot, file.output2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

### RE-LOAD FILES ----
file.dat1 <- file.path(dir.output, "meso_correlation_original_BAP1DEL_cna_mrna.tsv")
file.dat2 <- file.path(dir.output, "meso_correlation_original_BAP1DEL_cna_protein.tsv")

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
				Type2=dat2$Type,
				R2=dat2$R)



df$Gene <- as.character(df$Gene)
df$Type1 <- as.character(df$Type1)
df$Type2 <- as.character(df$Type2)

### REMOVE GENES WITH NO DATA ---
ind.del <- which((df$R1 == 0) | (df$R2 == 0))
df <- df[-ind.del,]

#> dim(df)
#[1] 6383    7

### COMPUTE ABSOLUTE DIFFERENCE BETWEEN R1-R2
df$Rdiff <- apply(df, 1, function(x) as.numeric(x[3]) - as.numeric(x[5]))
df$Rdiff.group <- ifelse(df$Rdiff >= 0.45, 1, 0)

#> table(df$Rdiff.group)
#
#   0    1
#4748 1635

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
file.output <- file.path(dir.output, "meso_correlation_original_BAP1DEL_cna-mrna_cna-protein.tsv")
write.table(df, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


### SUBSET PROTEIN ATTENUATED GENES ---
d <- subset(df, df$Rdiff.group == 1)
file.output <- file.path(dir.output, "meso_protein_original_BAP1DEL_attenuated_genes.tsv")
write.table(d$Gene, file.output, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

### PLOT ---
file.plot <- file.path(dir.plot, "meso_correlation_original_BAP1DEL_cna-mrna_cna-protein.pdf")
pdf(file.plot, height=4, width=4)
	#grid.arrange(get.plot(df), ncol=1, nrow=1)
	grid.arrange(get.densityplot(df, type="R1"), 
				get.densityplot(df, type="R2"),
				ncol=1, nrow=3)
dev.off()

