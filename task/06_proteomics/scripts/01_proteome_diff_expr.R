#### LOAD LIBRARIES -------------------
library("stringr")
library("gplots")
library("RColorBrewer")

#### DEFINE PATH ----------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.proteome <- file.path(dir.wrk, "data/proteome/processed_data")
dir.task <- file.path(dir.wrk, "task/06_proteomics")
dir.data <- file.path(dir.task, "data")
dir.output <- file.path(dir.task, "output")
dir.plot <- file.path(dir.task, "plot")
dir.enrichment <- file.path(dir.task, "enrichment")

#### DEFINE FILE ---------------------
file.expr <- file.path(dir.proteome, "meso_proteome_proteomediscover_all_log2_dqnorm.tsv.gz")
file.exprz <- file.path(dir.proteome, "meso_proteome_proteomediscover_all_log2_dqnorm_znorm.tsv.gz")
file.des <- file.path(dir.proteome, "design_table_proteome.tsv")
file.script <- file.path("/Data/Raunak/softwares/bdvtools/array_process/diff_expr.R")

#### LOAD DESIGN TABLE ----------------
des <- read.delim(file.des, header=T, stringsAsFactors=F)
ids.normal <- des$SequencingID[which(des$SampleType == "Normal")]
ids.tumor <- des$SequencingID[which(des$SampleType == "Tumor")]
ids.celline <- des$SequencingID[which(des$SampleType == "Celline")]

#### LOAD PROTEIN EXPRESSION DATA -----
expr <- read.delim(file.expr, header=T, stringsAsFactors=F, row.names=1)
colnames(expr) <- str_replace_all(colnames(expr), "[.]", "-")

write.table(rownames(expr), file.path(dir.data, "meso_proteome_expression_background_genelist.txt"), row.names=F, col.names=F, quote=F)

#### DIFFERENTIAL EXPRESSION ANALYSIS ----
source(file.script)
df <- get.diffexpr(expr, class1=ids.normal, class2=c(ids.tumor, ids.celline), var.equal=FALSE)

df <- subset(df, df$fdr <= 0.1)
genes.deg <- df$Gene

#### WRITE OUTPUT -----------
file.output <- file.path(dir.output, "meso_proteome_deg_summary_table.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)

file.deg <- file.path(dir.output, "meso_proteome_deg_genelist.txt")
write.table(genes.deg , file.deg, sep="\t", row.names=F, col.names=F, quote=F)

#### PLOT HEATMAP -----------
exprz <- read.delim(file.exprz, header=T, stringsAsFactors=F, row.names=1)
colnames(exprz) <- str_replace_all(colnames(exprz), "[.]", "-")

expr.deg <- subset(exprz, rownames(exprz) %in% genes.deg)
expr.deg <- subset(expr.deg, select=c(ids.normal, ids.tumor, ids.celline))

expr.deg[is.na(expr.deg)] <- 0

#### GROUP DATA LABEL ------------
label.groups <- c(rep("#4ECDC4", length(ids.normal)), rep("#C44D58",length(ids.tumor)),  rep("#fdae6b",length(ids.celline)))

# Generate Color Palette ------
color_values <- c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#ffffff","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695") #red-blue
color_palette <- colorRampPalette(color_values)(11)
jColFun <- colorRampPalette(color_palette)

# Generate Plot ------
file.plot <- file.path(dir.plot, "heatmap_proteome_expression_deg.pdf")
pdf(file.plot, height=4, width=10)
	heatmap.2(t(as.matrix(expr.deg)), 
		col = rev(jColFun(256)),
		Colv=TRUE, Rowv=TRUE, 
		dendrogram ="both", trace="none",  scale="none",
		cexCol=0.3, cexRow=0.5, margins = c(10,5), 
		hclustfun = function(x) hclust(x, method = "ward.D2"), 
		distfun = function(x) dist(x, method = "euclidean"),
		#colsep=c(1:1000), rowsep=c(1:1000),
		#sepcolor="black", sepwidth=c(0.000005,0.000005),
		RowSideColors=label.groups,
		key="TRUE", keysize=1, density.info="none", symkey=0,
		key.title=NA, # no title
		key.xlab=NA,  # no xlab
		key.par=list(mgp=c(1.5, 0.5, 0),
                       mar=c(2.5, 2.5, 1, 0)))
	legend("bottomleft", legend=c("Normal","Tumor","Celline"), fill=c("#4ECDC4","#C44D58","#fdae6b"), border=FALSE, bty="n", x.intersp = 1, y.intersp = 1, cex=0.5)
dev.off()


#### PATHWAY ENRICHMENT ----------------
file.script <- file.path("/Data/Raunak/softwares/bdvtools/genesets_enrichment/bin/enrichment_test.R")
source(file.script)

file.genelist <- file.path(dir.output, "meso_proteome_deg_genelist.txt")
file.bg <- file.path(dir.data, "meso_proteome_expression_background_genelist.txt")
get.enrichment.test(genelist=file.genelist, output=dir.enrichment, batch="PROTEOME_DEG",
					db="MSIGDB", background="CUSTOM", file.background=file.bg, threshold=0.01)
