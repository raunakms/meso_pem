#### LOAD LIBRARIES -------------------
library("stringr")

#### DEFINE PATH ----------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.proteome <- file.path(dir.wrk, "data/proteome/processed_data")
dir.task <- file.path(dir.wrk, "task/06_proteomics")
dir.data <- file.path(dir.task, "data")
dir.output <- file.path(dir.task, "output")
dir.plot <- file.path(dir.task, "plot")
dir.script <- file.path(dir.task, "scripts")
dir.enrichment <- file.path(dir.task, "enrichment")

#### DEFINE FILE ---------------------
file.expr <- file.path(dir.proteome, "meso_proteome_proteomediscover_all_log2_dqnorm_znorm.tsv.gz")
file.des <- file.path(dir.proteome, "design_table_proteome.tsv")
file.script <- file.path(dir.script, "function_get_pathway_zscore.R")

#### LOAD DESIGN TABLE ----------------
des <- read.delim(file.des, header=T, stringsAsFactors=F)
ids.normal <- des$SequencingID[which(des$SampleType == "Normal")]
ids.tumor <- des$SequencingID[which(des$SampleType == "Tumor")]
ids.celline <- des$SequencingID[which(des$SampleType == "Celline")]
ids <- c(ids.normal, ids.tumor, ids.celline)
group.names <- c(rep("Normal", length(ids.normal)), rep("Tumor",length(ids.tumor)), rep("Celline",length(ids.celline)))

#### LOAD PROTEIN EXPRESSION DATA -----
expr <- read.delim(file.expr, header=T, stringsAsFactors=F, row.names=1)
colnames(expr) <- str_replace_all(colnames(expr), "[.]", "-")
expr <- subset(expr, select=ids)
expr[is.na(expr)] <- 0

### Get Pathway Score --------

# REACTOME
source(file.script)
file.pathway <- file.path(dir.enrichment, "PROTEOME_DEG/enrichment_c2.cp.reactome.v6.0.symbols.txt")
file.plot <- file.path(dir.plot, "pathway_zscore_heatmap_reactome.pdf")
get.pathway.heatmap(expr, file.pathway=file.pathway, database="REACTOME", 
					group.names=group.names, 
					file.plot=file.plot, plot.height=15, plot.width=7)


# KEGG
source(file.script)
file.pathway <- file.path(dir.enrichment, "PROTEOME_DEG/enrichment_c2.cp.kegg.v6.0.symbols.txt")
file.plot <- file.path(dir.plot, "pathway_zscore_heatmap_kegg.pdf")
get.pathway.heatmap(expr, file.pathway=file.pathway, database="KEGG", 
					group.names=group.names, 
					file.plot=file.plot, plot.height=5, plot.width=7)


# GO:BP
source(file.script)
file.pathway <- file.path(dir.enrichment, "PROTEOME_DEG/enrichment_c5.bp.v6.0.symbols.txt")
file.plot <- file.path(dir.plot, "pathway_zscore_heatmap_gobp.pdf")
get.pathway.heatmap(expr, file.pathway=file.pathway, database="GOBP", 
					group.names=group.names, 
					file.plot=file.plot, plot.height=100, plot.width=7)



# HALLMARK
source(file.script)
file.pathway <- file.path(dir.enrichment, "PROTEOME_DEG/enrichment_h.all.v6.0.symbols.txt")
file.plot <- file.path(dir.plot, "pathway_zscore_heatmap_hallmark.pdf")
get.pathway.heatmap(expr, file.pathway=file.pathway, database="HALLMARK", 
					group.names=group.names, 
					file.plot=file.plot, plot.height=5, plot.width=7)





### REACTOME: NEW ----
file.dat <- file.path(dir.enrichment, "PROTEOME_DEG/enrichment_c2.cp.reactome.v6.0.symbols.txt")
dat <- read.delim(file.dat, header=TRUE, stringsAsFactors=FALSE)

file.pathways <- file.path(dir.output, "diff_pathway_reactome.tsv")
pathways <- read.table(file.pathways, header=FALSE, stringsAsFactors=FALSE)$V1

dat <- subset(dat, dat$Category %in% pathways)
write.table(dat, file.path(dir.enrichment, "PROTEOME_DEG/enrichment_reactome_curated.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# REACTOME
source(file.script)
file.pathway <- file.path(dir.enrichment, "PROTEOME_DEG/enrichment_reactome_curated.txt")
file.plot <- file.path(dir.plot, "pathway_zscore_heatmap_reactome_curated.pdf")
get.pathway.heatmap(expr, file.pathway=file.pathway, database="REACTOME", 
					group.names=group.names, 
					file.plot=file.plot, plot.height=5, plot.width=7)

