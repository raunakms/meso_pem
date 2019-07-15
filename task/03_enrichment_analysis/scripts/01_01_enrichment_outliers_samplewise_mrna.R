### SET PATHS ----------------------------------------------------------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/03_enrichment_analysis")
dir.data <- file.path(dir.wrk, "data")
dir.enrichment <- file.path(dir.wrk, "enrichment")

### DEFINE FILES -------------------------------------------------------------------------------------------
file.script <- file.path("/Data/Raunak/softwares/bdvtools/genesets_enrichment/bin/enrichment_test.R")
file.genelist <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/09_compare_expression/output/meso_outliers_samplewise_genelist.tsv")
source(file.script)


### ENRICHMENT TEST ----------------------------------------------------------------------------------------------
get.enrichment.test.module(genelist=file.genelist, output=dir.enrichment, batch="MESO_OUTLIERS_mRNA_p01",
					db="MSIGDB", background="BG_RNASEQ_MESO", file.background=NA, threshold=0.05, col.genes=2)
