### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/04_driver_stats")
dir.output <- file.path(dir.wrk, "output")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")
dir.enrichment <- file.path(dir.wrk, "enrichment")

#### PATHWAY ENRICHMENT ----------------
file.script <- file.path("/Data/Raunak/softwares/bdvtools/genesets_enrichment/bin/enrichment_test.R")
source(file.script)

file.background <- file.path(dir.output, "genelist_background_cnv_mrna_protein.txt")

### ANALYSIS - 1
file.genelist <- file.path(dir.output, "meso_driver_genelist.txt")
get.enrichment.test(genelist=file.genelist, output=dir.enrichment, batch="driver_genes", db="MSIGDB", background="CUSTOM", file.background, threshold=0.01)
