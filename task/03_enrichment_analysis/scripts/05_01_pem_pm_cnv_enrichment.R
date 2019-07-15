### SET PATHS ----------------------------------------------------------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/03_enrichment_analysis")
dir.data <- file.path(dir.wrk, "data")
dir.enrichment <- file.path(dir.wrk, "enrichment")


### DEFINE FUNCTION ---
file.func <- file.path("/Data/Raunak/softwares/bdvtools/genesets_enrichment/bin/enrichment_test.R")
source(file.func)

### DEFINE FILE ---
file.dat <- file.path(dir.data, "cnv_meso_pem_pm_diffexpr_features.txt")
file.bg <- file.path(dir.data, "cnv_meso_pem_pm_background_genelist.txt")

### CALL FUNCTION ---
get.enrichment.test(genelist=file.dat, output=dir.enrichment, 
				batch="cnv_meso_pem_pm", db="MSIGDB", 
				background="CUSTOM", file.background=file.bg, 
				threshold=0.001)
