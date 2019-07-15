### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/07_compare_mRNA_protein")
dir.output <- file.path(dir.wrk, "output")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")
dir.enrichment <- file.path(dir.wrk, "enrichment")

#### PATHWAY ENRICHMENT ----------------
file.script <- file.path("/Data/Raunak/softwares/bdvtools/genesets_enrichment/bin/enrichment_test.R")
source(file.script)

file.background <- file.path(dir.output, "genelist_background_cnv_mrna_protein.txt")

### ANALYSIS - 1 ---
file.genelist <- file.path(dir.output, "genelist_positive_mRNA-Protein.txt")
get.enrichment.test(genelist=file.genelist, output=dir.enrichment, batch="genelist_positive_mRNA-Protein", db="MSIGDB", background="CUSTOM", file.background, threshold=0.001)

### ANALYSIS - 2 ---
file.genelist <- file.path(dir.output, "genelist_positive_mRNA-CNV.txt")
get.enrichment.test(genelist=file.genelist, output=dir.enrichment, batch="genelist_positive_mRNA-CNV", db="MSIGDB", background="CUSTOM", file.background, threshold=0.001)

### ANALYSIS - 3 ---
file.genelist <- file.path(dir.output, "genelist_positive_common_mRNA-CNV_mRNA-Protein.txt")
get.enrichment.test(genelist=file.genelist, output=dir.enrichment, batch="genelist_positive_common_mRNA-CNV_mRNA-Protein", db="MSIGDB", background="CUSTOM", file.background, threshold=0.001)


### ANALYSIS - 4 ---
file.genelist <- file.path(dir.output, "correlation_summary_mrna_protein_pearson_negative_genelist.txt")
get.enrichment.test(genelist=file.genelist, output=dir.enrichment, batch="genelist_negative_mRNA-Protein", db="MSIGDB", background="CUSTOM", file.background, threshold=0.001)



### ANALYSIS - 5 ---
file.genelist <- file.path(dir.output, "meso_protein_attenuated_genes.tsv")
get.enrichment.test(genelist=file.genelist, output=dir.enrichment, batch="genelist_protein_attenuation", db="MSIGDB", background="CUSTOM", file.background, threshold=0.0001)
