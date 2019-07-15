### DEFINE LIBRARIES ---
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/home/collins/IonProton/Mesothelioma")
dir.vcf <- file.path(dir.wrk, "vcf")
dir.annovar <- file.path(dir.wrk, "annovar")
dir.script <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/mutation/scripts")

### GET FILES ---
files.input <- list.files(dir.annovar, pattern=".annovarinput", full.names=FALSE)
files.output <- str_replace_all(files.input, ".annovarinput", ".TableAnnovar")

### PREPARE ANNOVAR SCRIPT ----
list.cmd <- list()
for(i in 1:length(files.input)){
	file.input <- file.path(dir.annovar, files.input[i])
	file.output <- file.path(dir.annovar, files.output[i])

	file.annovar <- file.path("/mnt/enclosure/mofan/software/ANNOVAR/annovar_2017Jul16/table_annovar.pl")
	dir.db <- file.path("/mnt/enclosure/mofan/software/ANNOVAR/annovar_2017Jul16/humandb")

	cmd <- paste(file.annovar, 
				file.input,
				dir.db,
				"-buildver hg19",
				"-protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a",
				"-operation g,r,f,f,f",
				"-otherinfo",
				"-nastring .",
				"-outfile",
				file.output,
				sep=" ")

	list.cmd[[i]] <- cmd
}

### compile results ---
cmd <- unlist(list.cmd)
cmd <- c("#!/bin/sh", cmd)

### WRITE OUTPUT ---
file.script <- file.path(dir.script, "11_03_tableannovar.sh")
write.table(cmd, file.script, row.names=F, col.names=F, quote=F)
