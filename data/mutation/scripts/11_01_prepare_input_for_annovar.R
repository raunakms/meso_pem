### DEFINE LIBRARIES ---
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/home_oldcluster/collins/IonProton/Mesothelioma")
dir.vcf <- file.path(dir.wrk, "vcf")
dir.annovar <- file.path(dir.wrk, "annovar")
dir.script <- file.path("/collinsgroup/Raunak/projects/MESO_peritoneal/data/mutation/scripts")

### GET SAMPLEIDS ---
files <- list.files(dir.vcf, pattern="vcf", full.names=FALSE)
files <- str_replace(files, "LAGA-", "")
files <- str_replace(files, ".adj.vcf", "N")
files <- str_replace(files, ".vcf", "T")

### GET VCF FILES ---
files.vcf <- list.files(dir.vcf, pattern="vcf", full.names=FALSE)

#### LOOP BY SAMPLEIDS ----
list.cmd <- list()
for(i in 1:length(files.vcf)){
	sampleid <- files[i]

	# DEFINE FILES ----
	file.vcf <- file.path(dir.vcf, files.vcf[i])		
	file.input.annovar <- file.path(dir.annovar, paste(sampleid, ".annovarinput", sep=""))

	# PREPARE COMMANDS ---
	file.annovar <- "/mnt/enclosure/mofan/software/ANNOVAR/annovar_2017Jul16/convert2annovar.pl"
	cmd <- paste(file.annovar, "-format vcf4", file.vcf, "-outfile", file.input.annovar, "-include", sep=" ")

	list.cmd[[i]] <- cmd
}

### compile results ---
cmd <- unlist(list.cmd)
cmd <- c("#!/bin/sh", cmd)

### WRITE OUTPUT ---
file.output <- file.path(dir.script, "11_02_convert2annovar.sh")
write.table(cmd, file.output, row.names=F, col.names=F, quote=F)

