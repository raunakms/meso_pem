### Load LIBRARIES ----
library("stringr")

### DEFINE PATH ----
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/04_driver_stats")
dir.data <- file.path(dir.wrk, "data")
dir.output <- file.path(dir.wrk, "output")
dir.plot <- file.path(dir.wrk, "plot")

### DEFINE FILES -----
file.driver <- file.path(dir.data, "mesoCombo_patientdrivers.txt")

### GET DRIVER DATA ---
dat.driver <- read.delim(file.driver, header=T, stringsAsFactors=F)
dat.driver <- dat.driver[which(str_detect(dat.driver$patient, "M") == TRUE),]
dat.driver <- subset(dat.driver, dat.driver$gamma == "0.8")

### PREPARE DATA ---
dat <- data.frame(SampleID=dat.driver$patient)
dat$Gamma <- dat.driver$gamma
dat$Alpha <- dat.driver$alpha
dat$Beta <- dat.driver$beta
dat$Genes <- dat.driver$driver.genes

dat$SampleID <- str_replace(dat$SampleID, "M", "MESO-")
dat$Genes <- str_replace_all(dat$Genes, ",", ":")

### WRITE OUTPUT ---
file.output <- file.path(dir.output, "meso_driver_genes.tsv")
write.table(dat, file.output, sep="\t", row.names=F, col.names=T, quote=F)
