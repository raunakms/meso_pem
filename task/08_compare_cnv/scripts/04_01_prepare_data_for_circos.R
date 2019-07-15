### LOAD LIBRAIES ---
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/08_compare_cnv")
dir.data <- file.path(dir.wrk, "data")
dir.output <- file.path(dir.wrk, "output")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")

### DEFINE FILES ---
file.meso <- file.path(dir.data, "meso_peritoneal_nexus_aggregate_bedgraph.bgr")
file.ov <- file.path(dir.data, "TCGAOV_nexus_aggregate_bedgraph.bgr")

### LOAD MESO CNV ----
dat.meso <-  read.delim(file.meso, header=F, stringsAsFactors=F, skip=5)
dat.ov <-  read.delim(file.ov, header=F, stringsAsFactors=F, skip=5)


### REPLACE chr with hs ---
dat.meso$V1 <- str_replace(dat.meso$V1, "chr", "hs")
dat.ov$V1 <- str_replace(dat.ov$V1, "chr", "hs")

### REMOVE NA ---
dat.meso <- na.omit(dat.meso)
dat.ov <- na.omit(dat.ov)


### WRITE OUTPUT MESO ---
file.output <- file.path(dir.output, "cna_circos_meso.histo")
write.table(dat.meso, file.output, sep=" ", row.names=F, col.names=F, quote=F)

### WRITE OUTPUT OV ---
file.output <- file.path(dir.output, "cna_circos_ov.histo")
write.table(dat.ov, file.output, sep=" ", row.names=F, col.names=F, quote=F)

#> summary(dat.meso$V4)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-50.000 -11.000  -6.000  -1.768   6.000  28.000
#
#> summary(dat.ov$V4)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-99.000 -27.000   8.000   1.534  26.000  89.000
