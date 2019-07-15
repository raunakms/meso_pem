### LOAD LIBRAIES ---
library("stringr")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/08_compare_cnv")
dir.data <- file.path(dir.wrk, "data")
dir.output <- file.path(dir.wrk, "output")
dir.result <- file.path(dir.wrk, "results")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")
dir.function <- file.path("/Data/Raunak/softwares/bdvtools/machine_learning")

### DEFINE FILES ---
file.dat <- file.path(dir.data, "cnv_meso_ov_seg_mean_combined.tsv")
file.cls <- file.path(dir.data, "cnv_meso_ov_combined_class.txt")
file.features <- file.path(dir.output, "cnv_meso_ov_diff_summary_table_by_chromosome.tsv")

cat(paste(Sys.time()), "LOADING DATA MATRIX ...", "\n", sep=" ")

### LOAD DATA ---
dat <-  read.delim(file.dat, header=T, stringsAsFactors=F, row.names=1, nrow=5)
classes <- sapply(dat, class)
dat <- read.delim(file.dat, header=T, stringsAsFactors=F,  row.names=1, colClasses=classes)
colnames(dat) <- str_replace_all(colnames(dat), "[.]", "-")

cat(paste(Sys.time()), "ALL DONE ...", "\n", sep=" ")

### LOAD CLASS ---
cls <- read.table(file.cls, header=F, stringsAsFactors=F)$V1

### LOAD FEATURES ---
dat.features <- read.table(file.features, header=T, stringsAsFactors=F)

### SOURCE FUNCTIONS ---
source(file.path(dir.function, "run_kfold_validation.R"))

### CALL FUNCTION FOR k-FOLD VALIDATION ---
run.kfold.validation(univeral.matrix=dat, cls=cls, features=dat.features$Gene, nFolds=5, nReps=100,
						batch.name="CNV_MESO_OV_WCOX", dir.output=dir.result)


