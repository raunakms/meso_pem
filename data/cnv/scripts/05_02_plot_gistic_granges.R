### LOAD LIBRARIES ---
library("stringr")
library("ggplot2")
library("gridExtra")
library("ggbio")
library("GenomicRanges")

### DEFINE PATH ---
dir.wrk <- file.path("//jbrcsrv009/CollinsGroup/Raunak/projects/MESO_peritoneal/data/cnv")
#dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/cnv")
dir.data <- file.path(dir.wrk, "main_calls")
dir.plot <- file.path(dir.wrk, "plot")

### DEFINE FILES ---
file.dat <- file.path(dir.data, "MESO-Pem_gistic_cna_calls.txt")

### LOAD FILES ---
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)

### ADD ATTRIBUTES ---
dat$Region <- str_replace_all(dat$Region, ",", "")
dat$Chr <- str_replace(unlist(lapply(str_split(dat$Region, ":"), function(x) x[1])), "chr", "")
dat$Start <- unlist(lapply(str_split(unlist(lapply(str_split(dat$Region, ":"), function(x) x[2])), "-"), function(x) as.numeric(x[1])))
dat$End <- unlist(lapply(str_split(unlist(lapply(str_split(dat$Region, ":"), function(x) x[2])), "-"), function(x) as.numeric(x[2])))

### SUBSET DATA ---
df <- subset(dat, select=c("Chr","Start","End","Type","Q.Bound","G.Score"))
df$nlogp <- -log10(df$Q.Bound)

### SUBSET DATA BY TYPE ---
df.amp <- subset(df, df$Type == "CN Gain")
df.del <- subset(df, df$Type == "CN Loss")

### CONVERT TO GENOMICS RANGES ---
gr.amp <- GRanges(df.amp)
gr.del <- GRanges(df.del)

### FACTORIZE DATA ---
genome(gr.amp) <- "hg19"
gr.amp <- renameSeqlevels(gr.amp, c(1:22, "X", "Y"))

gr.amp <- keepSeqlevels(gr.amp, c(1:22, "X", "Y"))


### PLOT ---
autoplot(gr.amp, coord = "genome", geom = "line", aes(y = nlogp))






autoplot(gr.snp, coord = "genome", geom = "line", aes(y = pvalue, group = seqnames,
color = seqnames))