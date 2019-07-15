#### DEFINE PATH ------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/cnv")
dir.data <- file.path(dir.wrk, "seq_call_refseq_genes_meso")

### Define Files ------------------------------------------------------------
file.dat <- file.path(dir.data, "meso_seq_calls_matrix.txt")
file.header <- file.path(dir.data, "file_header.txt")

### Load Libraries ----------------------------------------------------------
library("stringr")

### Load Header Data --------------------------------------------------------
h <- read.table(file.header, header=F, stringsAsFactors=F)
index.prb <- which(h$V1 == "Probe")
index.seg <- which(h$V1 == "Segment")

h.prb <- h[index.prb,]
h.seg <- h[index.seg,]

genes <- unlist(lapply(str_split(h.prb$V4, "[(]"), function(x) x[1]))
chr <- unlist(lapply(str_split(unlist(lapply(str_split(h.prb$V4, "[(]"), function(x) x[2])), ":"), function(x) x[1]))
pos <- unlist(lapply(str_split(unlist(lapply(str_split(h.prb$V4, "[(]"), function(x) x[2])), ":"), function(x) x[2]))
pos.start <- as.numeric(str_replace_all(unlist(lapply(str_split(pos, "-"), function(x) x[1])),",",""))
pos.end <- as.numeric(str_replace_all(str_replace_all(unlist(lapply(str_split(pos, "-"), function(x) x[2])),",",""),"[)]",""))

df <- data.frame(Gene=genes, Chr=chr, Start=pos.start, End=pos.end)
write.table(df, file.path(dir.data, "refseq_gene_table.tsv"), sep="\t", row.names=F, col.names=T, quote=F)

### Cut Column Data from File --------------------------------------------
y.prb <- index.prb + 1
y.seg <- index.seg + 1

file.target.prb <- file.path(dir.data, "meso_cnv_probe_medians_calls.tsv")
file.target.seg <- file.path(dir.data, "meso_cnv_seg_values_calls.tsv")

cmd1 <- paste("cut -f1,1,", paste(paste(y.prb, collapse=","), " "), file.dat, " > ", file.target.prb, sep="")
cmd2 <- paste("cut -f1,1,", paste(paste(y.seg, collapse=","), " "), file.dat, " > ", file.target.seg, sep="")

cmd1 <- paste("cut -d, -f$")

write.table(cmd1, file.path(dir.data, "file_cmd_1.txt"), row.names=F, col.names=F, quote=F)

system(cmd1)
system(cmd2)

### Load CNA Data --------------------------------------------------------
dat <- read.delim(file.dat, header=T, stringsAsFactors=F, row.names=1, nrow=4)
classes <- sapply(dat, class)
dat <- read.delim(file.dat, header=T, stringsAsFactors=F, row.names=1, colClasses=classes)


