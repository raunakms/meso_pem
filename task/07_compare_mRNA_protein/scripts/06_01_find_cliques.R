### Load Libraries -----
library("stringr")
library("ggplot2")
library("gridExtra")
library("igraph")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.analysis <- file.path(dir.wrk, "task/07_compare_mRNA_protein")
dir.output <- file.path(dir.analysis, "output")
dir.plot <- file.path(dir.analysis, "plot")
dir.script <- file.path(dir.analysis, "scripts")
dir.prot <- file.path(dir.wrk, "data/proteome/processed_data")
dir.network <- file.path("/Data/Raunak/data_ref/interaction_networks/STRING10/data")

### DEFINE FILES ---
file.network <- file.path(dir.network, "string10_network.tsv")

### LOAD NETWORK ---
dat.network <- read.delim(file.network, header=F, stringsAsFactors=F)

### Create iGraph object ---
g <- graph.data.frame(dat.network, directed=FALSE)
g <- simplify(graph=g, remove.multiple=TRUE, remove.loops=TRUE)

### FIND CLIQUES ---
list.cq <- cliques(graph=g, min=10, max=30)
#vids <- which(V(g)$name %in% c( "SETD2","PBRM1","SMARCC1","BAP1"))
#list.cq <- max_cliques(graph=g, min=10, max=30, subset=10971, file=NULL)

