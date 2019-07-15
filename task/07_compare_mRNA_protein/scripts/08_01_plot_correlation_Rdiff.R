### LOAD LIBRARIES -----
library("stringr")
library("ggplot2")
library("gridExtra")
library("igraph")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.analysis <- file.path(dir.wrk, "task/07_compare_mRNA_protein")
dir.output <- file.path(dir.analysis, "output")
dir.plot <- file.path(dir.analysis, "plot/protein_attenuation_correlation")
dir.script <- file.path(dir.analysis, "scripts")
dir.network <- file.path("/Data/Raunak/data_ref/interaction_networks/STRING10/data")

### DEFINE FILES ---
file.network <- file.path(dir.network, "string10_network_all.tsv.gz")
file.dat11 <- file.path(dir.output, "meso_correlation_original_BAP1DEL_cna-mrna_cna-protein.tsv")
file.dat12 <- file.path(dir.output, "meso_correlation_original_BAP1INTACT_cna-mrna_cna-protein.tsv")

file.dat21 <- file.path(dir.output, "meso_correlation_ptm_BAP1DEL_cna-mrna_cna-protein.tsv")
file.dat22 <- file.path(dir.output, "meso_correlation_ptm_BAP1INTACT_cna-mrna_cna-protein.tsv")

### FUNCTION: get.friends ---
get.friends <- function(file.network, node){
	# LOAD NETWORK ---
	dat.network <- read.delim(file.network, header=FALSE, stringsAsFactors=FALSE)

	# CREATE igraph OBJECT ---
	g <- graph.data.frame(dat.network, directed=FALSE)
	g <- simplify(graph=g, remove.multiple=TRUE, remove.loops=TRUE)

	# GET NETWORK NEIGHBOURS OF BAP1 ---
	nodes.friends <- ego(graph=g, order=1, nodes=which(V(g)$name %in% node), mode="all", mindist=0)[[1]]

	# GET SUBGRAPH ---
	sg.nodes.friends <- induced.subgraph(g, nodes.friends)

	return(V(sg.nodes.friends)$name)
}

### FUNCTION: prepareData ----
prepareData <- function(file.dat, nodes.friends, analysis, subtype){
	# LOAD DATA --
	dat <- read.delim(file.dat, header=TRUE, stringsAsFactors=FALSE)

	# TRIM DATA ---
	dat <- subset(dat, dat$Gene %in% nodes.friends)
	dat <- subset(dat, select=c("Gene", "R1", "R2", "Rdiff"))
	
	# ORDER DATA ---
	dat <- dat[order(dat$Gene, decreasing=FALSE),]

	# ADD ANNOTATION ---
	dat$Analysis <- analysis
	dat$Subtype <- subtype

	return(dat)
}

### FUNCTION: combineData ---
combineData <- function(df1, df2){
	# GET COMMON GENES ---
	genes <- intersect(df1$Gene, df2$Gene)

	# SUBSET DATA ---
	df1 <- subset(df1, df1$Gene %in% genes)
	df2 <- subset(df2, df2$Gene %in% genes)

	# COMBINE DATA ---
	df <- rbind(df1, df2)

	# FACTORIZE DATA ---
	df$Subtype <- factor(df$Subtype, levels=c("BAP1INT","BAP1DEL"))
	df$Analysis <- as.factor(df$Analysis)
	df$Gene <- as.factor(df$Gene)

	return(df)
}

### FUNCTION: getBoxPlot ---
getBoxPlot <- function(df){
	# COLOR ---
	cbPalette <- c("#fffc00","blue")

	# PREPARE PLOT ---
	p <- ggplot(df, aes(x=Subtype, y=Rdiff)) +
			geom_boxplot(aes(fill=Subtype), lwd=0.3, color="#000000", alpha=0.5, outlier.size = 0.1, outlier.alpha = 0.3, notch=FALSE) +
			scale_fill_manual(values=cbPalette) +
			geom_line(aes(group = Gene), alpha = 0.3, colour = "darkgrey") +
			#facet_wrap(~Analysis, nrow=1, ncol=2) +
			theme(
				axis.text.x = element_text(size = 7, color="#000000"),
				axis.text.y = element_text(size = 7, color="#000000"),
				axis.title = element_text(size = 7, color="#000000"),
				plot.title = element_text(size = 10, color="#000000", hjust=0),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),			
				axis.ticks = element_line(size=0.1, color="#000000"),
				strip.text = element_text(size=3, color="#000000"),
				strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
				panel.background = element_rect(fill = "#FFFFFF", color = "#000000"),
				legend.position="none") + 
			ylab("Rdiff") +
			xlab("") + 
			ggtitle("") 

	return(p)
}

### GET BAP1 NETIGHBOURHOOD ---
nodes.bap1 <- get.friends(file.network, node="BAP1")

### PREPARE DATA ---
dat11 <- prepareData(file.dat=file.dat11, nodes.friends=nodes.bap1, analysis="NONPTM", subtype="BAP1DEL")
dat12 <- prepareData(file.dat=file.dat12, nodes.friends=nodes.bap1, analysis="NONPTM", subtype="BAP1INT")
dat21 <- prepareData(file.dat=file.dat21, nodes.friends=nodes.bap1, analysis="PTM", subtype="BAP1DEL")
dat22 <- prepareData(file.dat=file.dat22, nodes.friends=nodes.bap1, analysis="PTM", subtype="BAP1INT")

### COMBINE DATA ---
df.nonptm <- combineData(df1=dat11, df2=dat12)
df.ptm <- combineData(df1=dat21, df2=dat22)

### GENERATE PLOT ---
p.nonptm <- getBoxPlot(df=df.nonptm)
p.ptm <- getBoxPlot(df=df.ptm)

### PLOT DATA ---
file.plot <- file.path(dir.plot, "boxplot_Rdiff.pdf")
pdf(file.plot, height=2, width=3)
		grid.arrange(p.nonptm, p.ptm, ncol=2, nrow=1)
dev.off()


# https://stackoverflow.com/questions/36240695/connect-ggplot-boxplots-using-lines-and-multiple-factor
