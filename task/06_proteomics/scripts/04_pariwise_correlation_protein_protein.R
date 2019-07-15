#### LOAD LIBRARIES -------------------
library("stringr")
library("RColorBrewer")
library("reshape2")
library("ggplot2")
library("gridExtra")

#### DEFINE PATH ----------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.proteome <- file.path(dir.wrk, "data/proteome/processed_data")
dir.task <- file.path(dir.wrk, "task/06_proteomics")
dir.data <- file.path(dir.task, "data")
dir.output <- file.path(dir.task, "output")
dir.plot <- file.path(dir.task, "plot")

#### DEFINE FILE ---------------------
file.expr <- file.path(dir.proteome, "meso_proteome_proteomediscover_all_log2_dqnorm.tsv.gz")
file.des <- file.path(dir.proteome, "design_table_proteome.tsv")

#### LOAD DESIGN TABLE ----------------
des <- read.delim(file.des, header=T, stringsAsFactors=F)
ids.normal <- des$SequencingID[which(des$SampleType == "Normal")]
ids.tumor <- des$SequencingID[which(des$SampleType == "Tumor")]
ids.celline <- des$SequencingID[which(des$SampleType == "Celline")]
ids <- c(ids.normal, ids.tumor, ids.celline)

#### LOAD PROTEIN EXPRESSION DATA -----
expr <- read.delim(file.expr, header=T, stringsAsFactors=F, row.names=1)
colnames(expr) <- str_replace_all(colnames(expr), "[.]", "-")
expr <- subset(expr, select=ids)
expr[is.na(expr)] <- 0

### GET IDS ---
ids.test <- subset(colnames(expr),  !(colnames(expr) %in% c("MESO-18AT","MESO-18ET")))

### CORRELATION BY SAMPLEPAIR ----
getCorrelation <- function(expr, ids.test, id.ref, file.plot){
	list.plot <- list()
	for(i in 1:length(ids.test)){
		sampleids <- c(id.ref, ids.test[i])
		expr.new <- subset(expr, select=sampleids)
		df <- data.frame(Gene=rownames(expr.new),
					Expr1=expr.new[,1],
					Expr2=expr.new[,2])
		df$Gene <- as.character(df$Gene)
		df <- df[-which((df$Expr1 == 0) | (df$Expr2 == 0)),]

		list.plot[[i]] <- getCorrelationPlot(df, sampleids)
	}

	### PLOT ---
	#file.plot <- file.path(dir.plot, "correlation_marray_rnaseq.pdf")
	pdf(file.plot, height=4, width=4)
		grid.arrange(list.plot[[1]], list.plot[[2]], list.plot[[3]], list.plot[[4]], ncol=2, nrow=2)
		grid.arrange(list.plot[[5]], list.plot[[6]], list.plot[[7]], list.plot[[8]], ncol=2, nrow=2)
		grid.arrange(list.plot[[9]], list.plot[[10]], list.plot[[11]], list.plot[[12]], ncol=2, nrow=2)
		grid.arrange(list.plot[[13]], list.plot[[14]], list.plot[[15]], list.plot[[16]], ncol=2, nrow=2)
		grid.arrange(list.plot[[17]], list.plot[[18]], list.plot[[19]], list.plot[[20]], ncol=2, nrow=2)
		grid.arrange(list.plot[[21]], list.plot[[22]], list.plot[[23]], list.plot[[24]], ncol=2, nrow=2)
	dev.off()
}

### FUNCTION: getCorrelationPlot ---
getCorrelationPlot <- function(df, sampleids){
	# LOAD LIBRARIES ---
	require("ggplot2")

	#GENERATE PLOT -----
	p <- ggplot(df, aes(x=Expr1, y=Expr2)) + 
			geom_point(color="black", alpha=0.5, size=0.1) +
			#scale_color_manual(values=c("#bdbdbd","#bd0026")) +
			#geom_text(aes(x=Expr.rnaseq, y=Expr.marray, label=Label), size=1, color="black", hjust=0, vjust=0)+
			#stat_smooth(method="lm", geom = "smooth", formula = y ~ x, position = "identity", fullrange = FALSE, se=FALSE) +
			geom_density2d(stat = "density2d", lineend = "round", linejoin = "round", alpha=0.9, color="yellow", size=0.25) +
			#coord_cartesian(xlim=c(0, 25), ylim=c(0, 25)) +
			theme(
				axis.text = element_text(size = 5, color="black"),
				axis.title = element_text(size = 8, color="black"),
				#strip.text = element_text(size = 10, color="black"),
				#strip.background = element_rect(fill = "white", colour = "white"),
				plot.title = element_text(size = 5, color="black", hjust=0),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.ticks = element_line(size=0.4, color = "black"),	
				panel.background = element_rect(fill = "white", colour = "black"),
				legend.position="none") +
			ylab(sampleids[2]) + 
			xlab(sampleids[1]) + 
			ggtitle("") 
	
	return(p)			
}


### PLOT ---
file.plot <- file.path(dir.plot, "correlation_MESO-18A_proteome.pdf")
getCorrelation(expr, ids.test, id.ref="MESO-18AT", file.plot)

file.plot <- file.path(dir.plot, "correlation_MESO-18E_proteome.pdf")
getCorrelation(expr, ids.test, id.ref="MESO-18ET", file.plot)
