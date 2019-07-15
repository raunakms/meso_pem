### LOAD LIBRARIES ---
library("RColorBrewer")
library("ggplot2")
library("gridExtra")
library("reshape2")

### DEFINE PATH -----
dir.wrk <- file.path("//jbrcsrv009/CollinsGroup/Raunak/projects/MESO_peritoneal")
#dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal")
dir.task <- file.path(dir.wrk, "task/10_integrated_pathways")
dir.data <- file.path(dir.task, "data")
dir.analysis <- file.path(dir.task, "output")
dir.plot <- file.path(dir.task, "plots")

### DEFINE FILES ---
file.dat <- file.path(dir.data, "sample_annotation.tsv")
file.samples <- file.path(dir.data, "order_samples.txt")

### LOAD DATA ---
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
ids <- read.table(file.samples, header=F, stringsAsFactors=F)$V1

### FACTORIZE DATA ---
dat$CancerType <- as.factor(dat$CancerType)
dat$SampleID <- factor(dat$SampleID, levels=ids)
dat$Gender <- factor(dat$Gender, levels=c("M","F"))
dat$Group.3p21 <- factor(dat$Group.3p21, levels=c("G0","G1"))
dat$Group.CNA <- factor(dat$Group.CNA, levels=c("G1","G2","G3","G4"))
dat$Group.Expr.mrna <- factor(dat$Group.Expr.mrna, levels=c("G1","G2"))
dat$Group.Expr.protein <- factor(dat$Group.Expr.protein, levels=c("G1","G2","G3","G4"))

### PLOT: Sex ---
p1 <- ggplot(dat, aes(x=SampleID, y=CancerType)) + 
		geom_tile(aes(fill=Gender), color="black", size=0.3) +
		scale_fill_manual(values=c("#2171b5","#e7298a")) +
		theme(
			axis.text = element_blank(),
			#axis.text.x = element_text(size = 3, color="black", angle=90, hjust=1, vjust=0.5),
			axis.title = element_blank(),
			plot.title = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.ticks.x =element_blank(),
			axis.ticks.y =element_blank(),
			panel.background = element_rect(fill="white"),
			legend.text = element_text(size = 3, color="black"),
			legend.title = element_text(size = 3, color="black"),
			legend.key.size = unit(0.2, "cm"),
			legend.position="none") +
		xlab("") + ylab("") + ggtitle("") 

### PLOT: SNV.burden ---
p2 <- ggplot(dat, aes(x=SampleID, y=CancerType)) + 
		geom_tile(aes(fill=SNV.burden), color="black", size=0.3) +
		scale_fill_gradient(low="#FFCCD1", high="#D50015") +
		theme(
			axis.text = element_blank(),
			#axis.text.x = element_text(size = 3, color="black", angle=90, hjust=1, vjust=0.5),
			axis.title = element_blank(),
			plot.title = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.ticks.x =element_blank(),
			axis.ticks.y =element_blank(),
			panel.background = element_rect(fill="white"),
			legend.text = element_text(size = 3, color="black"),
			legend.title = element_text(size = 3, color="black"),
			legend.key.size = unit(0.2, "cm"),
			legend.position="none") +
		xlab("") + ylab("") + ggtitle("") 

### PLOT: CNV.burden ---
p3 <- ggplot(dat, aes(x=SampleID, y=CancerType)) + 
		geom_tile(aes(fill=CNV.burden), color="black", size=0.3) +
		scale_fill_gradient(low="#FFCCD1", high="#D50015") +
		theme(
			axis.text = element_blank(),
			#axis.text.x = element_text(size = 3, color="black", angle=90, hjust=1, vjust=0.5),
			axis.title = element_blank(),
			plot.title = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.ticks.x =element_blank(),
			axis.ticks.y =element_blank(),
			panel.background = element_rect(fill="white"),
			legend.text = element_text(size = 3, color="black"),
			legend.title = element_text(size = 3, color="black"),
			legend.key.size = unit(0.2, "cm"),
			legend.position="none") +
		xlab("") + ylab("") + ggtitle("") 

### PLOT: CNA-GROUP ---
p4 <- ggplot(dat, aes(x=SampleID, y=CancerType)) + 
		geom_tile(aes(fill=Group.CNA), color="black", size=0.3) +
		scale_fill_manual(values=c("#85d988","#c18377","#d3cdab", "#6f6d76")) +
		theme(
			axis.text = element_blank(),
			#axis.text.x = element_text(size = 3, color="black", angle=90, hjust=1, vjust=0.5),
			axis.title = element_blank(),
			plot.title = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.ticks.x =element_blank(),
			axis.ticks.y =element_blank(),
			panel.background = element_rect(fill="white"),
			legend.text = element_text(size = 3, color="black"),
			legend.title = element_text(size = 3, color="black"),
			legend.key.size = unit(0.2, "cm"),
			legend.position="none") +
		xlab("") + ylab("") + ggtitle("") 

### PLOT: EXPR-mRNA-GROUP ---
p5 <- ggplot(dat, aes(x=SampleID, y=CancerType)) + 
		geom_tile(aes(fill=Group.Expr.mrna), color="black", size=0.3) +
		scale_fill_manual(values=c("#c18377","#85d988")) +
		theme(
			axis.text = element_blank(),
			#axis.text.x = element_text(size = 3, color="black", angle=90, hjust=1, vjust=0.5),
			axis.title = element_blank(),
			plot.title = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.ticks.x =element_blank(),
			axis.ticks.y =element_blank(),
			panel.background = element_rect(fill="white"),
			legend.text = element_text(size = 3, color="black"),
			legend.title = element_text(size = 3, color="black"),
			legend.key.size = unit(0.2, "cm"),
			legend.position="none") +
		xlab("") + ylab("") + ggtitle("") 

### PLOT: EXPR-PROTEIN-GROUP ---
p6 <- ggplot(dat, aes(x=SampleID, y=CancerType)) + 
		geom_tile(aes(fill=Group.Expr.protein), color="black", size=0.3) +
		scale_fill_manual(values=c("#d3cdab","#6f6d76","#c18377","#85d988")) +
		theme(
			axis.text = element_blank(),
			#axis.text.x = element_text(size = 3, color="black", angle=90, hjust=1, vjust=0.5),
			axis.title = element_blank(),
			plot.title = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.ticks.x =element_blank(),
			axis.ticks.y =element_blank(),
			panel.background = element_rect(fill="white"),
			legend.text = element_text(size = 3, color="black"),
			legend.title = element_text(size = 3, color="black"),
			legend.key.size = unit(0.2, "cm"),
			legend.position="none") +
		xlab("") + ylab("") + ggtitle("") 

### PLOT: BAP1-GROUP ---
p7 <- ggplot(dat, aes(x=SampleID, y=CancerType)) + 
		geom_tile(aes(fill=Group.3p21), color="black", size=0.3) +
		scale_fill_manual(values=c("#fffc00","blue")) +
		theme(
			axis.text = element_blank(),
			#axis.text.x = element_text(size = 3, color="black", angle=90, hjust=1, vjust=0.5),
			axis.title = element_blank(),
			plot.title = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.ticks.x =element_blank(),
			axis.ticks.y =element_blank(),
			panel.background = element_rect(fill="white"),
			legend.text = element_text(size = 3, color="black"),
			legend.title = element_text(size = 3, color="black"),
			legend.key.size = unit(0.2, "cm"),
			legend.position="none") +
		xlab("") + ylab("") + ggtitle("") 


### PLOT ---
file.plot <- file.path(dir.plot, "plot_annotation_header.pdf")
pdf(file.plot, height=2, width=3)
	grid.arrange(p1, p2, p3, p4, p5, p6, p7, ncol=1, nrow=7)
dev.off()


