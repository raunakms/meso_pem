### DEFINE LIBRARIES ---
library("stringr")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.wrk <- file.path("/collinsgroup/Raunak/projects/MESO_peritoneal/data/mutation")
dir.data <- file.path(dir.wrk, "main_calls")
dir.plot <- file.path(dir.wrk, "plot")
dir.des <- file.path("/collinsgroup/Raunak/projects/MESO_peritoneal/data/annotation")

### DEFINE PATH ---
file.maf <- file.path(dir.data, "VPCLAGA.MESO_PeM.nonsilent_mutation_withVAF.maf")
file.des <- file.path(dir.des, "design_table_3p21genes.tsv")

### LOAD DESIGN TABLE ---
des <- read.delim(file.des, header=TRUE, stringsAsFactors=FALSE)
grp0 <- des$Sample.ID[which(des$Group.3p21 == 0)]
grp1 <- des$Sample.ID[which(des$Group.3p21 == 1)]

### LOAD DATA ---
dat.maf <- read.delim(file.maf, header=TRUE, stringsAsFactors=FALSE)
dat.maf$Tumor_Sample_Barcode <- factor(dat.maf$Tumor_Sample_Barcode, levels=c(grp1,grp0))

### ADD SUBTYPE ---
dat.maf$Subtype[which(dat.maf$Tumor_Sample_Barcode %in% grp0)] <- "BAP1-INTACT"
dat.maf$Subtype[which(dat.maf$Tumor_Sample_Barcode %in% grp1)] <- "BAP1-DEL"
dat.maf$Subtype <- factor(dat.maf$Subtype, levels=c("BAP1-INTACT","BAP1-DEL"))

### PLOT ---
p <- ggplot(dat.maf, aes(VAF, fill=Subtype)) + 
		geom_density(color="#000000", alpha=0.5, stat="density", position="identity") +
		facet_wrap(~Tumor_Sample_Barcode, nrow=4, ncol=5, scales="fixed") +
		scale_fill_manual(values=c("#fffc00","#0000ff")) +
		#coord_cartesian(xlim=c(0,1)) + 
		theme(
			axis.text = element_text(size=8, color="black"),
			axis.title = element_text(size=10, color="black"),
			plot.title = element_text(size=10, color="black", hjust=0),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			strip.text = element_text(size = 8, color="black"),
			strip.background = element_rect(fill="#d9d9d9", color="#000000"),
			axis.ticks = element_line(size=0.4, color="black"),	
			panel.background = element_rect(fill="white", color="black"),
			legend.position="none") +
		xlab("Variant Allele Frequency") + 
		ylab("Frequency") + 
		ggtitle("") 	


### PLOT ---
file.plot <- file.path(dir.plot, "meso_VAF_density_plot.pdf")
pdf(file.plot, height=6, width=6.8)
	grid.arrange(p, ncol=1, nrow=1)
dev.off()




### PLOT ---
p1 <- ggplot(dat.maf, aes(VAF, fill=Subtype)) + 
		geom_vline(xintercept=0.5, linetype=2, color="#969696") + 
		geom_vline(xintercept=0.25, linetype=2, color="#969696") + 
		geom_vline(xintercept=0.75, linetype=2, color="#969696") + 
		geom_density(color="#000000", alpha=0.5, stat="density", position="identity") +
		facet_wrap(~Tumor_Sample_Barcode, nrow=19, ncol=1, scales="fixed", strip.position="right") +
		scale_fill_manual(values=c("#FFFC00","#0000FF")) +
		coord_cartesian(xlim=c(0,1), ylim=c(0,8)) + 
		theme(
			axis.text = element_text(size=5, color="#000000"),
			axis.title = element_text(size=10, color="#000000"),
			plot.title = element_text(size=10, color="#000000", hjust=0),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			strip.text = element_text(size = 0.1, color="#FFFFFF"),
			strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
			axis.ticks = element_line(size=0.4, color="#000000"),	
			panel.background = element_rect(fill="#FFFFFF", color=NA),
			legend.position="none") +
		xlab("Variant Allele Frequency") + 
		ylab("Frequency") + 
		ggtitle("") 	


### PLOT ---
file.plot <- file.path(dir.plot, "meso_VAF_density_plot_oneColumn.pdf")
pdf(file.plot, height=8, width=5)
	grid.arrange(p1, ncol=1, nrow=1)
dev.off()

