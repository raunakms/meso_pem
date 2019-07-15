### LOAD LIBRARIES ---
library("stringr")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/cnv")
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

#y <- apply(dat, 1, function(x) paste(as.character(x[7]),as.numeric(x[8]),as.numeric(x[9]), sep=":"))
#write.table(y, file.path(dir.data, "MESO-Pem_gistic_cna_calls_chr_pos.tsv"), row.names=F, col.names=F, quote=F)

###
loc <- c("1p36","1q21","3p21","4p16",
		"5p15","7p11","7q11","8q21",
		"16p11","19p12","21p11","22q13","Xp11")

dat$Cytoband <- loc
dat$nlogq <- -log10(dat$Q.Bound)
dat <- dat[order(dat$nlogq, decreasing=F),]

### FACTORIZE DATA ---
dat$Cytoband <- factor(dat$Cytoband, levels=dat$Cytoband)
dat$Type <- factor(dat$Type, levels = c("CN Gain","CN Loss"))


### GENERATE PLOT ---
p1 <- ggplot(dat, aes(y=nlogq, x=Cytoband)) +
		geom_bar(aes(fill=Type), stat="Identity") +
		coord_flip() +
		scale_fill_manual(values=c("red","blue")) +
		theme(
			axis.text.x = element_text(size = 7, color="black"),
			axis.text.y = element_text(size = 7, color="black"),
			axis.title = element_text(size = 7, color="black"),
			plot.title = element_text(size = 10, color="black"),
			#axis.ticks = element_line(size=0.4),
			axis.ticks = element_line(size=0.4, color = "black"),
			panel.background = element_rect(fill = "white", colour = "black"),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			legend.text = element_text(size = 4, color="black"),
			legend.title = element_blank(),
			legend.key.size = unit(0.2, "cm"),			
			legend.position = "bottom") +
		xlab("") +
		ylab(expression( paste( "- log"[10] , " (qvalue)", sep="" ))) + 
		ggtitle("")

### PLOT ---
file.plot <- file.path(dir.plot, "gistic_calls.pdf")
pdf(file.plot, width=2, height=2.5)		
	grid.arrange(p1, ncol=1, nrow=1)
dev.off()

