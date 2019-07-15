### LOAD LIBRARIES ---
library("stringr")
library("ggplot2")
library("gridExtra")


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
df$arm <- c("1p","1q","3p","4p","5p","7p","7q","8q","16p","19p","21p","22q","Xp")


### CHROMOSOME DATA ---
chr <- c(1:22, "X", "Y")
armp <- paste(chr, "p", sep="")
armq <- paste(chr, "q", sep="")

chrarm <- rep(NA, 48)
chrarm[seq(from=1,to=48,by=2)] <- armp
chrarm[seq(from=2,to=48,by=2)] <- armq

### PREPARE DATA ---
dm <- data.frame(Chr=chrarm)
dm$Chr <- as.character(dm$Chr)

dm$nlogp <- 0
dm$Status <- ""

### FILL DATA ---
for(i in 1:nrow(df)){
	index <- which(dm$Chr == df$arm[i])
	dm$nlogp[index] <- df$nlogp[i]
	dm$Status[index] <- df$Type[i]
}

### FACTORIZE DATA --
dm$Chr <- factor(dm$Chr, levels=chrarm)
dm$Status <- factor(dm$Status, levels=c("CN Gain", "CN Loss", ""))

### PLOT
cbpallete <- c("#ff0000","#0000ff","#ffffff")
p <- ggplot(dm, aes(y=nlogp, x=Chr, fill=Status)) +
		geom_bar(stat="Identity", width=0.8) +
		scale_fill_manual(values=cbpallete) +
		theme(
			axis.text.x = element_text(size=3, color="black"),
			axis.text.y = element_text(size=5, color="black"),
			axis.title = element_text(size=3, color="black"),
			plot.title = element_text(size=10, color="black", hjust=0),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),			
			axis.ticks = element_line(size=0.2, color="black"),
			strip.text = element_text(size=6, color="black"),
			strip.background = element_rect(fill="white", color="white"),
			panel.background = element_rect(fill="white", color="black"),
			legend.position="none") + 
		ylab("") +
		xlab("") + 
		ggtitle("")


### OUTPUT ---
file.plot <- file.path(dir.plot, "GISTIC_chrarm_plot.pdf")
pdf(file.plot, height=1.5, width=4)
	grid.arrange(p, ncol=1, nrow=1)
dev.off()

