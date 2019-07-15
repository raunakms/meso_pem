### Load Libraries -----
library("stringr")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/task/07_compare_mRNA_protein")
dir.dat <- file.path(dir.wrk, "output")
dir.plot <- file.path(dir.wrk, "plot")
dir.script <- file.path(dir.wrk, "scripts")

### DEFINE FILE ---
file.gobp <- file.path(dir.dat, "enrichment_protein_attenuation_go_bp.tsv")
file.gomf <- file.path(dir.dat, "enrichment_protein_attenuation_go_mf.tsv")
file.gocc <- file.path(dir.dat, "enrichment_protein_attenuation_go_cc.tsv")

### LOAD FILES ---
dat.gobp <- read.delim(file.gobp, header=T, stringsAsFactors=F)[,1:3]
dat.gomf <- read.delim(file.gomf, header=T, stringsAsFactors=F)[,1:3]
dat.gocc <- read.delim(file.gocc, header=T, stringsAsFactors=F)[,1:3]

### ADD GROUP INFO ---
dat.gobp$Group <- "GO-BP"
dat.gomf$Group <- "GO-MF"
dat.gocc$Group <- "GO-CC"

### MERGE DATA ---
dat <- rbind(dat.gobp, dat.gomf, dat.gocc)

### ADD VALUE ---
dat$nlogp <- -log10(dat$fdr)

### TRIM NAME ---
dat$Category <- str_replace(dat$Category, "GO_", "")
dat$Category <- str_replace_all(dat$Category, "_", " ")
dat$Category <- tolower(dat$Category)

### ORDER DATA ---
dat <- dat[order(dat$nlogp, decreasing=FALSE),]

### FACTORIZE DATA ---
dat$Category <- factor(dat$Category, level=dat$Category)
dat$Group <- factor(dat$Group, level=c("GO-BP","GO-MF","GO-CC"))

### PLOT ---
cpallete <- c("#d07a46","#5d4334","#afa29d")
p1 <- ggplot(dat, aes(x=Category, y=nlogp, fill=Group)) + 
		geom_bar(stat="identity") + 
		coord_flip() +
		scale_fill_manual(values=cpallete) +
		theme(
			axis.text = element_text(size = 7, color="black"),
			axis.title = element_text(size = 7, color="black"),
			plot.title = element_text(size = 7, color="black"),
			axis.ticks = element_line(size=0.4, color="black"),	
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			#panel.border = element_rect(color="black"),
			panel.background = element_rect(fill="white", color="black"),
			legend.text = element_text(size = 3.5, color="black"),
			legend.title = element_blank(),
			legend.key.size = unit(0.3, "cm"),
			legend.position="right") +
		ylab("") + xlab("") + 
		ggtitle("") 

### PLOT ---
file.plot <- file.path(dir.plot, "meso_enrichment_go.pdf")
pdf(file.plot, height=3, width=3)
	grid.arrange(p1, ncol=1, nrow=1)
dev.off()


##### REACTOME -----
file.reactome <- file.path(dir.dat, "enrichment_protein_attenuation_reactome.tsv")
dat.reactome <- read.delim(file.reactome, header=T, stringsAsFactors=F)[,1:3]
dat.reactome$Group <- "REACTOME"
dat.reactome$nlogp <- -log10(dat.reactome$fdr)

### TRIM NAME ---
dat.reactome$Category <- str_replace(dat.reactome$Category, "REACTOME_", "")
dat.reactome$Category <- str_replace_all(dat.reactome$Category, "_", " ")
dat.reactome$Category <- tolower(dat.reactome$Category)

### ORDER DATA ---
dat.reactome <- dat.reactome[order(dat.reactome$nlogp, decreasing=FALSE),]

### FACTORIZE DATA ---
dat.reactome$Category <- factor(dat.reactome$Category, level=dat.reactome$Category)
dat.reactome$Group <- as.factor(dat.reactome$Group)

### PLOT ---
p2 <- ggplot(dat.reactome, aes(x=Category, y=nlogp, fill=Group)) + 
		geom_bar(stat="identity", fill="#5d4334") + 
		coord_flip() +
		#scale_fill_manual(values=cpallete) +
		theme(
			axis.text = element_text(size = 7, color="black"),
			axis.title = element_text(size = 7, color="black"),
			plot.title = element_text(size = 7, color="black"),
			axis.ticks = element_line(size=0.4, color="black"),	
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			#panel.border = element_rect(color="black"),
			panel.background = element_rect(fill="white", color="black"),
			legend.text = element_text(size = 3.5, color="black"),
			legend.title = element_blank(),
			legend.key.size = unit(0.3, "cm"),
			legend.position="right") +
		ylab("") + xlab("") + 
		ggtitle("") 



### PLOT ---
file.plot <- file.path(dir.plot, "meso_enrichment_reactome.pdf")
pdf(file.plot, height=3, width=2.5)
	grid.arrange(p2, ncol=1, nrow=1)
dev.off()

