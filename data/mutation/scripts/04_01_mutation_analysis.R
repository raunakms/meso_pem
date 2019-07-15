### DEFINE LIBRARIES --
library("stringr")
library("ggplot2")
library("gridExtra")
library("reshape2")

### DEFINe PATH ---
dir.wrk <- file.path("/collinsgroup/Raunak/projects/MESO_peritoneal/data/mutation")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")

### DEFINE FILES ---
file.maf <- file.path(dir.data, "VPCLAGA.MESO_PeM.nonsilent_mutation.maf")

# LOAD DATA ---
dat.maf <- read.delim(file.maf, header=TRUE, stringsAsFactors=FALSE)

############################## --------------
# SUBSET DATA ---
dat <- subset(dat.maf, select=c("Tumor_Sample_Barcode","Hugo_Symbol"))
dat <- dat[!duplicated(dat),]

df.freq <- data.frame(Gene=names(table(dat$Hugo_Symbol)), Freq=as.numeric(table(dat$Hugo_Symbol)))
df.freq <- df.freq[order(df.freq$Freq, decreasing=T),]

df.freq.mut <- subset(df.freq, df.freq$Freq >=2)
df.freq.mut$Gene <- factor(df.freq.mut$Gene, levels=as.character(df.freq.mut$Gene))


### PLOT ---
p <- ggplot(df.freq.mut, aes(y=Freq, x=Gene)) +
		geom_bar(stat="Identity") +
		theme(
			axis.text.x = element_text(size = 5, color="black", angle=90, hjust=1, vjust=0.5),
			axis.text.y = element_text(size = 5, color="black"),
			axis.title = element_text(size = 5, color="black"),
			plot.title = element_text(size = 5, color="black"),
			axis.ticks = element_line(size=0.4),
			legend.position = "none") +
		ylab("No. of Samples") +
		xlab("") + 
		ggtitle("Mutated Genes")

file.plot <- file.path(dir.plot, "mutation_freq_stat.pdf")
pdf(file.plot, width=4, height=2)		
	grid.arrange(p, ncol=1, nrow=1)
dev.off()
############################## --------------

### PLOT Variant Frequencies by Sample ---
df.freq <- data.frame(SampleID=names(table(dat.maf$Tumor_Sample_Barcode)), Freq=as.numeric(table(dat.maf$Tumor_Sample_Barcode)))
df.freq$SampleID <- as.character(df.freq$SampleID)
df.freq <- df.freq[order(df.freq$Freq, decreasing=T),]

### PREPARE DATA ----
ids <- df.freq$SampleID
var <- unique(dat.maf$Variant_Classification)
mat <- matrix(0, nrow=length(ids), ncol=length(var), dimnames=list(ids, var))
for(i in 1:length(ids)){
	id <- ids[i]
	dat.temp <- subset(dat.maf, dat.maf$Tumor_Sample_Barcode == id)

	for(j in 1:ncol(mat)){
		mat[i, j] <- length(which(dat.temp$Variant_Classification == colnames(mat)[j]))
	}
}

### MELT DATA ---
df <- melt(mat)
colnames(df) <- c("SampleID","Variant","Freq")
df$SampleID <- as.character(df$SampleID)
df$Variant <- as.character(df$Variant)

### FACTORIZE DATA ---
items <- c("Missense_Mutation","Nonsense_Mutation","Splice_Site",
			"Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins","In_Frame_Del")
df$SampleID <- factor(df$SampleID, levels=df.freq$SampleID)
df$Variant <- factor(df$Variant, levels=items)


# Missense_Mutation	#698B69
# Nonsense_Mutation	#8B4789
# Splice_Site		#00CED1
# Frame_Shift_Ins	#FF8C00
# Frame_Shift_Del	#E7B98A
# In_Frame_Ins		#003366
# In_Frame_Del		#475ECC

cbPalette <- c("#698B69","#8B4789","#00CED1",
				"#FF8C00","#E7B98A","#003366","#475ECC")

### PLOT ---
p1 <- ggplot(df, aes(y=Freq, x=SampleID)) +
		geom_bar(aes(fill=Variant),stat="Identity") +
		scale_fill_manual(values=cbPalette) +
		theme(
			axis.text.x = element_text(size = 7, color="black", angle=90, hjust=1, vjust=0.5),
			axis.text.y = element_text(size = 7, color="black"),
			axis.title = element_text(size = 8, color="black"),
			plot.title = element_text(size = 16, color="black"),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),		
			axis.ticks = element_line(size=0.4, color="black"),
			panel.background = element_rect(fill="white", color="black"),
			legend.background=element_blank(),
			legend.text = element_text(size = 4, color="black"),
			legend.title = element_blank(),
			legend.key = element_blank(),
			legend.key.size = unit(0.2, "cm"),			
			legend.position = c(0.85,0.6)) +
		ylab("No. of Variants") +
		xlab("") + 
		ggtitle("")


### GET SNV RATE ---
file.snv <- file.path(dir.analysis, "snv_rate_MESO-PeM.tsv")
dat.snv <- read.delim(file.snv, header=TRUE, stringsAsFactors=FALSE)
dat.snv <- dat.snv[order(dat.snv$Rate.SNV, decreasing=T),]
dat.snv$SampleID <- factor(dat.snv$SampleID, levels=dat.snv$SampleID)

### PLOT ---
p2 <- ggplot(dat.snv, aes(y=Rate.SNV, x=SampleID)) +
		geom_bar(stat="Identity") +
		theme(
			axis.text.x = element_text(size = 7, color="black", angle=90, hjust=1, vjust=0.5),
			axis.text.y = element_text(size = 7, color="black"),
			axis.title = element_text(size = 8, color="black"),
			plot.title = element_text(size = 16, color="black"),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),				
			axis.ticks = element_line(size=0.4, color="black"),
			panel.background = element_rect(fill="white", color="black"),
			legend.background=element_blank(),
			legend.text = element_text(size = 7, color="black"),
			legend.title = element_blank(),
			legend.key = element_blank(),
			legend.key.size = unit(0.45, "cm"),			
			legend.position = "none") +
		ylab("Mutation Frequency (per Mb)") +
		xlab("") + 
		ggtitle("")

### PLOT ---
file.plot <- file.path(dir.plot, "snv_variant_summary_revised_20180626.pdf")
pdf(file.plot, height=4, width=3)
	grid.arrange(p1, p2, ncol=1, nrow=2)
dev.off()

### SUBSET CHR 19 ---
genes <- c("C19orf73","GYS1","OSCAR","PPP1R15A","PTOV1","RSPH6A","S1PR2","SBSN","ZNF28","ZNF419","ZNF600","ZNF813","ZNF814","ABCA7","APC2","CCDC130","CD209","CD320","CLEC4G","EVI5L","FBN3","KIAA1683","ZNF626","ZNF878")
d <- subset(dat.maf, dat.maf$Hugo_Symbol %in% genes)

write.table(d, file.path(dir.data, "VPCLAGA_MESO_PeM_chr19.maf"), sep="\t", row.names=F, col.names=T, quote=F)

d.bed <- subset(d, select=c("Chromosome","Start_Position","End_Position","Hugo_Symbol"))
d.bed <- d.bed[!duplicated(d.bed),]
d.bed$Chromosome <- paste("chr", d.bed$Chromosome, sep="")
d.bed <- d.bed[order(d.bed$Start_Position, decreasing=FALSE),]
write.table(d.bed, file.path(dir.data, "VPCLAGA_MESO_PeM_chr19.bed"), sep="\t", row.names=F, col.names=F, quote=F)
