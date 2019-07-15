### SetDirectories ----------------------------------------------------------
dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/mutation")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")
dir.cosmic <- file.path("/Data/Raunak/data_ref/cancer_gene_census")
dir.drug <- file.path("/Data/Raunak/data_ref/drugs_data/druggability")

### Define Files ------------------------------------------------------------
file.dat <- file.path(dir.data, "somatic_mutation_hg19_non_silent_non_dbsnp.tsv")
file.dat2 <- file.path(dir.analysis, "germline_mutation_samplewise_filterd_by_normals.tsv")
file.dat3 <- file.path(dir.analysis, "TCGA_MESO_mutation_calls.tsv")
file.cgc <- file.path(dir.cosmic, "cancer_gene_census_genelist.txt")
file.cosmic <- file.path(dir.cosmic, "cosmic_genelist.txt")
file.targetdb <- file.path(dir.drug, "target_db.tsv")
file.dgidb <- file.path(dir.drug, "dgidb_druggable_gene_categories.tsv")
file.oncoprint <- file.path("/Data/Raunak/softwares/bdvtools/oncoprint/get_oncoprint_plot.R")

### Load Libraries ----------------------------------------------------------
library("ggplot2")
library("gridExtra")
library("stringr")

### Load COSMIC & DRUG Data --------------------------------------------------
genes.cgc <- read.table(file.cgc, header=F, as.is="V1")$V1	
genes.cosmic <- read.table(file.cosmic, header=F, as.is="V1")$V1
genes.targetdb <- unique(read.delim(file.targetdb, header=T, as.is="Gene")$Gene)
genes.dgidb <- unique(read.delim(file.dgidb, header=T, as.is="entrez_gene_symbol")$entrez_gene_symbol)
			
### Load Mutation Data --------------------------------------------------------
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
colnames(dat)[126:141] <- str_replace_all(colnames(dat)[126:141], "[.]", "-")
dat <- dat[,c(1:12, 126:141)]

### Parse Files ---------------------------------------------------------------
dat.list <- list()
for(i in 13:28){
	id <- colnames(dat)[i]
	dat.list[[i]] <- cbind(SampleID=id, dat[which(dat[,i] != "."),1:12])
}
df <- do.call(rbind.data.frame, dat.list)
df$SampleID <- as.character(df$SampleID)
df <- df[,-10]


### Load data for Samples wof which Mutations were filtered agains collective normals ---------------
dat2 <- read.delim(file.dat2, header=T, stringsAsFactors=F)

### Match Germline Filtered Mutated Genes -----------------------------------------------------------
dat2$Match <- ""
for(i in 1:nrow(dat2)){
	gene.id <- dat2$Name[i]
	index <- which(df$Name == gene.id)
	dat2$Match[i] <- ifelse(length(index) == 0, "N", "Y")
}

### Load data for Pleural MESO TCGA ----------------------------------------------------------------
dat3 <- read.delim(file.dat3, header=T, stringsAsFactors=F)

dat2$TCGAMatch <- ""
for(i in 1:nrow(dat2)){
	gene.id <- dat2$Name[i]
	index <- which(dat3$Gene == gene.id)
	dat2$TCGAMatch[i] <- ifelse(length(index) == 0, "N", "Y")
}

index <- which((dat2$Match == "Y") | (dat2$TCGAMatch == "Y"))
dat2 <- dat2[index,]

### Merge data
df <- rbind(df, dat2[,-c(13,14)])

file.output <- file.path(dir.analysis, "somatic_mutation_non_silent_samplewise.tsv")
write.table(df, file.output, sep="\t", row.names=F, col.names=T, quote=F)







### Stat: Mutation --------------------------------------------------------------
dat.temp <- df[,c(1,9)]
dat.id.mutgene  <- dat.temp [!duplicated(dat.temp),]

dstat1 <- data.frame(SampleID=names(table(df$SampleID)), 
					FreqMut=as.numeric(table(df$SampleID)), 
					FreqGene=as.numeric(table(dat.id.mutgene$SampleID)))
dstat1$SampleID <- as.character(dstat1$SampleID)

### Stat: COSMIC/DRUGS -----------------------------------------------------------
id <- dstat1$SampleID
dstat1$FreqNOHIT <- dstat1$FreqCOSMIC <- dstat1$FreqCGC <- 0
dstat1$FreqNODRUGNHIT <- dstat1$FreqDGIDB <- dstat1$FreqTARGETDB <- 0

for(i in 1:length(id)){
	dat.temp <- subset(dat.id.mutgene, dat.id.mutgene$SampleID == id[i])
	dstat1$FreqCGC[i] <- length(intersect(dat.temp$Name, genes.cgc))
	dstat1$FreqCOSMIC[i] <- length(setdiff(intersect(dat.temp$Name, genes.cosmic), genes.cgc))
	dstat1$FreqNOHIT[i] <- dstat1$FreqGene[i] - (dstat1$FreqCGC[i] + dstat1$FreqCOSMIC[i])
	
	dstat1$FreqTARGETDB[i] <- length(intersect(dat.temp$Name, genes.targetdb))
	dstat1$FreqDGIDB[i] <- length(setdiff(intersect(dat.temp$Name, genes.dgidb), genes.targetdb))
	dstat1$FreqNODRUGNHIT[i] <- dstat1$FreqGene[i] - (dstat1$FreqTARGETDB[i] + dstat1$FreqDGIDB[i])	
}

file.output <- file.path(dir.analysis, "somatic_mutation_non_silent_samplewise_stats.tsv")
write.table(dstat1, file.output, sep="\t", row.names=F, col.names=T, quote=F)

### Plot Graph --------------------------------------------------------------
dstat1 <- dstat1[order(dstat1$FreqMut, decreasing=T),]
dstat1$SampleID <- factor(dstat1$SampleID, levels=dstat1$SampleID)

p1 <- ggplot(dstat1, aes(y=FreqMut, x=SampleID)) +
		geom_bar(stat="Identity") +
		theme(
			axis.text.x = element_text(size = 7, color="black", angle=90, hjust=1, vjust=0.5),
			axis.text.y = element_text(size = 7, color="black"),
			axis.title = element_text(size = 7, color="black"),
			plot.title = element_text(size = 10, color="black"),
			axis.ticks = element_line(size=0.4),
			legend.position = "none") +
		ylab("No. of Mutations") +
		xlab("") + 
		ggtitle("Distribution of Somatic Mutation Events")
		
p2 <- ggplot(dstat1, aes(y=FreqGene, x=SampleID)) +
		geom_bar(stat="Identity") +
		theme(
			axis.text.x = element_text(size = 7, color="black", angle=90, hjust=1, vjust=0.5),
			axis.text.y = element_text(size = 7, color="black"),
			axis.title = element_text(size = 7, color="black"),
			plot.title = element_text(size = 10, color="black"),
			axis.ticks = element_line(size=0.4),
			legend.position = "none") +
		ylab("No. of Mutated Genes") +
		xlab("") + 
		ggtitle("Distribution of Somatic Mutated Genes")

file.plot <- file.path(dir.plot, "mutation_samplewise_stat.pdf")
pdf(file.plot, width=7, height=4)		
	grid.arrange(p1, p2, ncol=1, nrow=2)
dev.off()


### Stat: Mutated Genes ------------------------------------------------------------
dstat2 <- data.frame(Gene=names(table(dat.id.mutgene$Name)), Freq=as.numeric(table(dat.id.mutgene$Name)))
dstat2 <- dstat2[order(dstat2$Freq, decreasing=T),]
dstat2$Gene <- factor(dstat2$Gene , levels=dstat2$Gene)

p3 <- ggplot(dstat2, aes(y=Freq, x=Gene)) +
		geom_bar(stat="Identity") +
		theme(
			axis.text.x = element_text(size = 7, color="black", angle=90, hjust=1, vjust=0.5),
			axis.text.y = element_text(size = 7, color="black"),
			axis.title = element_text(size = 7, color="black"),
			plot.title = element_text(size = 10, color="black"),
			axis.ticks = element_line(size=0.4),
			legend.position = "none") +
		ylab("No. of Samples") +
		xlab("") + 
		ggtitle("No. of Samples with Mutated Genes")

file.plot <- file.path(dir.plot, "mutation_freq_stat.pdf")
pdf(file.plot, width=25, height=3.5)		
	grid.arrange(p3, ncol=1, nrow=1)
dev.off()

### Mutation Oncoprint ------------------------------------------------------------		
dat.oncp <- dat.id.mutgene
dat.oncp$Status <- "MUT;"
colnames(dat.oncp) <- c("SampleID","Gene","Status")

source(file.oncoprint)
file.plot <- file.path(dir.plot, "mutation_oncoprint.pdf")  
pdf(file.plot, width=3, height=15)
	oncoprint(dat.oncp, orderGenes="memo.sort", file.genes=NA, file.samples=NA, circularize=FALSE, cex.text=4, cex.axis=5, cex.title=5)
dev.off()


### Mutation Mapper ------------------------------------------------------------		
z <- unlist(lapply(str_split(df$AAChange.ensGene, ":"), function(x) x[length(x)]))
aa.change <- unlist(lapply(str_split(z, "[.]"), function(x) x[2]))
aa.change[which(aa.change == "")] <- "."

dat.mutmap <- data.frame(Hugo_Symbol=df$Name,
						Protein_Change=aa.change,
						Sample_ID=df$SampleID,
						Mutation_Type=df$ExonicFunc.ensGene,
						Chromosome=df$Chr,
						Start_Position=df$Start,
						End_Position=df$End,
						Reference_Allele=df$Ref,
						Variant_Allele=df$Alt)

file.output <- file.path(dir.analysis, "meso_mutation_mapper.tsv")
write.table(dat.mutmap , file.output, sep="\t", row.names=F, col.names=T, quote=F)

### Prepare Data for HIT'nDRIVE -----------------------------------------------------------
df.alt <- df[,c(1,9)]
df.alt  <- df.alt[!duplicated(df.alt),]
colnames(df.alt) <- c("SampleID","Gene")
df.alt$Status <- "MUT;"

file.output <- file.path(dir.analysis, "alterations_samplewise_mutation.tsv")
write.table(df.alt, file.output, sep="\t", row.names=F, col.names=T, quote=F)
