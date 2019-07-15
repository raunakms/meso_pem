### DEFINE LIBRARIES --
library("stringr")
#library("BiocGenerics")
library("copynumber")

### DEFINe PATH ---
dir.wrk <- file.path("//jbrcsrv009/CollinsGroup/Raunak/projects/MESO_peritoneal/data/cnv")
#dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/cnv")
dir.data <- file.path(dir.wrk, "main_calls")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")

### DEFINE FUNCTION FILE ---
#file.func <- file.path("//jbrcsrv009/CollinsGroup/Raunak/softwares/bdvtools/copynumber/main.r")
#file.func <- file.path("/Data/Raunak/softwares/bdvtools/copynumber/main.r")

### DEFINE FILES ---
file.des <- file.path(dir.data, "design_table_3p21genes.tsv")
file.dat <- file.path(dir.data, "MESO-PeM_cnv_seg_calls_for_chromFreq_plot.tsv")

### LOAD DESIGN TABLE ---
des <- read.delim(file.des, header=TRUE, stringsAsFactors=FALSE)
grp0 <- des$Sample.ID[which(des$Group.3p21 == 0)]
grp1 <- des$Sample.ID[which(des$Group.3p21 == 1)]

#> grp0
# [1] "MESO-01"  "MESO-03"  "MESO-04"  "MESO-08"  "MESO-11"  "MESO-12"  "MESO-13"  "MESO-15"  "MESO-18A" "MESO-18E"
#[11] "MESO-19" 
#> grp1
#[1] "MESO-02" "MESO-05" "MESO-06" "MESO-07" "MESO-09" "MESO-10" "MESO-14" "MESO-17"

### Load CNA DATA ----
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
#dat$sampleID <- factor(dat$sampleID, levels=c(grp0, grp1))

### TO EDIT R-script ----
trace("interpolate.pcf", edit=TRUE)

### PLOT ---
#source(file.func)
file.plot <- file.path(dir.plot, "cnv_freq_plot_test.pdf")
pdf(file.plot, width=6.8, height=3)
	plotFreq(segments=dat, thres.gain=0.2, thres.loss = -0.2, pos.unit = "bp")
dev.off()


### PLOT ---
file.plot <- file.path(dir.plot, "cnv_freq_plot_revised.pdf")
pdf(file.plot, width=6.8, height=3)
  plotFreq(segments=dat, thres.gain=0.357, thres.loss = -0.357, pos.unit = "bp")
dev.off()

file.plot <- file.path(dir.plot, "cnv_freq_plot_chr3.pdf")
pdf(file.plot, width=3, height=2)
  plotFreq(segments=dat, thres.gain=0.2, thres.loss = -0.2, pos.unit = "bp", chrom=3)
dev.off()


### PLOTFreq COSTOM ---
file.plot <- file.path(dir.plot, "cnv_freq_plot_costom.pdf")
pdf(file.plot, width=6.8, height=3)
  plotFreq(segments=dat, thres.gain=0.357, thres.loss = -0.357, pos.unit = "bp", 
          main="", ylim=c(-100,100), continuous=TRUE)
dev.off()



### CIRCLE ---
plotCircle(segments=dat, thres.gain=0.2, thres.loss = -0.2, pos.unit = "bp", assembly = "hg19", freq.colors = c("red","blue"))

### HEATMAP ---
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)
dat$sampleID[which(dat$sampleID == "MESO-11")] <- "S-11"
dat$sampleID[which(dat$sampleID == "MESO-01")] <- "R-01"
dat$sampleID[which(dat$sampleID == "MESO-08")] <- "Q-08"
dat$sampleID[which(dat$sampleID == "MESO-15")] <- "P-15"
dat$sampleID[which(dat$sampleID == "MESO-04")] <- "O-04"
dat$sampleID[which(dat$sampleID == "MESO-12")] <- "N-12"
dat$sampleID[which(dat$sampleID == "MESO-19")] <- "M-19"
dat$sampleID[which(dat$sampleID == "MESO-13")] <- "L-13"
dat$sampleID[which(dat$sampleID == "MESO-03")] <- "K-03"
dat$sampleID[which(dat$sampleID == "MESO-18E")] <- "J-18E"
dat$sampleID[which(dat$sampleID == "MESO-18A")] <- "I-18A"


dat$sampleID[which(dat$sampleID == "MESO-17")] <- "H-17"
dat$sampleID[which(dat$sampleID == "MESO-02")] <- "G-02"
dat$sampleID[which(dat$sampleID == "MESO-07")] <- "F-07"
dat$sampleID[which(dat$sampleID == "MESO-10")] <- "E-10"
dat$sampleID[which(dat$sampleID == "MESO-05")] <- "D-05"
dat$sampleID[which(dat$sampleID == "MESO-06")] <- "C-06"
dat$sampleID[which(dat$sampleID == "MESO-09")] <- "B-09"
dat$sampleID[which(dat$sampleID == "MESO-14")] <- "A-14"

dat <- dat[order(dat$sampleID, decreasing=TRUE),]



file.plot <- file.path(dir.plot, "cnv_heatmap_per_sample_new.pdf")
pdf(file.plot, width=6.8, height=3)

  plotHeatmap(segments=dat, upper.lim=0.2, layout=c(1,1), 
              sep.samples=0, sample.line=0, main="", 
              colors=c("blue","white","red"), chrom.col="#bbbbbb")

dev.off()





### plotAberration ---
file.plot <- file.path(dir.plot, "cnv_aberration_per_sample.pdf")
pdf(file.plot, width=6.8, height=3)

  plotAberration(segments=dat, thres.gain=0.2, thres.loss = -0.2, 
                sep.samples=0, sample.line=0,
                pos.unit = "bp", main="", colors=c("blue","red"))

dev.off()




  plotHeatmap(segments=dat, upper.lim=0.2, layout=c(1,1), 
              sep.samples=0, sample.line=0, main="", chrom=c("1","2"),
              colors=c("blue","white","red"), chrom.col="black")



  plotHeatmap(segments=dat, thres.gain=0.2, thres.loss = -0.2, 
                sep.samples=0, sample.line=0, chrom="1",
                pos.unit = "bp", main="", colors=c("blue","red"))

chrom=c(1,3,21),



### CHECK ---
file.df <- file.path(dir.data, "MESO-PeM_bgr.bgr")
df <- read.delim(file.df, header=F, stringsAsFactors=F, skip=5)



##### FOR CLUSTER ----
#file.dat <- file.path(dir.data, "MESO-PeM_cnv_seg_calls_forcluster.tsv")
#dat <- read.delim(file.dat, header=T, stringsAsFactors=F)





####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

### Requires:
### pullOutContent

#Function with interpolates pcf-segments
interpolate.pcf = function(segments, x) {
  
  #Make sure segments is a data frame
  segments <- pullOutContent(res=segments,what="segments")
  
  usamp = unique(segments$sampleID)
  nsamp = length(usamp) 
  chrom = unique(x[,1])
  z = data.frame(cbind(x[,c(1:2)],matrix(0,nrow(x),nsamp)))
  #z = data.frame(x[,c(1,2)], matrix(0, nrow(x), nsamp))
  names(z) = c("chr","pos",usamp)

  for (i in 1:nsamp) {
    for (j in 1:length(chrom)) {
      
      index <- which((segments$sampleID==usamp[i]) & (segments$chrom==chrom[j]))
      if(length(index) == 0) next
      fitij = segments[index,]  

      #fitij = segments[segments$sampleID==usamp[i] & segments$chrom==chrom[j],] # BUG FOUND -- @RAUNAK
      
      #fitij = segments[index,]
      v = (c(fitij$start.pos[-1],10^9)+fitij$end.pos)/2
      xj = x[x[,1]==chrom[j],2]
      kj = rep(0,length(xj))
      for (k in rev(1:length(v))) {
        kj[xj <= v[k]] = k
      }
      z[z$chr==chrom[j],2+i] = fitij$mean[kj]

      #cat("PROCESSED: i=", i, "j =", j, "\n", sep=" ")
    }
  }

  return(z)
}

