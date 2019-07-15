### DEFINE LIBRARIES --
library("stringr")
#library("BiocGenerics")
library("copynumber")

### DEFINe PATH ---
dir.wrk <- file.path("//jbrcsrv009/CollinsGroup/Raunak/projects/MESO_peritoneal/data/cnv")
#dir.wrk <- file.path("/Data/Raunak/projects/MESO_peritoneal/data/cnv")
dir.data <- file.path(dir.wrk, "meso_wdpm")
dir.analysis <- file.path(dir.wrk, "analysis")
dir.plot <- file.path(dir.wrk, "plot")

### DEFINE FUNCTION FILE ---
#file.func <- file.path("//jbrcsrv009/CollinsGroup/Raunak/softwares/bdvtools/copynumber/main.r")
#file.func <- file.path("/Data/Raunak/softwares/bdvtools/copynumber/main.r")

### DEFINE FILES ---
file.dat <- file.path(dir.data, "WDPM-segment.tsv")

### Load CNA DATA ----
dat <- read.delim(file.dat, header=T, stringsAsFactors=F)

### TO EDIT R-script ----
trace("interpolate.pcf", edit=TRUE)

### PLOT ---
#source(file.func)
file.plot <- file.path(dir.plot, "cnv_freq_plot_wdpm.pdf")
pdf(file.plot, width=6.8, height=3)
	plotFreq(segments=dat[1:10,], thres.gain=0.2, thres.loss = -0.2, pos.unit = "bp")
dev.off()

file.plot <- file.path(dir.plot, "cnv_freq_plot_chr3.pdf")
pdf(file.plot, width=3, height=2)
  plotFreq(segments=dat, thres.gain=0.2, thres.loss = -0.2, pos.unit = "bp", chrom=3)
dev.off()


### PLOTFreq COSTOM ---
file.plot <- file.path(dir.plot, "cnv_freq_plot_costom.pdf")
pdf(file.plot, width=6.8, height=3)
  plotFreq(segments=dat, thres.gain=0.2, thres.loss = -0.2, pos.unit = "bp", 
          main="", percentLines=FALSE, chrom.col="white", ylim=c(-100,100),
          continuous=TRUE, col.gain="#652736", col.loss="#652736")
dev.off()



### CIRCLE ---
plotCircle(segments=dat, thres.gain=0.2, thres.loss = -0.2, pos.unit = "bp", assembly = "hg19", freq.colors = c("red","blue"))

### HEATMAP ---
file.plot <- file.path(dir.plot, "cnv_heatmap_per_sample.pdf")
pdf(file.plot, width=6.8, height=3)

  plotHeatmap(segments=dat, upper.lim=0.2, layout=c(1,1), 
              sep.samples=0, sample.line=0, main="",
              colors=c("blue","white","red"), chrom.col="black")

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

