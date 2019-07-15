### DEFINE LIBRARIES --
library("copynumber")

### DEFINe PATH ---
dir.pm <- file.path("//jbrcsrv009/CollinsGroup/Raunak/HITnDRIVE/datasets/TCGA_MESO/data/cnv/main_calls")
dir.wrk <- file.path("//jbrcsrv009/CollinsGroup/Raunak/projects/MESO_peritoneal/data/cnv")
dir.pem <- file.path(dir.wrk, "main_calls")
dir.plot <- file.path(dir.wrk, "plot")

### DEFINE FILES ---
file.pm <- file.path(dir.pm, "TCGA_MESO_cnv_seg_calls_for_chromFreq_plot.tsv")
file.pem <- file.path(dir.pem, "MESO-PeM_cnv_seg_calls_for_chromFreq_plot.tsv")

### Load CNA DATA ----
dat.pm <- read.delim(file.pm, header=T, stringsAsFactors=F)
dat.pem <- read.delim(file.pem, header=T, stringsAsFactors=F)

### TO EDIT R-script ----
trace("interpolate.pcf", edit=TRUE)

### PLOTFreq COSTOM ---
file.plot <- file.path(dir.plot, "cnv_freq_plot_pem.pdf")
pdf(file.plot, width=6.8, height=3)
  plotFreq(segments=dat.pem, thres.gain=0.2, thres.loss = -0.2, pos.unit = "bp", 
          main="", ylim=c(-100,100), continuous=TRUE, col.gain="red", col.loss="blue")
dev.off()

file.plot <- file.path(dir.plot, "cnv_freq_plot_pm.pdf")
pdf(file.plot, width=6.8, height=3)
  plotFreq(segments=dat.pm, thres.gain=0.2, thres.loss = -0.2, pos.unit = "bp", 
          main="", ylim=c(-100,100), continuous=TRUE, col.gain="red", col.loss="blue")
dev.off()

