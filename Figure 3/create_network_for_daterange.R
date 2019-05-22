##Create haplotype network for subset of villages for given date range (Figure 3D)

#### --------- load packages ----------------- ####
library(ggplot2)
library(readr)
library(dplyr)
library(gridExtra) 
library(ggthemes)
library(wesanderson)
library(tidyr)
library(lubridate)
library(Biostrings)
library(ape)
library(haplotypes)
library(pegas)



haps_village = read.csv("csp_allhaps_demog_person.csv")

fastaLines = c()
for (rowNum in 1:nrow(haps_village)){
  if (!is.na(haps_village[rowNum,"location"])){
    if (haps_village[rowNum,"month"]=="6-2013"||haps_village[rowNum,"month"]=="5-2013"||haps_village[rowNum,"month"]=="4-2013"){
      if (haps_village[rowNum, "location"]=="Maraka"||haps_village[rowNum, "location"]=="Miendo"||haps_village[rowNum, "location"]=="Misikhu"||haps_village[rowNum, "location"]=="Sitikho"||haps_village[rowNum, "location"]=="Webuye"){
    
        fastaLines = c(fastaLines, as.character(paste(">",haps_village[rowNum,"location"], sep = "")))
        fastaLines = c(fastaLines,as.character(haps_village[rowNum,"haplotype"]))
      }
    }
  }
 }
fileConn<-file('test.fasta')
writeLines(fastaLines, fileConn)
close(fileConn)

##  load data. 
haps = read.dna('test.fasta', format='fasta')
haps
count =0

compute_haps = haplotype(haps)
compute_haps

haps_net = haploNet (compute_haps)
plot(haps_net, size = attr(haps_net, "freq"), fast = FALSE)

ind_hap<-with (stack(setNames(attr(compute_haps, "index"), rownames(compute_haps))),table(hap=ind, individuals=rownames(haps)[values]))
ind_hap


plot(haps_net, size =log2(attr(haps_net, "freq")), col = "red", bg = "white",
     col.link = "black", lwd = 1, lty = 1, pie = ind_hap,
     labels = FALSE, font = 2, cex = 1, scale.ratio = 1,
     asp = 2, legend = FALSE, fast = FALSE, show.mutation = 0,
     threshold = c(1, 2))



