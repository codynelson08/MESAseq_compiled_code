#### --------- load packages ----------------- ####
library(ggplot2)
library(readr)
library(dplyr)
library(gridExtra) 
library(ggthemes)
library(wesanderson)
library(tidyr)
library(lubridate)
library(phangorn)
library(pegas)


haps_village = read.csv ('allhaps_demog.csv')
locations = read.csv('ama_locations.csv')


for (iteration in 1:nrow(locations)){
  
  fastaLines = c()
  for (rowNum in 1:nrow(haps_village)){
    if (haps_village[rowNum,"location"]==locations$locations[iteration]&& !is.na (haps_village[rowNum,'location']){
        fastaLines = c(fastaLines, as.character(paste(">",haps_village$ID[rowNum], sep = "")))
        fastaLines = c(fastaLines,as.character(haps_village[rowNum,"haplotype"]))
      
    }
  }
  fileConn<-file('test.fasta')
  writeLines(fastaLines, fileConn)
  close(fileConn)
  
  haps = read.dna('test.fasta', format='fasta')
  haps
  count =0
  
  compute_haps = haplotype(haps)
  compute_haps
  
  haps_net = haploNet (compute_haps)
  
  ind_hap<-with (stack(setNames(attr(compute_haps, "index"), rownames(compute_haps))),table(hap=ind, individuals=rownames(haps)[values]))
  ind_hap  
  
  shared_count=0
  
  for (bootstrap in 1:10000){
    random_nums=sample(1:ncol(ind_hap), 2, replace=FALSE)
  
    ind_hap[,random_nums[1]]
    ind_hap[,random_nums[2]]
    
    shared = FALSE
    for (i in 1:nrow(ind_hap)){
      if ((ind_hap[i,random_nums[1]]==1)&&(ind_hap[i, random_nums[2]]==1)){
        shared=TRUE
      }
    }
    
    if (shared == TRUE){
      shared_count = shared_count+1
    }
  }
  
  shared_count
  shared_count_freq = shared_count/10000*100
  
  locations$share_score[iteration]=shared_count_freq
}

write.csv(locations, 'ama location share score.csv')


