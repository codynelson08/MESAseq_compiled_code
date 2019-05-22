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

households = read.csv('ama_households.csv')
haps_village = read.csv('ama_allhaps_demog_person_age_decimaldate.csv')


for (iteration in 1:nrow(households)){

  fastaLines = c()
  for (rowNum in 1:nrow(haps_village)){
    if (haps_village[rowNum,"household"]==households$households[iteration]){
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
  
  freq_sum=0
  
  for (bootstrap in 1:10000){
    random_nums=sample(1:ncol(ind_hap), 2, replace=FALSE)
  
    ind_hap[,random_nums[1]]
    ind_hap[,random_nums[2]]
    
    shared=total=0
    for (i in 1:nrow(ind_hap)){
      if ((ind_hap[i,random_nums[1]]==1)&&(ind_hap[i, random_nums[2]]==0)){
        total = total+1
      }
      if ((ind_hap[i,random_nums[2]]==1)&&(ind_hap[i, random_nums[1]]==0)){
        total = total+1
      }
      if ((ind_hap[i,random_nums[1]]==1)&&(ind_hap[i, random_nums[2]]==1)){
        shared=shared+1
        total = total +1
      }
    }
    
   
    freq_sum = freq_sum+ (shared/total)
    
  }
  
 
  weighted_share_score = freq_sum/10000*100
  
  households$weighted_share_score[iteration]=weighted_share_score
}

write.csv(households, 'ama household weighted share score.csv')


