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

households = read.csv('ama_case_households.csv')
haps_village = read.csv('ama_allhaps_demog_person.csv')
case_children = read.csv('ama_casechild_share.csv')


for (case_iteration in 1:nrow(households)){
  tablecount = 2
  print (case_iteration)

  
  for (household_iteration in 1:nrow(households)){
    
    fastaLines = c()
    
    for (rowNum in 1:nrow(haps_village)){ #pull out case child first
      if (haps_village[rowNum,"household"]==households$households[case_iteration] && haps_village[rowNum, 'person']=='case child'){
          fastaLines = c(fastaLines, as.character(paste(">case")))
          fastaLines = c(fastaLines,as.character(haps_village[rowNum,"haplotype"]))
      }
    }
    
    for (rowNum in 1:nrow(haps_village)){ #pull out household members of matching house
      if (haps_village[rowNum,"household"]==households$households[household_iteration] && haps_village[rowNum, 'person']=='case household member'){
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
    haps_net = haploNet (compute_haps)
    ind_hap<-with (stack(setNames(attr(compute_haps, "index"), rownames(compute_haps))),table(hap=ind, individuals=rownames(haps)[values]))
    ind_hap
   
    
    freq_sum=0
    
    for (bootstrap in 1:10000){
      random_household_member=sample(1:(ncol(ind_hap)-1), 1)
      
      ind_hap[,random_household_member]
      ind_hap[,ncol(ind_hap)]
      
      shared=total=0
      for (i in 1:nrow(ind_hap)){
        if ((ind_hap[i,random_household_member]==1)&&(ind_hap[i,ncol(ind_hap)]==0)){
          total = total+1
        }
        if ((ind_hap[i,random_household_member]==0)&&(ind_hap[i,ncol(ind_hap)]==1)){
          total = total+1
        }
        if ((ind_hap[i,random_household_member]==1)&&(ind_hap[i,ncol(ind_hap)]==1)){
          shared=shared+1
          total = total +1
        }
      }
      
      freq_sum = freq_sum+ (shared/total)
    }
    
    weighted_share_score = freq_sum/10000*100
    
    case_children[case_iteration,tablecount] = weighted_share_score
    tablecount=tablecount+1
    
  }
}


write.csv(case_children, 'ama weighted share score case_child matrix.csv')


