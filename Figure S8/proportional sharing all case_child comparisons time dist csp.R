#csp proportional sharing between all case children

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
library (rlist)
library(Biostrings)


haps_village = read.csv('csp_allhaps_demog_person_age_decimaldate_latlong.csv')
case_children = read.csv('csp_case_children.csv') 
households = read.csv('csp_case_households.csv')

case_children_edited = case_children
#for (i in 1:nrow(case_childs)){ #select out household members from our 50 households
#  for (j in 1:nrow(households)){
#    if (households$households[j]==case_childs$household[i]){
#      case_children_edited = rbind (case_children_edited, case_childs[i,])
#    }
#  }
#}
size = nrow(case_children_edited)^2
cc_comps = data.frame(case_child = character(size), comp_cc = character(size), cc_date = double(size), comp_cc_date=double(size),temporal_dist = double(size), geog_dist= double(size), share_score=double(size), stringsAsFactors=FALSE)

for (cc_iteration in 1:nrow(case_children_edited)){
  
  print(paste0(Sys.time(), " --> starting case_child ", cc_iteration, " of ", nrow(case_children_edited) ))
  cat("\n")
  
  for (comp_cc_iteration in 1:nrow(case_children_edited)){
    fastaLines = c()
    
    for (rowNum in 1:nrow(haps_village)){ #pull out first case child
      if (haps_village[rowNum,"ID"]==case_children_edited$case_children[cc_iteration]){
        fastaLines = c(fastaLines, as.character(paste(">h",haps_village$ID[rowNum], sep="")))
        fastaLines = c(fastaLines,as.character(haps_village[rowNum,"haplotype"]))
      }
    }
    
    for (rowNum in 1:nrow(haps_village)){ #pull comparison case child
      if (haps_village[rowNum,"ID"]==case_children_edited$case_children[comp_cc_iteration]){
        fastaLines = c(fastaLines, as.character(paste(">ccc",haps_village$ID[rowNum], sep = "")))
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
    ind_hap<-with (stack(setNames(attr(compute_haps, "index"), rownames(compute_haps))),table(hap=ind, individuals=rownames(haps)[values]))
    ind_hap
    
    shared=total=0
    for (i in 1:nrow(ind_hap)){
      if ((ind_hap[i,1]==1)&&(ind_hap[i,2]==0)){
        total = total+1
      }
      if ((ind_hap[i,1]==0)&&(ind_hap[i,2]==1)){
        total = total+1
      }
      if ((ind_hap[i,1]==1)&&(ind_hap[i,2]==1)){
        shared=shared+1
        total = total +1
      }
    }
    
    weighted_share_score = (shared/total)*100
    
    line = (cc_iteration-1)*nrow(case_children_edited) + comp_cc_iteration
    
    for (i in 1:nrow(haps_village)){ #pick out lat/long and date for household
      if (haps_village$ID[i] == case_children_edited$case_children[cc_iteration]){
        cc_lat = haps_village$latitude[i]
        cc_long = haps_village$longitude[i]
        cc_date = haps_village$date[i]
      }
    }
    for (i in 1:nrow(haps_village)){ #pick out lat/long and date for comp_household
      if (haps_village$ID[i] == case_children_edited$case_children[comp_cc_iteration]){
        comp_cc_lat = haps_village$latitude[i]
        comp_cc_long = haps_village$longitude[i]
        comp_cc_date = haps_village$date[i]
      }
    }
    
    cc_comps$case_child[line] = as.character(case_children_edited$case_children[cc_iteration])
    cc_comps$comp_cc[line] = as.character(case_children_edited$case_children[comp_cc_iteration])
    cc_comps$geog_dist[line] = sqrt(((comp_cc_lat-cc_lat)*(comp_cc_lat-cc_lat)) + ((comp_cc_long-cc_long)*(comp_cc_long-cc_long)))
    cc_comps$cc_date[line]=cc_date
    cc_comps$comp_cc_date[line]=comp_cc_date
    cc_comps$temporal_dist[line] = abs(comp_cc_date - cc_date)
    cc_comps$share_score[line] = weighted_share_score
    
  }
}


write.csv(cc_comps, 'csp weighted share score all283 case_child comp time dist.csv')


