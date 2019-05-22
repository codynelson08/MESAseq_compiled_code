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
library(Biostrings)


haps_village = read.csv('csp_allhaps_demog_person_age_decimaldate_latlong.csv')
case_children = read.csv('csp_case_children.csv') 
households = read.csv('csp_case_households.csv')

case_children_edited = case_children

#for (i in 1:nrow(case_children)){ #select out household members from our 50 households
#  for (j in 1:nrow(households)){
#    if (households$households[j]==case_children$household[i]){
#      case_children_edited = rbind (case_children_edited, case_children[i,])
#    }
#  }
#}
size = nrow(case_children_edited)^2
cc_comps = data.frame(household_member = character(size), comp_cc = character(size), cc_date = double (size),comp_cc_date = double(size), temporal_dist = double(size), geog_dist= double(size), L1norm=double(size),L2norm=double(size), stringsAsFactors=FALSE)

for (cc_iteration in 152:nrow(case_children_edited)){

  print(paste0(Sys.time(), " --> starting case_child ", cc_iteration, " of ", nrow(case_children_edited) ))
  cat("\n")
  
  for (comp_cc_iteration in 1:nrow(case_children_edited)){
    cc_seqs = data.frame(haplotype=character(), location=character(), village=character(),freq=double(), ID=character(), household=character(), month=character())
    comp_cc_seqs = data.frame(haplotype=character(), location=character(), village=character(),freq=double(), ID=character(), household=character(), month=character())
    
    for (rowNum in 1:nrow(haps_village)){
      if (haps_village[rowNum,"ID"]==case_children_edited$case_ID[cc_iteration]){
          cc_seqs <- rbind(cc_seqs, haps_village[rowNum,])
      }
      if (haps_village[rowNum,"ID"]==case_children_edited$case_ID[comp_cc_iteration]){
          comp_cc_seqs <- rbind(comp_cc_seqs, haps_village[rowNum,])
      }
    }
    
  
    
    L1_D_sum =0
    L2_D_sum =0
    
  
    
      ID1=cc_seqs$ID[1] #ID1 will be case
      ID2=comp_cc_seqs$ID[1] #ID2 will be comp_cc random member
      
      ID1freqs = data.frame (A = numeric(288), C = numeric (288), G = numeric(288), Tn = numeric(288))
      ID2freqs = data.frame (A = numeric(288), C = numeric (288), G = numeric(288), Tn = numeric(288))
      
      
      for (nucleotide in 1:288){ # compute table of frequencies for ID1 and ID2
        for (rowNum in 1:nrow(cc_seqs)){ #compute for ID1-case seqs first
            seq=DNAString(cc_seqs$haplotype[rowNum])
            bp = as.character(seq [nucleotide])
            if (bp == 'T'){
              ID1freqs$Tn[nucleotide] = ID1freqs$Tn[nucleotide] + cc_seqs$freq[rowNum]
            }
            if (bp == 'A'){
              ID1freqs$A[nucleotide] = ID1freqs$A[nucleotide] + cc_seqs$freq[rowNum]
            }
            if (bp == 'G'){
              ID1freqs$G[nucleotide] = ID1freqs$G[nucleotide] + cc_seqs$freq[rowNum]
            }
            if (bp == 'C'){
              ID1freqs$C[nucleotide] = ID1freqs$C[nucleotide] + cc_seqs$freq[rowNum]
            }
          }
        
        for (rowNum in 1:nrow(comp_cc_seqs)){
            seq=DNAString(comp_cc_seqs$haplotype[rowNum])
            bp = as.character(seq [nucleotide])
            if (bp == 'T'){
              ID2freqs$Tn[nucleotide] = ID2freqs$Tn[nucleotide] + comp_cc_seqs$freq[rowNum]
            }
            if (bp == 'A'){
              ID2freqs$A[nucleotide] = ID2freqs$A[nucleotide] + comp_cc_seqs$freq[rowNum]
            }
            if (bp == 'G'){
              ID2freqs$G[nucleotide] = ID2freqs$G[nucleotide] + comp_cc_seqs$freq[rowNum]
            }
            if (bp == 'C'){
              ID2freqs$C[nucleotide] = ID2freqs$C[nucleotide] + comp_cc_seqs$freq[rowNum]
            }
          
          }
      }
      
      Dk = data.frame (A = numeric(288), C = numeric (288), G = numeric(288), Tn = numeric(288), L1Dk = numeric (288), L2Dk = numeric (288))
      for (nucleotide in 1:288){ #Compute L1 and L2 Dk for each nucleotide position
        Dk$A[nucleotide] = abs(ID1freqs$A[nucleotide]-ID2freqs$A[nucleotide])
        Dk$C[nucleotide] = abs(ID1freqs$C[nucleotide]-ID2freqs$C[nucleotide])
        Dk$G[nucleotide] = abs(ID1freqs$G[nucleotide]-ID2freqs$G[nucleotide])
        Dk$Tn[nucleotide] = abs(ID1freqs$T[nucleotide]-ID2freqs$Tn[nucleotide])
        Dk$L1Dk[nucleotide] = Dk$A[nucleotide] + Dk$C[nucleotide] + Dk$G[nucleotide] + Dk$T[nucleotide]
        Dk$L2Dk[nucleotide] = sqrt((Dk$A[nucleotide]*Dk$A[nucleotide])+(Dk$C[nucleotide]*Dk$C[nucleotide])+(Dk$G[nucleotide]*Dk$G[nucleotide])+(Dk$T[nucleotide]*Dk$T[nucleotide]))
      }
      
      L1_D = L2_D = 0
      
      for (nucleotide in 1:288){ #compute L1 and L2 statistics
        L1_D = L1_D + Dk$L1Dk[nucleotide]
        L2_D = L2_D + Dk$L2Dk[nucleotide]
      }
   
      
      
      line = (cc_iteration-1)*nrow(case_children_edited) + comp_cc_iteration
      
      for (i in 1:nrow(haps_village)){ #pick out lat/long and date for household
        if (haps_village$ID[i] == case_children_edited$case_ID[cc_iteration]){
          cc_lat = haps_village$latitude[i]
          cc_long = haps_village$longitude[i]
          cc_date = haps_village$date[i]
        }
      }
      for (i in 1:nrow(haps_village)){ #pick out lat/long and date for comp_household
        if (haps_village$ID[i] == case_children_edited$case_ID[comp_cc_iteration]){
          comp_cc_lat = haps_village$latitude[i]
          comp_cc_long = haps_village$longitude[i]
          comp_cc_date = haps_village$date[i]
        }
      }
      
      cc_comps$household_member[line] = as.character(case_children_edited$case_ID[cc_iteration])
      cc_comps$comp_cc[line] = as.character(case_children_edited$case_ID[comp_cc_iteration])
      cc_comps$geog_dist[line] = sqrt(((comp_cc_lat-cc_lat)*(comp_cc_lat-cc_lat)) + ((comp_cc_long-cc_long)*(comp_cc_long-cc_long)))
      cc_comps$cc_date[line] = cc_date
      cc_comps$comp_cc_date [line] = comp_cc_date
      cc_comps$temporal_dist[line] = abs(comp_cc_date - cc_date)
      cc_comps$L1norm[line] = L1_D
      cc_comps$L2norm[line] = L2_D
    }
    
  }

  
write.csv(cc_comps, 'csp L1L2norm all283 household_member comp time dist')
##-----------------------------------




