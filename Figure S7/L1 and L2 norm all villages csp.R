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


haps_village = read.csv('allhaps_demog.csv')
villages = read.csv('csp_villages.csv')

for (iteration in 1:nrow(villages)){
  village_seqs = data.frame(haplotype=character(), village=character(), village=character(),freq=double(), ID=character(), village=character(), month=character())
  print(paste0(Sys.time(), " --> starting iteration ", iteration, " of ", nrow(villages) ))
 
  
  
  for (rowNum in 1:nrow(haps_village)){
    if (haps_village[rowNum,"village"]==villages$villages[iteration]&&!is.na(haps_village[rowNum,"village"])){
        village_seqs <- rbind(village_seqs, haps_village[rowNum,])
    }
  }
  
  unique_ids <- village_seqs[,'ID']
  unique_ids <- unique (unique_ids)
  
  L1_D_sum =0
  L2_D_sum =0
  
  for (bootstrap in 1:100){

    random=sample(1:length(unique_ids), 2, replace=FALSE)
    ID1=unique_ids[random[1]]
    ID2=unique_ids[random[2]]
    
    ID1freqs = data.frame (A = numeric(288), C = numeric (288), G = numeric(288), Tn = numeric(288))
    ID2freqs = data.frame (A = numeric(288), C = numeric (288), G = numeric(288), Tn = numeric(288))
    
    
    for (nucleotide in 1:288){ # compute table of frequencies for ID1 and ID2
      for (rowNum in 1:nrow(village_seqs)){
        if (village_seqs$ID[rowNum]==ID1){
          seq=DNAString(village_seqs$haplotype[rowNum])
          bp = as.character(seq [nucleotide])
          if (bp == 'T'){
            ID1freqs$Tn[nucleotide] = ID1freqs$Tn[nucleotide] + village_seqs$freq[rowNum]
          }
          if (bp == 'A'){
            ID1freqs$A[nucleotide] = ID1freqs$A[nucleotide] + village_seqs$freq[rowNum]
          }
          if (bp == 'G'){
            ID1freqs$G[nucleotide] = ID1freqs$G[nucleotide] + village_seqs$freq[rowNum]
          }
          if (bp == 'C'){
            ID1freqs$C[nucleotide] = ID1freqs$C[nucleotide] + village_seqs$freq[rowNum]
          }
        }
        if (village_seqs$ID[rowNum]==ID2){
          seq=DNAString(village_seqs$haplotype[rowNum])
          bp = as.character(seq [nucleotide])
          if (bp == 'T'){
            ID2freqs$Tn[nucleotide] = ID2freqs$Tn[nucleotide] + village_seqs$freq[rowNum]
           }
           if (bp == 'A'){
             ID2freqs$A[nucleotide] = ID2freqs$A[nucleotide] + village_seqs$freq[rowNum]
           }
          if (bp == 'G'){
            ID2freqs$G[nucleotide] = ID2freqs$G[nucleotide] + village_seqs$freq[rowNum]
          }
          if (bp == 'C'){
            ID2freqs$C[nucleotide] = ID2freqs$C[nucleotide] + village_seqs$freq[rowNum]
          }
        }
      }
    }
    
    Dk = data.frame (A = numeric(288), C = numeric (288), G = numeric(288), Tn = numeric(288), L1Dk = numeric (288), L2Dk = numeric (288))
    for (nucleotide in 1:288){ #Compute L1 and L2 Dk for each nucleotide position
      Dk$A[nucleotide] = abs(ID1freqs$A[nucleotide]-ID2freqs$A[nucleotide])
      Dk$C[nucleotide] = abs(ID1freqs$C[nucleotide]-ID2freqs$C[nucleotide])
      Dk$G[nucleotide] = abs(ID1freqs$G[nucleotide]-ID2freqs$G[nucleotide])
      Dk$T[nucleotide] = abs(ID1freqs$T[nucleotide]-ID2freqs$T[nucleotide])
      Dk$L1Dk[nucleotide] = Dk$A[nucleotide] + Dk$C[nucleotide] + Dk$G[nucleotide] + Dk$T[nucleotide]
      Dk$L2Dk[nucleotide] = sqrt((Dk$A[nucleotide]*Dk$A[nucleotide])+(Dk$C[nucleotide]*Dk$C[nucleotide])+(Dk$G[nucleotide]*Dk$G[nucleotide])+(Dk$T[nucleotide]*Dk$T[nucleotide]))
    }
    
    L1_D = L2_D = 0
    
    for (nucleotide in 1:288){ #compute L1 and L2 statistics
      L1_D = L1_D + Dk$L1Dk[nucleotide]
      L2_D = L2_D + Dk$L2Dk[nucleotide]
    }
    
    L1_D_sum = L1_D_sum + L1_D
    L2_D_sum = L2_D_sum + L2_D
  }
  
  villages$L1_D[iteration]=L1_D_sum/bootstrap
  villages$L2_D[iteration]=L2_D_sum/bootstrap
  
}
  
write.csv(villages, 'csp L1L2norm village.csv')
##-----------------------------------




