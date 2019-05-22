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

haps_location = read.csv('ama allhaps_demog.csv')
locations = read.csv('ama_locations.csv')

for (iteration in 1:(nrow(locations))){
  location_seqs = data.frame(haplotype=character(), location=character(), village=character(),freq=double(), ID=character(), household=character(), month=character())
  print(paste0(Sys.time(), " --> starting iteration ", iteration, " of ", nrow(locations) ))
 
  
  
  for (rowNum in 1:nrow(haps_location)){
    if (haps_location[rowNum,"location"]==locations$locations[iteration]&&!is.na(haps_location[rowNum,"location"])){
        location_seqs <- rbind(location_seqs, haps_location[rowNum,])
    }
  }
  
  unique_ids <- location_seqs[,'ID']
  unique_ids <- unique (unique_ids)
  
  L1_D_sum =0
  L2_D_sum =0
  
  for (bootstrap in 1:100){

    random=sample(1:length(unique_ids), 2, replace=FALSE)
    ID1=unique_ids[random[1]]
    ID2=unique_ids[random[2]]
    
    ID1freqs = data.frame (A = numeric(300), C = numeric (300), G = numeric(300), Tn = numeric(300))
    ID2freqs = data.frame (A = numeric(300), C = numeric (300), G = numeric(300), Tn = numeric(300))
    
    
    for (nucleotide in 1:300){ # compute table of frequencies for ID1 and ID2
      for (rowNum in 1:nrow(location_seqs)){
        if (location_seqs$ID[rowNum]==ID1){
          seq=DNAString(location_seqs$haplotype[rowNum])
          bp = as.character(seq [nucleotide])
          if (bp == 'T'){
            ID1freqs$Tn[nucleotide] = ID1freqs$Tn[nucleotide] + location_seqs$freq[rowNum]
          }
          if (bp == 'A'){
            ID1freqs$A[nucleotide] = ID1freqs$A[nucleotide] + location_seqs$freq[rowNum]
          }
          if (bp == 'G'){
            ID1freqs$G[nucleotide] = ID1freqs$G[nucleotide] + location_seqs$freq[rowNum]
          }
          if (bp == 'C'){
            ID1freqs$C[nucleotide] = ID1freqs$C[nucleotide] + location_seqs$freq[rowNum]
          }
        }
        if (location_seqs$ID[rowNum]==ID2){
          seq=DNAString(location_seqs$haplotype[rowNum])
          bp = as.character(seq [nucleotide])
          if (bp == 'T'){
            ID2freqs$Tn[nucleotide] = ID2freqs$Tn[nucleotide] + location_seqs$freq[rowNum]
           }
           if (bp == 'A'){
             ID2freqs$A[nucleotide] = ID2freqs$A[nucleotide] + location_seqs$freq[rowNum]
           }
          if (bp == 'G'){
            ID2freqs$G[nucleotide] = ID2freqs$G[nucleotide] + location_seqs$freq[rowNum]
          }
          if (bp == 'C'){
            ID2freqs$C[nucleotide] = ID2freqs$C[nucleotide] + location_seqs$freq[rowNum]
          }
        }
      }
    }
    
    Dk = data.frame (A = numeric(300), C = numeric (300), G = numeric(300), Tn = numeric(300), L1Dk = numeric (300), L2Dk = numeric (300))
    for (nucleotide in 1:300){ #Compute L1 and L2 Dk for each nucleotide position
      Dk$A[nucleotide] = abs(ID1freqs$A[nucleotide]-ID2freqs$A[nucleotide])
      Dk$C[nucleotide] = abs(ID1freqs$C[nucleotide]-ID2freqs$C[nucleotide])
      Dk$G[nucleotide] = abs(ID1freqs$G[nucleotide]-ID2freqs$G[nucleotide])
      Dk$T[nucleotide] = abs(ID1freqs$T[nucleotide]-ID2freqs$T[nucleotide])
      Dk$L1Dk[nucleotide] = Dk$A[nucleotide] + Dk$C[nucleotide] + Dk$G[nucleotide] + Dk$T[nucleotide]
      Dk$L2Dk[nucleotide] = sqrt((Dk$A[nucleotide]*Dk$A[nucleotide])+(Dk$C[nucleotide]*Dk$C[nucleotide])+(Dk$G[nucleotide]*Dk$G[nucleotide])+(Dk$T[nucleotide]*Dk$T[nucleotide]))
    }
    
    L1_D = L2_D = 0
    
    for (nucleotide in 1:300){ #compute L1 and L2 statistics
      L1_D = L1_D + Dk$L1Dk[nucleotide]
      L2_D = L2_D + Dk$L2Dk[nucleotide]
    }
    
    L1_D_sum = L1_D_sum + L1_D
    L2_D_sum = L2_D_sum + L2_D
  }
  
  locations$L1_D[iteration]=L1_D_sum/bootstrap
  locations$L2_D[iteration]=L2_D_sum/bootstrap
  
}
  
write.csv(locations, 'ama location L1L2norm.csv')
##-----------------------------------




