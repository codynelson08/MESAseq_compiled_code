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


haps_village = read.csv('csp_allhaps_demog_person_age_decimaldate.csv')
households = read.csv('csp_case_households.csv')

for (iteration in 1:nrow(households)){
  household_seqs = data.frame(haplotype=character(), location=character(), village=character(),freq=double(), ID=character(), household=character(), month=character())
  case_seqs = data.frame(haplotype=character(), location=character(), village=character(),freq=double(), ID=character(), household=character(), month=character())
  print(paste0(Sys.time(), " --> starting iteration ", iteration, " of ", nrow(households) ))
 
  for (rowNum in 1:nrow(haps_village)){
    if (haps_village[rowNum,"household"]==households$households[iteration] && haps_village[rowNum,"person"]=='case child'){
        case_seqs <- rbind(case_seqs, haps_village[rowNum,])
    }
    if (haps_village[rowNum,"household"]==households$households[iteration] && haps_village[rowNum,"person"]=='case household member'){
        household_seqs <- rbind(household_seqs, haps_village[rowNum,])
    }
  }
  
  
  unique_ids_household <- household_seqs[,'ID']
  unique_ids_case <- case_seqs[,'ID']
  unique_ids_household <- unique (unique_ids_household)
  unique_ids_case <- unique (unique_ids_case)
  
  L1_D_sum =0
  L2_D_sum =0
  
  for (bootstrap in 1:100){

    random=sample(1:length(unique_ids_household), 1)
    ID1=unique_ids_case[1] #ID1 will be case
    ID2=unique_ids_household[random] #ID2 will be household random member
    
    ID1freqs = data.frame (A = numeric(288), C = numeric (288), G = numeric(288), Tn = numeric(288))
    ID2freqs = data.frame (A = numeric(288), C = numeric (288), G = numeric(288), Tn = numeric(288))
    
    
    for (nucleotide in 1:288){ # compute table of frequencies for ID1 and ID2
      for (rowNum in 1:nrow(case_seqs)){ #compute for ID1-case seqs first
        if (case_seqs$ID[rowNum]==ID1){
          seq=DNAString(case_seqs$haplotype[rowNum])
          bp = as.character(seq [nucleotide])
          if (bp == 'T'){
            ID1freqs$Tn[nucleotide] = ID1freqs$Tn[nucleotide] + case_seqs$freq[rowNum]
          }
          if (bp == 'A'){
            ID1freqs$A[nucleotide] = ID1freqs$A[nucleotide] + case_seqs$freq[rowNum]
          }
          if (bp == 'G'){
            ID1freqs$G[nucleotide] = ID1freqs$G[nucleotide] + case_seqs$freq[rowNum]
          }
          if (bp == 'C'){
            ID1freqs$C[nucleotide] = ID1freqs$C[nucleotide] + case_seqs$freq[rowNum]
          }
        }
      }
      for (rowNum in 1:nrow(household_seqs)){
        if (household_seqs$ID[rowNum]==ID2){
          seq=DNAString(household_seqs$haplotype[rowNum])
          bp = as.character(seq [nucleotide])
          if (bp == 'T'){
            ID2freqs$Tn[nucleotide] = ID2freqs$Tn[nucleotide] + household_seqs$freq[rowNum]
          }
          if (bp == 'A'){
            ID2freqs$A[nucleotide] = ID2freqs$A[nucleotide] + household_seqs$freq[rowNum]
          }
          if (bp == 'G'){
            ID2freqs$G[nucleotide] = ID2freqs$G[nucleotide] + household_seqs$freq[rowNum]
          }
          if (bp == 'C'){
            ID2freqs$C[nucleotide] = ID2freqs$C[nucleotide] + household_seqs$freq[rowNum]
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
  
  households$L1_D[iteration]=L1_D_sum/bootstrap
  households$L2_D[iteration]=L2_D_sum/bootstrap
  
}
  
write.csv(households, 'csp household L1L2norm casechild.csv')
##-----------------------------------




