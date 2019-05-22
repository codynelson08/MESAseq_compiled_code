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
case_childrenL1 = case_childrenL2 = read.csv('csp_casechild_share.csv')

for (case_iteration in 24:nrow(households)){
  tablecount = 2
  print(paste0(Sys.time(), " --> starting case_child ", case_iteration, " of ", nrow(households) ))
  cat("\n")
  
  for (household_iteration in 1:nrow(households)){
    
    matched_household_seqs = data.frame(haplotype=character(), location=character(), village=character(),freq=double(), ID=character(), household=character(), month=character())
    case_seqs = data.frame(haplotype=character(), location=character(), village=character(),freq=double(), ID=character(), household=character(), month=character())
    print(paste0(Sys.time(), " --> starting household ", household_iteration, " of ", nrow(households) ))
    print(paste0("Case_child is ", households$households[case_iteration], " and household is ", households$households[household_iteration]) )
    cat("\n")
    
    for (rowNum in 1:nrow(haps_village)){
      if (haps_village[rowNum,"household"]==households$households[case_iteration] && haps_village[rowNum,"person"]=='case child'){
        case_seqs <- rbind(case_seqs, haps_village[rowNum,])
      }
      if (haps_village[rowNum,"household"]==households$households[household_iteration] && haps_village[rowNum,"person"]=='case household member'){
        matched_household_seqs <- rbind(matched_household_seqs, haps_village[rowNum,])
      }
    }
    
    
    unique_ids_matched_household <- matched_household_seqs[,'ID']
    unique_ids_case <- case_seqs[,'ID']
    unique_ids_matched_household <- unique (unique_ids_matched_household)
    unique_ids_case <- unique (unique_ids_case)
    
    L1_D_sum =0
    L2_D_sum =0
    
    for (bootstrap in 1:100){
      
      random=sample(1:length(unique_ids_matched_household), 1)
      ID1=unique_ids_case[1] #ID1 will be case
      ID2=unique_ids_matched_household[random] #ID2 will be matched_household random member
      
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
        for (rowNum in 1:nrow(matched_household_seqs)){
          if (matched_household_seqs$ID[rowNum]==ID2){
            seq=DNAString(matched_household_seqs$haplotype[rowNum])
            bp = as.character(seq [nucleotide])
            if (bp == 'T'){
              ID2freqs$Tn[nucleotide] = ID2freqs$Tn[nucleotide] + matched_household_seqs$freq[rowNum]
            }
            if (bp == 'A'){
              ID2freqs$A[nucleotide] = ID2freqs$A[nucleotide] + matched_household_seqs$freq[rowNum]
            }
            if (bp == 'G'){
              ID2freqs$G[nucleotide] = ID2freqs$G[nucleotide] + matched_household_seqs$freq[rowNum]
            }
            if (bp == 'C'){
              ID2freqs$C[nucleotide] = ID2freqs$C[nucleotide] + matched_household_seqs$freq[rowNum]
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
    
    case_childrenL1[case_iteration,tablecount] = L1_D_sum/bootstrap
    case_childrenL2[case_iteration,tablecount] = L2_D_sum/bootstrap
    tablecount=tablecount+1
  }
}

write.csv(case_childrenL1, 'csp household L1norm casechild matrix.csv')
write.csv(case_childrenL2, 'csp household L2norm casechild matrix.csv')
##-----------------------------------




