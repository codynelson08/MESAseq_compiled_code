#Create table of csp L1 norm between all months (randomly sampled individuals)


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
library(ComplexHeatmap)
library(circlize)

haps_village = read.csv('csp_allhaps_demog_person.csv')
months = read.csv('csp_months.csv')

for (iteration in 5:nrow(months)){
  tablecount = 2
  month_vector=double( length=nrow(months))
  print(paste0(Sys.time(), " --> starting month ", iteration, " of ", nrow(months) ))
  
  for (cycle in 1:nrow(months)){
    print(paste0(Sys.time(), " --> starting cycle ", cycle, " of ", nrow(months) ))
    month1_seqs = data.frame(haplotype=character(), location=character(), village=character(),freq=double(), ID=character(), household=character(), month=character())
    month2_seqs = data.frame(haplotype=character(), location=character(), village=character(),freq=double(), ID=character(), household=character(), month=character())
    
    for (rowNum in 1:nrow(haps_village)){
      if (haps_village[rowNum,"month"]==months$months[iteration]){
          month1_seqs <- rbind(month1_seqs, haps_village[rowNum,])
      }
      if (haps_village[rowNum,"month"]==months$months[cycle]){
          month2_seqs <- rbind(month2_seqs, haps_village[rowNum,])
      }
    }
    
    
    unique_ids_month1 <- month2_seqs[,'ID']
    unique_ids_month2 <- month1_seqs[,'ID']
    unique_ids_month1 <- unique (unique_ids_month1)
    unique_ids_month2 <- unique (unique_ids_month2)
    
    L1_D_sum =0
    L2_D_sum =0
    
    for (bootstrap in 1:100){
  
      random=sample(1:length(unique_ids_month1), 1)
      ID1=unique_ids_month2[1] #ID1 will be case
      ID2=unique_ids_month1[random] #ID2 will be matched_household random member
      
      ID1freqs = data.frame (A = numeric(288), C = numeric (288), G = numeric(288), Tn = numeric(288))
      ID2freqs = data.frame (A = numeric(288), C = numeric (288), G = numeric(288), Tn = numeric(288))
      
      
      for (nucleotide in 1:288){ # compute table of frequencies for ID1 and ID2
        for (rowNum in 1:nrow(month1_seqs)){ #compute for ID1-case seqs first
          if (month1_seqs$ID[rowNum]==ID1){
            seq=DNAString(month1_seqs$haplotype[rowNum])
            bp = as.character(seq [nucleotide])
            if (bp == 'T'){
              ID1freqs$Tn[nucleotide] = ID1freqs$Tn[nucleotide] + month1_seqs$freq[rowNum]
            }
            if (bp == 'A'){
              ID1freqs$A[nucleotide] = ID1freqs$A[nucleotide] + month1_seqs$freq[rowNum]
            }
            if (bp == 'G'){
              ID1freqs$G[nucleotide] = ID1freqs$G[nucleotide] + month1_seqs$freq[rowNum]
            }
            if (bp == 'C'){
              ID1freqs$C[nucleotide] = ID1freqs$C[nucleotide] + month1_seqs$freq[rowNum]
            }
          }
        }
        for (rowNum in 1:nrow(month2_seqs)){
          if (month2_seqs$ID[rowNum]==ID2){
            seq=DNAString(month2_seqs$haplotype[rowNum])
            bp = as.character(seq [nucleotide])
            if (bp == 'T'){
              ID2freqs$Tn[nucleotide] = ID2freqs$Tn[nucleotide] + month2_seqs$freq[rowNum]
            }
            if (bp == 'A'){
              ID2freqs$A[nucleotide] = ID2freqs$A[nucleotide] + month2_seqs$freq[rowNum]
            }
            if (bp == 'G'){
              ID2freqs$G[nucleotide] = ID2freqs$G[nucleotide] + month2_seqs$freq[rowNum]
            }
            if (bp == 'C'){
              ID2freqs$C[nucleotide] = ID2freqs$C[nucleotide] + month2_seqs$freq[rowNum]
            }
          }
        }
      }
      
      Dk = data.frame (A = numeric(288), C = numeric (288), G = numeric(288), Tn = numeric(288), L1Dk = numeric (288), L2Dk = numeric (288))
      for (nucleotide in 1:288){ #Compute L1 and L2 Dk for each nucleotide position
        Dk$A[nucleotide] = abs(ID1freqs$A[nucleotide]-ID2freqs$A[nucleotide])
        Dk$C[nucleotide] = abs(ID1freqs$C[nucleotide]-ID2freqs$C[nucleotide])
        Dk$G[nucleotide] = abs(ID1freqs$G[nucleotide]-ID2freqs$G[nucleotide])
        Dk$Tn[nucleotide] = abs(ID1freqs$Tn[nucleotide]-ID2freqs$Tn[nucleotide])
        Dk$L1Dk[nucleotide] = Dk$A[nucleotide] + Dk$C[nucleotide] + Dk$G[nucleotide] + Dk$Tn[nucleotide]
        #Dk$L2Dk[nucleotide] = sqrt((Dk$A[nucleotide]*Dk$A[nucleotide])+(Dk$C[nucleotide]*Dk$C[nucleotide])+(Dk$G[nucleotide]*Dk$G[nucleotide])+(Dk$T[nucleotide]*Dk$T[nucleotide]))
      }
      
      L1_D = L2_D = 0
      
      for (nucleotide in 1:288){ #compute L1 and L2 statistics
        L1_D = L1_D + Dk$L1Dk[nucleotide]
       # L2_D = L2_D + Dk$L2Dk[nucleotide]
      }
      
      L1_D_sum = L1_D_sum + L1_D
     # L2_D_sum = L2_D_sum + L2_D
    }
    
    
    month_vector[cycle]=L1_D_sum/bootstrap
  }
  months=cbind(months, month_vector)
}
  
months=months[,5:ncol(months)]

write.csv(months, 'L1_all_month_comparison.csv')

months=read.csv("L1_all_month_comparison.csv")
months=months[,2:ncol(months)]


months[upper.tri(months, diag = FALSE)]=NA
Heatmap(as.matrix(months),col = colorRamp2(c(18, 7, 2),c("blue", "yellow", "red")), cluster_rows=FALSE,cluster_columns=FALSE, na_col="white")

##-----------------------------------




