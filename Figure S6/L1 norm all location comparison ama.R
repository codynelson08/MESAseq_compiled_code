#Create table of ama L1 norm between locations 
#within a specific time window (randomly sampled individuals)

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



haps_village = read.csv('ama_allhaps_demog_person_age_decimaldate_latlong.csv')
districts = read.csv('ama location subset.csv')

for (iteration in 1:nrow(districts)){
  tablecount = 2
  district_vector=double( length=nrow(districts))
  print(paste0(Sys.time(), " --> starting district ", iteration, " of ", nrow(districts) ))
  
  for (cycle in 1:nrow(districts)){
    print(paste0(Sys.time(), " --> starting cycle ", cycle, " of ", nrow(districts) ))
    district1_seqs = data.frame(haplotype=character(), location=character(), village=character(),freq=double(), ID=character(), household=character(), district=character())
    district2_seqs = data.frame(haplotype=character(), location=character(), village=character(),freq=double(), ID=character(), household=character(), district=character())
    
    for (rowNum in 1:nrow(haps_village)){
      if (haps_village$location[rowNum]==districts$locations[iteration]&& !is.na(haps_village$location[rowNum])){
        if (haps_village$date[rowNum]>2013.293 &&haps_village$date[rowNum]<2013.488){
          district1_seqs <- rbind(district1_seqs, haps_village[rowNum,])
        }
      }
      if (haps_village$location[rowNum]==districts$locations[cycle]&& !is.na(haps_village$location[rowNum])){
        if (haps_village$date[rowNum]>2013.293 &&haps_village$date[rowNum]<2013.488){
          district2_seqs <- rbind(district2_seqs, haps_village[rowNum,])
        }
      }
    }
    
    
    unique_ids_district1 <- district2_seqs[,'ID']
    unique_ids_district2 <- district1_seqs[,'ID']
    unique_ids_district1 <- unique (unique_ids_district1)
    unique_ids_district2 <- unique (unique_ids_district2)
    
    L1_D_sum =0
    L2_D_sum =0
    
    for (bootstrap in 1:100){
  
      random=sample(1:length(unique_ids_district1), 1)
      ID1=unique_ids_district2[1] #ID1 will be case
      ID2=unique_ids_district1[random] #ID2 will be matched_household random member
      
      ID1freqs = data.frame (A = numeric(300), C = numeric (300), G = numeric(300), Tn = numeric(300))
      ID2freqs = data.frame (A = numeric(300), C = numeric (300), G = numeric(300), Tn = numeric(300))
      
      
      for (nucleotide in 1:300){ # compute table of frequencies for ID1 and ID2
        for (rowNum in 1:nrow(district1_seqs)){ #compute for ID1-case seqs first
          if (district1_seqs$ID[rowNum]==ID1){
            seq=DNAString(district1_seqs$haplotype[rowNum])
            bp = as.character(seq [nucleotide])
            if (bp == 'T'){
              ID1freqs$Tn[nucleotide] = ID1freqs$Tn[nucleotide] + district1_seqs$freq[rowNum]
            }
            if (bp == 'A'){
              ID1freqs$A[nucleotide] = ID1freqs$A[nucleotide] + district1_seqs$freq[rowNum]
            }
            if (bp == 'G'){
              ID1freqs$G[nucleotide] = ID1freqs$G[nucleotide] + district1_seqs$freq[rowNum]
            }
            if (bp == 'C'){
              ID1freqs$C[nucleotide] = ID1freqs$C[nucleotide] + district1_seqs$freq[rowNum]
            }
          }
        }
        for (rowNum in 1:nrow(district2_seqs)){
          if (district2_seqs$ID[rowNum]==ID2){
            seq=DNAString(district2_seqs$haplotype[rowNum])
            bp = as.character(seq [nucleotide])
            if (bp == 'T'){
              ID2freqs$Tn[nucleotide] = ID2freqs$Tn[nucleotide] + district2_seqs$freq[rowNum]
            }
            if (bp == 'A'){
              ID2freqs$A[nucleotide] = ID2freqs$A[nucleotide] + district2_seqs$freq[rowNum]
            }
            if (bp == 'G'){
              ID2freqs$G[nucleotide] = ID2freqs$G[nucleotide] + district2_seqs$freq[rowNum]
            }
            if (bp == 'C'){
              ID2freqs$C[nucleotide] = ID2freqs$C[nucleotide] + district2_seqs$freq[rowNum]
            }
          }
        }
      }
      
      Dk = data.frame (A = numeric(300), C = numeric (300), G = numeric(300), Tn = numeric(300), L1Dk = numeric (300), L2Dk = numeric (300))
      for (nucleotide in 1:300){ #Compute L1 and L2 Dk for each nucleotide position
        Dk$A[nucleotide] = abs(ID1freqs$A[nucleotide]-ID2freqs$A[nucleotide])
        Dk$C[nucleotide] = abs(ID1freqs$C[nucleotide]-ID2freqs$C[nucleotide])
        Dk$G[nucleotide] = abs(ID1freqs$G[nucleotide]-ID2freqs$G[nucleotide])
        Dk$Tn[nucleotide] = abs(ID1freqs$Tn[nucleotide]-ID2freqs$Tn[nucleotide])
        Dk$L1Dk[nucleotide] = Dk$A[nucleotide] + Dk$C[nucleotide] + Dk$G[nucleotide] + Dk$Tn[nucleotide]
        #Dk$L2Dk[nucleotide] = sqrt((Dk$A[nucleotide]*Dk$A[nucleotide])+(Dk$C[nucleotide]*Dk$C[nucleotide])+(Dk$G[nucleotide]*Dk$G[nucleotide])+(Dk$T[nucleotide]*Dk$T[nucleotide]))
      }
      
      L1_D = L2_D = 0
      
      for (nucleotide in 1:300){ #compute L1 and L2 statistics
        L1_D = L1_D + Dk$L1Dk[nucleotide]
       # L2_D = L2_D + Dk$L2Dk[nucleotide]
      }
      
      L1_D_sum = L1_D_sum + L1_D
     # L2_D_sum = L2_D_sum + L2_D
    }
    
    
    district_vector[cycle]=L1_D_sum/bootstrap
  }
  districts=cbind(districts, district_vector)
}
  
districts=districts[,2:ncol(districts)]

write.csv(districts, 'L1_2013_district_comparison.csv')

districts=read.csv("L1_2013_district_comparison.csv")
districts=districts[,2:ncol(districts)]
diag = nondiag = double()
for (i in 1:nrow(districts)){
  for (j in 1:ncol(districts)){
    
    if (i==j){
      diag = c(diag, districts[i,j])
    }
    if (i>j){
      nondiag = c (nondiag, districts[i,j])
    }
  }
}

median (diag)
median (nondiag)
wilcox.test(diag, nondiag, mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
            conf.int = FALSE, conf.level = 0.95)


median(districts)
districts[upper.tri(districts, diag = FALSE)]=NA
Heatmap(as.matrix(districts),col = colorRamp2(c(26, 18, 14),c("blue", "yellow", "red")), cluster_rows=FALSE,cluster_columns=FALSE, na_col="white")

##-----------------------------------




