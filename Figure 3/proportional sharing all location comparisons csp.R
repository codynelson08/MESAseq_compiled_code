#Create table of csp proportional sharing between locations 
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
library(ComplexHeatmap)
library(circlize)


haps_village <- read.csv("csp_allhaps_demog_person_age_decimaldate_latlong.csv")
districts = read.csv('csp location subset.csv')


for (iteration in 1:nrow(districts)){
  village_vector=double( length=nrow(districts))

  for (cycle in 1:nrow(districts)){
  
    fastaLines= c()
    for (i in 1:nrow(haps_village)){
      if (haps_village$location[i]==districts$locations[iteration]&& !is.na(haps_village$location[i])){
        if (haps_village$date[i]>2013.293 &&haps_village$date[i]<2013.488){
          fastaLines = c(fastaLines, as.character(paste(">",haps_village$location[i],"_", haps_village$ID[i], sep = "")))
          fastaLines = c(fastaLines,as.character(haps_village[i,"haplotype"]))
        }
      }
    }
  
    for (rowNum in 1:nrow(haps_village)){
      if (haps_village[rowNum,"location"]==districts$locations[cycle] && !is.na(haps_village$location[rowNum])){
        if (haps_village$date[i]>2013.293 &&haps_village$date[i]<2013.488){
          fastaLines = c(fastaLines, as.character(paste(">",haps_village$location[rowNum],"_", haps_village$ID[rowNum], sep = "")))
          fastaLines = c(fastaLines,as.character(haps_village[rowNum,"haplotype"]))
        }
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
    
    
    count=1
    query = substr(colnames(ind_hap)[1],1, 3)
    for (j in 2:ncol(ind_hap)){
      if (substr(colnames(ind_hap)[j],1, 3)==query){
        count=count+1
      }
    }
    
    freq_sum=0
    
    for (bootstrap in 1:10000){
      
      random_v1=sample(1:count, 1, replace=FALSE)
      random_v2=sample((count+1):ncol(ind_hap), 1, replace=FALSE)
      
      if (count>=ncol(ind_hap)){
        random = sample(1:ncol(ind_hap),2, replace=FALSE)
        random_v1=random[1]
        random_v2=random[2]
      }
      
     # ind_hap = ind_hap[,""]
     # ind_hap[,random_v1]
      #ind_hap[,random_v2]
      
      shared=total=0
      for (k in 1:nrow(ind_hap)){
        if ((ind_hap[k,random_v1]>=1)&&(ind_hap[k, random_v2]==0)){
          total = total+1
        }
        if ((ind_hap[k,random_v1]==0)&&(ind_hap[k, random_v2]>=1)){
          total = total+1
        }
        if ((ind_hap[k,random_v1]>=1)&&(ind_hap[k, random_v2]>=1)){
          shared=shared+1
          total = total +1
        }
      }
      
      freq_sum = freq_sum+ (shared/total)
    }
  
  
      
      weighted_share_score = freq_sum/10000*100
      
      village_vector[cycle]=weighted_share_score
  }
  districts=cbind(districts, village_vector)
}

districts=districts[,2:ncol(districts)]

write.csv(districts, 'weighted_share_score_all_village_comparison_2013.csv')
district_mat=read.csv("weighted_share_score_all_village_comparison_2013.csv")
district_mat=district_mat[,2:6]
diag = nondiag = double()
for (i in 1:nrow(district_mat)){
  for (j in 1:ncol(district_mat)){
    
    if (i==j){
      diag = c(diag, district_mat[i,j])
    }
    if (i>j){
      nondiag = c (nondiag, district_mat[i,j])
    }
  }
}

median (diag)
median (nondiag)
wilcox.test(diag, nondiag, mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
            conf.int = FALSE, conf.level = 0.95)

#districts = read.csv('csp_locations.csv')
colnames(district_mat)=districts$locations
rownames(district_mat)=districts$locations
district_mat[upper.tri(district_mat, diag = FALSE)]=NA
Heatmap(as.matrix(district_mat), col = colorRamp2(c(0, 45, 65),c("blue", "yellow", "red")),cluster_rows=FALSE,cluster_columns=FALSE, na_col="white")





