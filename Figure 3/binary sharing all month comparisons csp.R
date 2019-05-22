#Create table of csp binary sharing between all months (randomly sampled individuals)

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


haps_village <- read.csv("csp_allhaps_demog_person.csv")
months = read.csv('csp_months.csv')

for (iteration in 1:nrow(months)){
  month_vector=double( length=nrow(months))
  
  for (cycle in 1:nrow(months)){
  
    fastaLines= c()
    for (i in 1:nrow(haps_village)){
      if (haps_village$month[i]==months$months[iteration]&& !is.na(haps_village$month[i])){
        fastaLines = c(fastaLines, as.character(paste(">",haps_village$month[i],"_", haps_village$ID[i], sep = "")))
        fastaLines = c(fastaLines,as.character(haps_village[i,"haplotype"]))
      }
    }
  
    for (rowNum in 1:nrow(haps_village)){
      if (haps_village[rowNum,"month"]==months$months[cycle] && !is.na(haps_village$month[rowNum])){
          fastaLines = c(fastaLines, as.character(paste(">",haps_village$month[rowNum],"_", haps_village$ID[rowNum], sep = "")))
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
    compute_haps
    
    haps_net = haploNet (compute_haps)
    
    ind_hap<-with (stack(setNames(attr(compute_haps, "index"), rownames(compute_haps))),table(hap=ind, individuals=rownames(haps)[values]))
    ind_hap  
    
    shared_count=0
    
    count=1
    query = substr(colnames(ind_hap)[1],1, 3)
    for (j in 2:ncol(ind_hap)){
      if (substr(colnames(ind_hap)[j],1, 3)==query){
        count=count+1
      }
    }
    
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
      
      shared = FALSE
      for (k in 1:nrow(ind_hap)){
        if ((ind_hap[k,random_v1]>=1)&&(ind_hap[k, random_v2]>=1)){
          shared=TRUE
        }
      }
      
      if (shared == TRUE){
        shared_count = shared_count+1
      }
    }
  
  
      shared_count
      shared_count_freq = shared_count/10000*100
      
      month_vector[cycle]=shared_count_freq
  }
  months=cbind(months, month_vector)
}
  
months=months[,5:ncol(months)]

write.csv(months, 'csp share_score_all_month_comparison.csv')

months=read.csv("csp share_score_all_month_comparison.csv")

#statistical test for same month vs different
months=months[,2:ncol(months)]
diag = nondiag = double()
for (i in 1:nrow(months)){
  for (j in 1:ncol(months)){

    if (i==j){
      diag = c(diag, months[i,j])
    }
    if (i>j){
      nondiag = c (nondiag, months[i,j])
    }
  }
}

median (diag)
median (nondiag)
wilcox.test(diag, nondiag, mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
            conf.int = FALSE, conf.level = 0.95)

#plot heat map
months[upper.tri(months, diag = FALSE)]=NA
Heatmap(as.matrix(months),col = colorRamp2(c(0, 45, 65),c("blue", "yellow", "red")), cluster_rows=FALSE,cluster_columns=FALSE, na_col="white")



