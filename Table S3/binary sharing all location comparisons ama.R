#Create table of ama binary sharing between locations 
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

# read in the ama haplotype merged data set you made for Cody
#ama_data <- read.csv('ama_haps_demog_merged_nocontrols.csv', header = T)

# read in the merged data set
#merged_data <- read.csv('mesa_data_clean (1).csv', header = T)

# merge in the interview dates from the merged_data set
#ama_merged = left_join(ama_data,merged_data,by="labid")

# only keep the first 95 columns
#ama_merged = ama_merged[,c(1:95,155)]

# create a new variable that is just the month
str(ama_merged$interview_date)
ama_merged$interview_date = mdy(ama_merged$interview_date)
str(ama_merged$interview_date)
ama_merged$month = paste0(month(ama_merged$interview_date),"-",year(ama_merged$interview_date))
table(ama_merged$month, useNA = "always")

haps_village <- data.frame (haplotype=character(), village=character(),freq=double(), ID=character(), household=character(), month=character()) 
for(i in 1:nrow(ama_merged)){
  for(j in 1:(ncol(ama_merged)-9)){
    if (!is.na(ama_merged[i,j])){
      temp <- data.frame (colnames(ama_merged[j]), ama_merged[i,"location.x"], ama_merged[i,j], ama_merged[i,'labid'], ama_merged[i,'studyid_case_control_new'], ama_merged[i,'month'])
      haps_village <- rbind(haps_village,temp)
    }}}


colnames (haps_village) <- c('haplotype', 'village', 'freq', 'ID', 'household', 'month')
#----------------------------------------------------
haps_village <- read.csv("ama_allhaps_demog_person_age_decimaldate_latlong.csv")
districts = read.csv('ama location subset.csv')


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
      
      village_vector[cycle]=shared_count_freq
  }
  districts=cbind(districts, village_vector)
}

districts=districts[,2:ncol(districts)]

write.csv(districts, 'ama share_score_all_village_comparison_2013.csv')
district_mat=read.csv("ama share_score_all_village_comparison_2013.csv")
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

#districts = read.csv('ama_locations.csv')
colnames(district_mat)=districts$locations
rownames(district_mat)=districts$locations
district_mat[upper.tri(district_mat, diag = FALSE)]=NA
Heatmap(as.matrix(district_mat), col = colorRamp2(c(0, 30, 70),c("blue", "yellow", "red")),cluster_rows=FALSE,cluster_columns=FALSE, na_col="white")






