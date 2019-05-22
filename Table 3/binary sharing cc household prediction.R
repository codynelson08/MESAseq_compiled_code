library(som)
library(ComplexHeatmap)
library(circlize)
#program to identify the maximum binary share score for each household

csp_matrix = read.csv ("csp share score case_child matrix edited.csv")
ama_matrix = read.csv ("ama share score case_child matrix edited.csv")
households = colnames(ama_matrix)
households = households [2:39]

csp_matrix_norm = normalize(as.matrix(csp_matrix[,2:ncol(csp_matrix)]))
ama_matrix_norm = normalize(as.matrix(ama_matrix[,2:ncol(ama_matrix)]))

combined_matrix = as.data.frame(as.matrix(csp_matrix_norm) + as.matrix(ama_matrix_norm))

sort_matrix = data.frame()

for (i in 1:length(households)){ #fill table with ranks
  ranks = rev(order(rank(combined_matrix[i,], ties.method = "random")))
  sort_matrix = rbind(sort_matrix, ranks)
}
colnames(sort_matrix) = 1:38
sort_matrix_ids = sort_matrix

for (j in 1:length(households)){ #swap ranks with household ids
  for(k in 1:length(households)){
    sort_matrix_ids[j,k] = households[sort_matrix[j,k]]
  }
}


##---------------------------

allhaps = read.csv('csp_allhaps_demog_person_age_decimaldate_latlong.csv')

household_df = data.frame()
district_df = data.frame()
time_df = data.frame()

for (j in 1:nrow(sort_matrix_ids)){
  
  case_child = households[j]
  ranked_households = as.vector(sort_matrix_ids[j,])
  household_vector = district_vector = time_vector = integer(38)
  
  for (k in 1:length(ranked_households)){
    if (case_child == as.character(ranked_households[k])){
      household_vector [k] = 1
    }
    cc_lat = comp_cc_lat = double()
    cc_long = comp_cc_long = double()
    cc_time = comp_cc_time = double()
    
    for (l in 1:nrow(allhaps)){
      if (case_child==allhaps$household[l]){
        cc_lat = allhaps$latitude[l]
        cc_long = allhaps$longitude[l]
        cc_time = allhaps$date[l]
      }
      if (as.character(ranked_households[k])==allhaps$household[l]){
        comp_cc_lat = allhaps$latitude[l]
        comp_cc_long = allhaps$longitude[l]
        comp_cc_time = allhaps$date[l]
      }
    }
    
    distance = sqrt((cc_lat-comp_cc_lat)*(cc_lat-comp_cc_lat)+(cc_long-comp_cc_long)*(cc_long-comp_cc_long))
    
    if(distance<0.0405){
      district_vector[k] = 4.5
    }
    if(distance<0.0202){
      district_vector[k] = 2.25
    }
    if(distance<0.00675){
      district_vector[k] = 0.75
    }
   
   
    if(abs(cc_time-comp_cc_time)<(60/365)){
      time_vector [k] = 60
    }
    if(abs(cc_time-comp_cc_time)<(30/365)){
      time_vector [k] = 30
    }
    if(abs(cc_time-comp_cc_time)<(10/365)){
      time_vector [k] = 10
    }
    
  }
  
  
  household_df = rbind(household_df, household_vector)
  district_df = rbind (district_df, district_vector)
  time_df = rbind (time_df, time_vector)
}

write.csv(district_df, "binary geog sum.csv")
write.csv(time_df, "binary time sum.csv")
group_by(district_df)




mat = matrix (c(13, 25, 2, 36), nrow=2)
fisher.test(mat)
