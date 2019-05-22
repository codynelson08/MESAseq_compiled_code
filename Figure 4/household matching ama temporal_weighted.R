##Purpose: demographic matching of case households for ama
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


haps = read.csv('ama_allhaps_demog_person_age_decimaldate.csv')
households = read.csv('ama_case_households.csv')

household_match = data.frame(household=households$households, date=double(44), matched_household=character(44),stringsAsFactors = FALSE)

for (iteration in 1:nrow(household_match)){
  for (row in 1:nrow(haps)){
    if (household_match[iteration,'household']==as.character(haps[row,'household']) && haps[row,'person'] == 'case child'){
      household_match$date[iteration] = haps$date[row]
    }
  }
}

match_points = array (dim=c(44,44))

for (iteration in 1:nrow(household_match)){
  
  for (comp_household in 1:nrow(household_match)){
    count = 0
    if (abs(household_match$date[iteration]-(household_match$date[comp_household]))<=(60/365)){
      count=1
    }
    if (abs(household_match$date[iteration]-(household_match$date[comp_household]))<=(30/365)){
      count=2
    }
    if (abs(household_match$date[iteration]-(household_match$date[comp_household]))<=(15/365)){
      count=3
    }
    if (abs(household_match$date[iteration]-(household_match$date[comp_household]))<=(10/365)){
      count=4
    }
    if (abs(household_match$date[iteration]-(household_match$date[comp_household]))<=(5/365)){
      count=5
    }
    match_points [iteration,comp_household] = count
  }
  
  
}
match_points

maximum = array(data=NA, 44)
match = array (data=NA, 44)
max_array = array(data=NA, dim=c(44,44))


#define maximum for each
maximum[1]=max(match_points[1,2:44])
maximum[44]=max(match_points[44,1:43])
for (i in 2:43){
  maximum[i]=max(match_points[i,1:(i-1)], match_points[i,(i+1):44])
}

maximum
#Pick random household with max matching score
for (i in 1:44){
  for (j in 1:44){
    if (match_points[i,j]==maximum[i]&& i!=j){
      max_array[i,j]=maximum[i]
    }
  }
  max_array
  locns = which(!is.na(max_array[i,]))
  random_num=sample(1:length(locns), 1, replace=FALSE)
  match[i]=locns[random_num]
}


match

for (i in 1:nrow(household_match)){
  household_match$matched_household[i] = as.character(household_match$household[match[i]])
}


write.csv(household_match, 'ama matched households time_control.csv')


