##Purpose: demographic matching of case households for csp

haps = read.csv('csp_allhaps_demog_person_age_decimaldate.csv')
households = read.csv('csp_case_households.csv')

household_match = data.frame(household=households$households, date=double(41), matched_household=character(41),stringsAsFactors = FALSE)

for (iteration in 1:nrow(household_match)){
  for (row in 1:nrow(haps)){
    if (household_match[iteration,'household']==as.character(haps[row,'household']) && haps[row,'person'] == 'case child'){
      household_match$date[iteration] = haps$date[row]
    }
  }
}

match_points = array (dim=c(41,41))

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

maximum = array(data=NA, 41)
match = array (data=NA, 41)
max_array = array(data=NA, dim=c(41,41))


#define maximum for each
maximum[1]=max(match_points[1,2:41])
maximum[41]=max(match_points[41,1:40])
for (i in 2:40){
  maximum[i]=max(match_points[i,1:(i-1)], match_points[i,(i+1):41])
}

maximum
#Pick random household with max matching score
for (i in 1:41){
  for (j in 1:41){
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


write.csv(household_match, 'csp matched households.csv')


