#Fisher p analysis of haplotype difference in cc and household for csp

#### --------- load packages ----------------- ####
library(ggplot2)
library(readr)
library(dplyr)
library(gridExtra) 
library(ggthemes)
library(wesanderson)
library(tidyr)
library(lubridate)


#### ----- creating the figure of haplotype presence over time ----- ####

## for CSP


# # read in the csp haplotype merged data set you made for Cody, only keep the first 125 columns
csp_merged = read.csv('csp_haps_demog_merged_nocontrols.csv')
csp_merged = csp_merged[,c(2:126,250)]

# create a new variable that is just the month
str(csp_merged$interview_date)
csp_merged$interview_month = NA
tmp=dmy(csp_merged$interview_date)
csp_merged$interview_month = paste0(month(tmp),"-",year(tmp))
write.csv(csp_merged, 'csp_haps_demog_merged_nocontrols_month.csv')

#convert to number of reads
csp_merged_freqs = csp_merged
for (i in 1:nrow(csp_merged)){
  for (j in 1:120){
    csp_merged[i,j]=as.integer(csp_merged_freqs[i,j]*csp_merged_freqs$totalcount[i])
  }
}

# only keep the first 120 columns and month
csp_merged = csp_merged[,c(1:120,126:127)]
colnames (csp_merged) <- c(as.integer(1:120),"person",'interview_month')

csp_case_child = data.frame()
csp_household = data.frame()

for (i in 1:nrow(csp_merged)){
  if(csp_merged$person[i]=="case household member"&&!is.na(csp_merged$person[i])){
    csp_household=rbind(csp_household, csp_merged[i,])
  }
  if(csp_merged$person[i]=="case child"&&!is.na(csp_merged$person[i])){
    csp_case_child=rbind(csp_case_child, csp_merged[i,])
  }
}

# create a data frame summarizing each haplotype and the months it is present
# trying gathering the code to long format
long_csp_merged_case_child = long_csp_merged_household = data.frame()
long_csp_merged_case_child = gather(csp_case_child, H, reads, -interview_month) 
long_csp_merged_household = gather(csp_household, H, reads, -interview_month) 

# rename columns in csp_merged long file
names(long_csp_merged_case_child)[names(long_csp_merged_case_child) == 'reads'] <- 'reads_present'
names(long_csp_merged_case_child)[names(long_csp_merged_case_child) == 'H'] <- 'haplotype'
names(long_csp_merged_household)[names(long_csp_merged_household) == 'reads'] <- 'reads_present'
names(long_csp_merged_household)[names(long_csp_merged_household) == 'H'] <- 'haplotype'
  
# remove all rows with reads_present equal to 0
long_csp_merged_case_child = long_csp_merged_case_child[-which(long_csp_merged_case_child$reads_present == 0),]
long_csp_merged_case_child = long_csp_merged_case_child[-which(long_csp_merged_case_child$haplotype == "person"),]
long_csp_merged_household = long_csp_merged_household[-which(long_csp_merged_household$reads_present == 0),]
long_csp_merged_household = long_csp_merged_household[-which(long_csp_merged_household$haplotype == "person"),]

# summarize the new data set by month
month_summary_case_child = long_csp_merged_case_child %>% 
  group_by(interview_month,haplotype) %>%
  summarize(case_child_n=n())
month_summary_household = long_csp_merged_household %>% 
  group_by(interview_month,haplotype) %>%
  summarize(household_n=n())
 
# set order for x-axis for months
table(month_summary_case_child$interview_month, useNA = "always")
month_order = c("4-2013","5-2013","6-2013","7-2013","8-2013","9-2013","10-2013","11-2013","12-2013","1-2014","2-2014","3-2014",
                "4-2014","5-2014","6-2014")
month_summary_case_child$haplotype <- as.numeric(as.character(month_summary_case_child$haplotype))
month_summary_case_child = month_summary_case_child[  with(month_summary_case_child, order(haplotype, interview_month)), ]
month_summary_household$haplotype <- as.numeric(as.character(month_summary_household$haplotype))
month_summary_household = month_summary_household[  with(month_summary_household, order(haplotype, interview_month)), ]

# create a plot of the presence and abundance of each haplotype over time
month_summary_case_child$interview_month = factor (month_summary_case_child$interview_month, levels=month_order)
month_summary_case_child=month_summary_case_child[order(month_summary_case_child$interview_month),]
#month_summary_case_child$case_child_n = (month_summary_case_child$case_child_n)/283

month_summary_household$interview_month = factor (month_summary_household$interview_month, levels=month_order)
month_summary_household=month_summary_household[order(month_summary_household$interview_month),]
#month_summary_household$household_n = (month_summary_household$household_n)/253

month_summary_merge <- month_summary_case_child %>% full_join(month_summary_household, by=c("interview_month","haplotype"))
month_summary_merge[is.na(month_summary_merge)] <- 0
month_summary_merge$n_difference = month_summary_merge$case_child_n/283 - month_summary_merge$household_n/253
month_summary_merge$FisherP = NA

for (i in 1:nrow(month_summary_merge)){ #compute fisher p for each value
  mat = matrix(c(month_summary_merge$case_child_n[i], month_summary_merge$household_n[i],283-month_summary_merge$case_child_n[i], 253-month_summary_merge$household_n[i]),
             nrow = 2)
  testresult=fisher.test(mat)
  month_summary_merge$FisherP[i] = -1*log10(as.double(testresult[1]))
 
}

csp_month_plot = ggplot(month_summary_merge, aes(as.numeric(x=interview_month), y=FisherP, color = n_difference)) +
 #geom_line()+
  geom_jitter(height=0.007, width=0.2)+
 # geom_smooth(method = loess, span=10, size = 0, level = 0.99)+
  theme_bw()+
  labs (color="PDH")+
  scale_color_gradient(low="blue", high="red")+
  theme(axis.text.x = element_text(angle = 90, hjust=1))
 
csp_month_plot

ggsave(csp_month_plot, filename="relprev graph pvalue ci csp.png", device="png",
       height=4, width=6, units="in", dpi=500)

