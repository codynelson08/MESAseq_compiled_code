#Visualize csp haplotype presence over time (Figure 2B)

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

# # read in the csp haplotype merged data set you made for Cody, only keep the first 125 columns
csp_merged = read.csv('csp_haps_demog_merged_nocontrols.csv')
csp_merged = csp_merged[,2:126]

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
csp_merged = csp_merged[,c(1:120,126)]
colnames (csp_merged) <- c(as.integer(1:120),'interview_month')

# create a data frame summarizing each haplotype and the months it is present
# trying gathering the code to long format
long_csp_merged = gather(csp_merged, H, reads, -interview_month) 

# rename columns in csp_merged long file
names(long_csp_merged)[names(long_csp_merged) == 'reads'] <- 'reads_present'
names(long_csp_merged)[names(long_csp_merged) == 'H'] <- 'haplotype'
  
# remove all rows with reads_present equal to 0
long_csp_merged = long_csp_merged[-which(long_csp_merged$reads_present == 0),]

# summarize the new data set by month
month_summary = long_csp_merged %>% 
  group_by(interview_month,haplotype) %>%
  summarize(n=n())

month_summary$total_reads = 0

for (i in 1:nrow(month_summary)){
  for (j in 1:nrow(long_csp_merged)){
    if (month_summary$interview_month[i]==long_csp_merged$interview_month[j] && month_summary$haplotype[i]==long_csp_merged$haplotype[j]){
      month_summary$total_reads[i] = month_summary$total_reads [i] + long_csp_merged$reads_present[j]
    }
  }
}
  
# set order for x-axis for months
table(month_summary$interview_month, useNA = "always")
month_order = c("4-2013","5-2013","6-2013","7-2013","8-2013","9-2013","10-2013","11-2013","12-2013","1-2014","2-2014","3-2014",
                "4-2014","5-2014","6-2014")
month_summary <- within(month_summary, interview_month <- factor(interview_month, levels=month_order))
month_summary$haplotype <- as.integer(as.character(month_summary$haplotype))


case_children = read.csv('csp_case_children.csv')
case_children
case_children = merge(case_children)
x=table(csp_merged$interview_month)


# create a plot of the presence and abundance of each haplotype over time
csp_month_plot = ggplot(month_summary, aes(x=interview_month, y=haplotype, size=n, color=total_reads)) +
  geom_point() +
  labs(title = "Pfcsp haplotype presence and abundance over time across MESA samples (1/1/18 sequencing run)",
       y = "Month and year",x="Haplotype", color = "Total number reads", size = "Number of participants")+
  scale_y_reverse(breaks = c(1:120))+
  theme_bw()+
  #scale_x_continuous(breaks=c(2013.4, 2014.6, 0.1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

csp_month_plot

ggsave(csp_month_plot, filename="csp_month_plot.png", device="png",
       height=20, width=6, units="in", dpi=500)

