#Visualize ama1 haplotype presence over time (Figure S5B)

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

## for ama


# # read in the ama haplotype merged data set you made for Cody, only keep the first 125 columns
ama_merged = read.csv('ama_haps_demog_merged_nocontrols_censored.csv')
ama_merged = ama_merged[,2:186]

# create a new variable that is just the month
str(ama_merged$interview_date)
ama_merged$interview_month = NA
tmp=dmy(ama_merged$interview_date)
ama_merged$interview_month = paste0(month(tmp),"-",year(tmp))
write.csv(ama_merged, 'ama_haps_demog_merged_nocontrols_month.csv')

#convert to number of reads
ama_merged_freqs = ama_merged
for (i in 1:nrow(ama_merged)){
  for (j in 1:180){
    ama_merged[i,j]=as.integer(ama_merged_freqs[i,j]*ama_merged_freqs$totalcount[i])
  }
}

# only keep the first 180 columns and month
ama_merged = ama_merged[,c(1:180,186)]
colnames (ama_merged) <- c(as.integer(1:180),'interview_month')

# create a data frame summarizing each haplotype and the months it is present
# trying gathering the code to long format
long_ama_merged = gather(ama_merged, H, reads, -interview_month) 

# rename columns in ama_merged long file
names(long_ama_merged)[names(long_ama_merged) == 'reads'] <- 'reads_present'
names(long_ama_merged)[names(long_ama_merged) == 'H'] <- 'haplotype'
  
# remove all rows with reads_present equal to 0
long_ama_merged = long_ama_merged[-which(long_ama_merged$reads_present == 0),]

# summarize the new data set by month
month_summary = long_ama_merged %>% 
  group_by(interview_month,haplotype) %>%
  summarize(n=n())
  
month_summary$total_reads = 0

for (i in 1:nrow(month_summary)){
  for (j in 1:nrow(long_ama_merged)){
    if (month_summary$interview_month[i]==long_ama_merged$interview_month[j] && month_summary$haplotype[i]==long_ama_merged$haplotype[j]){
      month_summary$total_reads[i] = month_summary$total_reads [i] + long_ama_merged$reads_present[j]
    }
  }
}

# check the output
length(which(ama_merged$interview_month == "1-2014" & ama_merged$1 > 0))
length(which(ama_merged$interview_month == "1-2014" & ama_merged$10 > 0))
length(which(ama_merged$interview_month == "1-2014" & ama_merged$11 > 0))

# set order for x-axis for months
table(month_summary$interview_month, useNA = "always")
month_order = c("4-2013","5-2013","6-2013","7-2013","8-2013","9-2013","10-2013","11-2013","12-2013","1-2014","2-2014","3-2014",
                "4-2014","5-2014","6-2014")
month_summary <- within(month_summary, interview_month <- factor(interview_month, levels=month_order))
month_summary$haplotype <- as.integer(as.character(month_summary$haplotype))

case_children = read.csv('ama_case_children.csv')
case_children
case_children = merge(case_children)
x=table(ama_merged$interview_month)

# create a plot of the presence and abundance of each haplotype over time
ama_month_plot = ggplot(month_summary, aes(x=interview_month, y=haplotype, size=n, color=total_reads)) +
  geom_point() +
  labs(title = "Pfama haplotype presence and abundance over time across MESA samples (1/1/18 sequencing run)",
       y = "Month and year",x="Haplotype", color = "Total number reads", size = "Number of participants")+
  scale_y_reverse(breaks = c(1:180))+
  theme_bw()+
  #scale_x_continuous(breaks=c(2013.4, 2014.6, 0.1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ama_month_plot

ggsave(ama_month_plot, filename="ama_month_plot.png", device="png",
       height=30, width=5, units="in", dpi=500)

