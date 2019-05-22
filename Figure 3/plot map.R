#Function to plot maps for csp (Figure 3A-C)

library (ggmap)
library(ggplot2)
library (tidyverse)

allhaps = read.csv ("csp_allhaps_demog_person_age_decimaldate_latlong.csv")

#remove all haplotypes in compiled data frame except for haplotype of interest
allhaps = allhaps[which(allhaps$haplotype =="TTATTAATCCTATTGAACTATTTACGACATTAAACACACTGGAACATTTTTCCATTTTACAAATTTTTTTTTCAATATCATTTTCATAATCTAATTGGTCTTTAGGTTTATTAGCAGAGCCAGGCTTTATTCTAACTTGAATACCATTTCCACAAGTTACACTACATGGGGACCATTCAGTTGAAATAGAATTTTTTATTTTCGTTAAATATTCTTTTATGTGCTTATCACTTGGTTCTTCGTTATTATTATTTTTTACAGCATTGTTGGCATTAGCATTTTCATCTA"),]

#create map using stamen database
map = get_stamenmap( bbox = c(left = 34.66, bottom = 0.47, right = 34.82, top = 0.75),
        zoom = 11, scale = "auto", maptype = "toner-background", source = "stamen",
         messaging = FALSE,
        urlonly = FALSE, filename = NULL, crop = TRUE,  language = "en-EN", force=TRUE)

#plot map
ggmap(map) + 
  geom_jitter  (data = allhaps, height=0.005, width=0.005, aes(x = longitude, y = latitude, color=date))+
  scale_color_gradient(low="blue", high="red", limits=c(2013.2, 2014.7))+
  theme_bw()

ggsave("cc_share_score_2013_csp.jpg", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 5, height = 10, 
       dpi = 300, limitsize = TRUE)

      