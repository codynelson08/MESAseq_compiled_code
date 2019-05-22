#Function to plot maps

library (ggmap)
library(ggplot2)
library (tidyverse)
library(lubridate)

outbreak = read.csv ("may_july_cc.csv")
outbreak=outbreak[1:77,]

#register_google(key = "AIzaSyB0cSXIZlBtn2xyvM_x5pvI6ve6vlr_A_Q", write = TRUE)

outbreak$interview_date=mdy(outbreak$interview_date)
outbreak$Outbreak = as.character(outbreak$Outbreak)

outbreak$interview_date=decimal_date(outbreak$interview_date)
map = get_stamenmap( bbox = c(left = 34.66, bottom = 0.47, right = 34.82, top = 0.75),
        zoom = 11, scale = "auto", maptype = "toner-background", source = "stamen",
         messaging = FALSE,
        urlonly = FALSE, filename = NULL, crop = TRUE,  language = "en-EN", force=TRUE)
ggmap(map) + 
  geom_jitter  (data = outbreak, height=0.005, width=0.005, aes(x = longitude, y = latitude, color=Outbreak))+
  scale_color_manual(values = c('0' = "#999998", '1' = "#E90006"))+
  theme_bw()

ggsave("cc_outbreak.jpg", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 5, height = 5, 
       dpi = 300, limitsize = TRUE)

      