cc_L1norm = read.csv("csp L1L2norm all283 case_child comp time dist.csv")
cc_L1norm_edited = data.frame()

cc_L1norm_edited=cc_L1norm[!cc_L1norm$case_child == cc_L1norm$comp_cc, ]


cc_L1norm_2013 = cc_L1norm_2014 = data.frame()
cc_L1norm_2013 = cc_L1norm_edited[cc_L1norm_edited$cc_date <2014, ]
cc_L1norm_2013 = cc_L1norm_2013[cc_L1norm_2013$comp_cc_date <2014  ,]

cc_L1norm_2014 = cc_L1norm_edited[cc_L1norm_edited$cc_date >2014,]
cc_L1norm_2014 = cc_L1norm_2014[cc_L1norm_2014$comp_cc_date >2014  ,]

require (ggplot2)
require (dplyr)

ggplot(cc_L1norm_edited, aes(x=temporal_dist, y=L1norm))+
  geom_point(size=0.05, alpha=0.1)+
  geom_smooth(method = "loess",  se=FALSE, aes(size=5))+
  xlab("Distance (years)")+
  ylab("Probability of sharing")+
  scale_x_continuous (limits = c(0,1.125),breaks = seq(0,1.125,0.25))   +
  scale_y_continuous (limits = c(0,20),breaks = seq(0,20,5))+
  theme_bw()


ggsave("csp_cc_L1norm_geog_all.jpg", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 5, height = 5, 
       dpi = 300, limitsize = TRUE)
