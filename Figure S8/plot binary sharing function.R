cc_share_score = read.csv("csp share score all283 case_child comp time dist.csv")
cc_share_score_edited = data.frame()

cc_share_score_edited=cc_share_score[!cc_share_score$case_child == cc_share_score$comp_cc, ]

cc_share_score_2013 = cc_share_score_2014 = data.frame()
cc_share_score_2013 = cc_share_score_edited[cc_share_score_edited$cc_date <2014, ]
cc_share_score_2013 = cc_share_score_2013[cc_share_score_2013$comp_cc_date <2014  ,]

cc_share_score_2014 = cc_share_score_edited[cc_share_score_edited$cc_date >2014,]
cc_share_score_2014 = cc_share_score_2014[cc_share_score_2014$comp_cc_date >2014  ,]

require (ggplot2)
require (dplyr)


ggplot(cc_share_score_edited, aes(x=temporal_dist, y=share_score))+
  geom_point(size=0.5, alpha=1)+
  geom_smooth(method = "loess",  se=FALSE, aes(size=5))+
  xlab("Distance (years)")+
  ylab("Probability of sharing")+
  scale_x_continuous (limits = c(0,1.125),breaks = seq(0,1.25,0.25))   +
  scale_y_continuous (limits = c(0,1.0),breaks = seq(0,1.0,0.25))+
  theme_bw()

ggsave("csp_cc_binary_score_all.jpg", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 5, height = 5, 
       dpi = 300, limitsize = TRUE)
