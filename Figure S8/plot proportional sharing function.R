cc_weighted_share_score = read.csv("csp weighted share score all283 case_child comp time dist.csv")
cc_weighted_share_score_edited = data.frame()

cc_weighted_share_score_edited=cc_weighted_share_score[!cc_weighted_share_score$case_child == cc_weighted_share_score$comp_cc, ]
cc_weighted_share_score_edited$share_score = cc_weighted_share_score_edited$share_score / 100

cc_weighted_share_score_2013 = cc_weighted_share_score_2014 = data.frame()
cc_weighted_share_score_2013 = cc_weighted_share_score_edited[cc_weighted_share_score_edited$cc_date <2014, ]
cc_weighted_share_score_2013 = cc_weighted_share_score_2013[cc_weighted_share_score_2013$comp_cc_date <2014  ,]

cc_weighted_share_score_2014 = cc_weighted_share_score_edited[cc_weighted_share_score_edited$cc_date >2014,]
cc_weighted_share_score_2014 = cc_weighted_share_score_2014[cc_weighted_share_score_2014$comp_cc_date >2014  ,]

require (ggplot2)
require (dplyr)

ggplot(cc_weighted_share_score_edited, aes(x=as.numeric(temporal_dist), y=as.numeric(weighted_share_score)))+
  geom_point(size=0.5, alpha=1)+
  geom_smooth(method = "loess",  se=FALSE, aes(size=5))+
  xlab("Distance (years)")+
  ylab("Probability of sharing")+
  scale_x_continuous (limits = c(0, 1.125),breaks = seq(0,1.125,0.25))   +
  scale_y_continuous (limits = c(0,100),breaks = seq(0,100,25))+
  theme_bw()

ggsave("csp_cc_weighted_share_score_all.jpg", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 5, height = 5, 
       dpi = 300, limitsize = TRUE)
