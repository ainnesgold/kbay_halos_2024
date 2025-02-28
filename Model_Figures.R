##Temp and Nutrient sensitivity
source('Model_Sensitivity_SST.R')
#Top panel of Figure 3 (model methods figure)
figure3top<-ggarrange(sst_plot, sst_function_plot, q_gg_sst, nrow=1, ncol = 3, common.legend = TRUE, legend="right")

ggsave("~/Desktop/figure3top.png", figure3top, width = 14, height=5, bg="transparent")

#Figure 5
figure5<-ggarrange(a_gg, h_gg, a_gg2, h_gg2, nrow=2, ncol=2, common.legend = TRUE, legend="right")

ggsave("~/Desktop/figure5.png", figure5, width = 10, height=8, bg="transparent")




source('Model_Sensitivity_Nutrients.R')

#Bottom panel of Figure 3 (model methods figure)
figure3bottom<-ggarrange(nutrient_plot, nutrient_function_plot, q_gg_nutrients, nrow=1, ncol=3, common.legend = TRUE, legend="right")

ggsave("~/Desktop/figure3bottom.png", figure3bottom, width = 14, height=5, bg="transparent")

#Figure S3
figureS3<-ggarrange(a_gg, h_gg, a_gg2, h_gg2, nrow=2, ncol=2, common.legend = TRUE, legend="right")

ggsave("~/Desktop/figureS3.png", figureS3, width = 10, height=8, bg="transparent")