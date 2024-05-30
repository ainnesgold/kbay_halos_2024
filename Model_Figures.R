##Temp and Nutrient sensitivity
source('Model_Sensitivity_SST.R')
#Top panel of Figure 3 (model methods figure)
ggarrange(sst_plot, sst_function_plot, q_gg_sst, nrow=1, ncol = 3, common.legend = TRUE, legend="right") #saved 14 x 5

#Figure 5
ggarrange(a_gg, h_gg, a_gg2, h_gg2, nrow=2, ncol=2, common.legend = TRUE, legend="right") #saved 12 x 8




source('Model_Sensitivity_Nutrients.R')

#Bottom panel of Figure 3 (model methods figure)
ggarrange(nutrient_plot, nutrient_function_plot, q_gg_nutrients, nrow=1, ncol=3, common.legend = TRUE, legend="right") #save 14 x 5

#Figure S3
ggarrange(a_gg, h_gg, a_gg2, h_gg2, nrow=2, ncol=2, common.legend = TRUE, legend="right")
