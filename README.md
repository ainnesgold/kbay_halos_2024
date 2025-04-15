# kbay_halos_2024

#This repo contains all the data and code needed to reproduce figures and analyses for the manuscript titled "Herbivory and temperature mediate coral reef halo dynamics". This study explores how herbivore biomass, temperature, and nutrients contribute to fluctuations in halo presence and size. We conducted a field study using artificial reefs, and used a consumer-resource model to simulate these effects as well.

#Data:
#All collected field data are in csv files in the data folder. 2022_noaa surveys is the NOAA buoy temperature data used for the two months of data lost by our loggers, fish surveys.csv is our fish survey data, three halo_data.csv files contain the halo measurements, length-weight.csv contains the length-weight relationships we used to convert abundance and length to biomass, nutrients_combined.csv is our collected nutrient data, and Temp_combined2.csv is our collected temperature data.


#Field data analysis:

#field_data_analysis.R runs all the statistical models and creates field data figures (Figures 4, S1, S2, all supplemental tables of model results). This script calls reading_field_data.R, which reads in temperature, nutrient, and fish data.


#Consumer-resource model analysis:

#Model_Figures.R creates the modeling figures (Figure 3, 5, S3). This script calls Model_Sensitivity_SST.R, which runs the temperature dependent version of the mechanisic model, and Model_Sensitivty_Nutrients.R, which runs the nutrient dependent version of the mechanistic model. Both those scripts call Stable_Cycle_Baseline.R, which is the baseline verison of the model from Ong et al. (2025). 

#Simulated_Sampling.R creates the simulated sampling figures (Figure 6, S4)

#Simulated_Sampling.alt_parameters.R creates the simulate sampling alternative parameter figure  (Figure S5)

#reef_distances.R has the calculations for the two versions of Clark Evans R used in this study

#density_dependence_exercise.R creates Figure S3