# kbay_halos_2024

#This repo contains all the data and code needed to reproduce figures and analyses for the manuscript titled "Herbivory and temperature mediate coral reef halo dynamics". 
#This study explores how herbivore biomass, temperature, and nutrients contribute to fluctuations in halo presence and size. We conducted a field study using artificial reefs, and used a consumer-resource model to simulate these effects as well.

#Data:
#All collected field data are in csv files in the data folder. 
#A description of each CSV file and the variables used:
#2022_noaa.csv is the NOAA buoy temperature data used for the two months of data lost by our loggers. 
  #It contains columns year (YY), month (MM), day (DD), hour (hh), minute (mm), and water temperature (WTMP). 
  #It is used in field_data_analysis.R.
#fish surveys.csv contains our field data on fish biomass. 
  #It contains columns Date, Month, Year, Time, Site, Scientific_Name, Common_Name, Family, Size_cm (referring to fish length in cm), Number (referring to abundance of species), and Notes (a column of field notes). 
  #It is used in reading_env_data.R, which is sourced in field_data_analysis.R.
#There are three CSV files that contain halo measurements, each performed by a different coauthor of the paper. 
  #These are halo_data_ann.csv, halo_data_dava.csv, and halo_data_josh.csv. 
  #They each contain the same columns which are: ImageName (referring to the file name), Date, Month, Year, Day, Reef (referring to the artificial reef number), Transect (referring to the transect number, 1-3), DistanceTransect (referring to the distance from the central structure in cm), BinNumber (referring distance bin number, starting with 1 closest to the structure), NumberAlgae (referring to the number of random points landing on algae), NumberSand (referring to the number of random points landing on sand ), QuadratSize (referring to the size of the quadrat, always 20x20 cm), Notes, and DominantVegetation (referring to the dominant vegetation type; only done for some photos and not used in this analysis).
  #They are used in field_data_analysis.R.
#length-weight.csv contains the length-weight relationships we used to convert abundance and length to biomass. 
  #It contains columns: Scientific name, Common name, Herbivore (referring to herbivore classification; yes or no), Intercept (referring to the length-weight equation intercept), Slope (referring to the length-weight equation slope), Source, and Notes.
  #It is used in reading_env_data.R, which is sourced in field_data_analysis.R.
#nutrients_combined.csv is our field data on water nutrient levels. 
  #It contains the columns: Month, Year, Reef_Number (referring to the artificial reef number, 1-3), Total_N (referring to total nitrogen), Total_P (referring to total phosphorous), Phosphate, Silicate, N+N (referring to Nitrate + Nirite), Ammonia, and Chlorophyll. 
  #All nutrients are in umol/L units, and chlorophyll is in ug/L.
  #It is used in reading_env_data.R,  which is sourced in field_data_analysis.R.
#Temp_combined2.csv is our field data on water temperature. 
  #It contains columns Date, Temp1 (referring to water temperature at sensor 1), Intensity1 (referring to light intensity at sensor 1), Temp2 (referring to water temperature at sensor 2), Intensity2 (referring to light intensity at sensor 2), Temp3 (referring to water temperature at sensor 3), and Intensity3 (referring to light intensity at sensor 3). 
  #Water temperature is in degrees F, intensity measured in lux. Both measured by hobo loggers. 
  #It is used in reading_env_data.R,  which is sourced in field_data_analysis.R.
#Temp_combined.csv is the same field data as above, but with additional Month and Year columns. 
  #It is used as an input to the consumer-resource modeling portion of the study (see Model_Sensitivity_SST.R)



#Field data analysis:

#field_data_analysis.R runs all the statistical models and creates field data figures (Figures 4, S1, S2, all supplemental tables of model results). This script calls reading_field_data.R, which reads in temperature, nutrient, and fish data.

#Consumer-resource model analysis:

#Model_Figures.R creates the modeling figures (Figure 3, 5, S3). This script calls Model_Sensitivity_SST.R, which runs the temperature dependent version of the mechanisic model, and Model_Sensitivty_Nutrients.R, which runs the nutrient dependent version of the mechanistic model. Both those scripts call Stable_Cycle_Baseline.R, which is the baseline verison of the model from Ong et al. (2025). 

#Simulated_Sampling.R creates the simulated sampling figures (Figure 6, S4)

#Simulated_Sampling.alt_parameters.R creates the simulate sampling alternative parameter figure  (Figure S5)

#reef_distances.R has the calculations for the two versions of Clark Evans R used in this study

#density_dependence_exercise.R creates Figure S3


#R and package versons:
#This code was generated using R Version 2025.05.0+496.
#Package versions used: deSolve_1.40, kableExtra_1.4.0, broom.mixed_0.2.9.5, reshape2_1.4.4, lmerTest_3.1-3, lme4_1.1-35.3, Matrix_1.7-0, RcppRoll_0.3.0, ggpubr_0.6.0, lubridate_1.9.3, forcats_1.0.0, stringr_1.5.1, dplyr_1.1.4, purrr_1.0.2, readr_2.1.5, tidyr_1.3.1, tibble_3.2.1, ggplot2_3.5.2, tidyverse_2.0.0, mgcv_1.9-1, nlme_3.1-164       
