
##Field data analysis

#This script contains all the statistical models and creates field data figures (Figures 4, S1, S2)


#load packages
library(mgcv)
library(tidyverse)
library(ggpubr)
library(lubridate)
library(RcppRoll)
library(lmerTest)
library(reshape2)
library(broom.mixed)
library(kableExtra)


#Read in Vegetation Data
#data sheet 1
raw <- read.csv("Data/halo_data_dava.csv", na.strings=c("","NA"))

#data sheet 2
raw2 <- read.csv("Data/halo_data_josh.csv", na.strings=c("","NA"))

#data sheet 3
raw3 <- read.csv("Data/halo_data_ann.csv", na.strings=c("","NA"))

#Combine them
raw <- rbind(raw, raw2, raw3)
raw$Date <- as.Date(paste(raw$Year, raw$Month, raw$Day, sep = "-"), format = "%Y-%m-%d")




########################################### VEGETATION COVER AND HALO SIZE #####################################################
#Get total points and vegetation cover
raw$total_points <- raw$NumberAlgae + raw$NumberSand
raw$algae_cover <- raw$NumberAlgae / raw$total_points


#First get mean vegetation cover for each bin - avg of transects
data_sum_bin <- raw %>% 
  group_by(Year, Month, Day, Date, Reef, BinNumber) %>%
  summarise(veg_cover = mean(algae_cover, na.rm = TRUE))


#Then mean vegetation cover for each reef - overall reef avg
data_sum_reef <- raw %>%
  group_by(Year, Month, Day, Reef) %>% #
  summarize(mean_cover = mean(algae_cover, na.rm = TRUE),
            sd_veg = sd(algae_cover, na.rm = TRUE))


mean_veg <- merge(data_sum_bin, data_sum_reef, by = c("Year", "Month", "Day", "Reef"))
mean_veg <- mean_veg %>% arrange(Year, Month, Day, Reef)

mean_veg$NumericPart <- as.integer(sapply(strsplit(mean_veg$BinNumber, "B"), function(x) x[2]))
mean_veg <- mean_veg[order(mean_veg$NumericPart), ]
mean_veg <- mean_veg[order(mean_veg$Month, mean_veg$Year, mean_veg$Reef),]

# Calculate the threshold
mean_veg$threshold <- mean_veg$mean_cover - sd(mean_veg$mean_cover)
mean_veg$threshold2 <- mean_veg$mean_cover - (2*sd(mean_veg$mean_cover))

# Create the new column IsHalo
#two versions based on 1 or 2 SDs different from the mean
mean_veg$IsHalo <- ifelse(mean_veg$veg_cover > mean_veg$threshold, FALSE, TRUE)
mean_veg$IsHalo2 <- ifelse(mean_veg$veg_cover > mean_veg$threshold2, FALSE, TRUE)

# Loop through the data until a FALSE is encountered
count_true_until_false <- function(data) {
  count <- 0
  for (value in data) {
    if (value) {
      count <- count + 1
    } else {
      break
    }
  }
  return(count)
}

#This gives us halo size - 2 versions, one using 1SD as the cutoff, one using 2SD as the cutoff
result <- mean_veg %>%
  group_by(Month, Year, Day, Reef, mean_cover, sd_veg) %>% 
  summarize(HaloBins1SD = count_true_until_false(IsHalo),
            HaloBins2SD = count_true_until_false(IsHalo2))

result$halosize <- result$HaloBins1SD * 20
result$halosize2 <- result$HaloBins2SD * 20

hist(result$halosize, breaks = seq(min(result$halosize), max(result$halosize) + 10, by = 10))

#write.csv(mean_veg, "rawdata.csv")

############################# TEMPERATURE #################################################
source('reading_env_data.R')
temp_mean$Date <- as.Date(paste(temp_mean$Year, temp_mean$Month, temp_mean$Day, sep = "-"), format = "%Y-%m-%d")
filtered_temp_mean <- temp_mean %>%
  filter(Date %in% raw$Date)
filtered_temp_mean$rolling_average_14 <- (5/9) * (filtered_temp_mean$rolling_average_14 - 32)
filtered_temp_mean$rolling_average_7 <- (5/9) * (filtered_temp_mean$rolling_average_7 - 32)
filtered_temp_mean$rolling_max_14 <- (5/9) * (filtered_temp_mean$rolling_max_14 - 32)
filtered_temp_mean$rolling_max_7 <- (5/9) * (filtered_temp_mean$rolling_max_7 - 32)
filtered_temp_mean$daily_temp <- (5/9) * (filtered_temp_mean$daily_temp - 32)


#noaa buoy data for missing months
noaa_2022 <- read.csv("Data/2022_noaa.csv")
names(noaa_2022)[names(noaa_2022) == 'YY'] <- 'Year'
names(noaa_2022)[names(noaa_2022) == 'MM'] <- 'Month'
names(noaa_2022)[names(noaa_2022) == 'DD'] <- 'Day'

noaa_2022 <- noaa_2022[order(noaa_2022$Year, noaa_2022$Month, noaa_2022$Day), ]
noaa_2022 <- noaa_2022 %>%
  group_by(Year, Month, Day) %>%
  summarise(mean_WTMP = mean(WTMP))

noaa_2022$rolling_average_14 <- roll_mean(noaa_2022$mean_WTMP, n = 14, align = "right", fill = NA)
noaa_2022$rolling_average_7 <- roll_mean(noaa_2022$mean_WTMP, n = 7, align = "right", fill = NA)

noaa_2022$rolling_max_14 <- roll_max(noaa_2022$mean_WTMP, n = 14, align = "right", fill = NA)
noaa_2022$rolling_max_7 <- roll_max(noaa_2022$mean_WTMP, n = 7, align = "right", fill = NA)


noaa_2022$Date <- as.Date(paste(noaa_2022$Year, noaa_2022$Month, noaa_2022$Day, sep = "-"), format = "%Y-%m-%d")

noaa_2022 <- noaa_2022 %>%
  filter(Month == 8 | Month == 9)

noaa_2022 <- noaa_2022 %>%
  filter(Date %in% raw$Date)

noaa_2022 <- noaa_2022 %>%
  rename(daily_temp = mean_WTMP)

full_temp <- merge(filtered_temp_mean, noaa_2022, by = c("Date", "Day", "Month", "Year", "daily_temp", "rolling_average_14",
                                                         "rolling_average_7", "rolling_max_14", "rolling_max_7"), all = TRUE)
full_temp <- full_temp %>%
  arrange(Date, daily_temp) %>%
  distinct(Date, .keep_all = TRUE)   

#merge data
full_data <- merge(result, full_temp, by = c("Month", "Year", "Day"), all = TRUE)

full_data <- full_data %>%
  arrange(Year, Month, Day)


############################# NUTRIENTS #################################################


full_data <- merge(full_data, nutrient_raw, by = c("Month", "Year", "Reef"), all = TRUE)



############################# FISH #################################################
full_data <- merge(full_data,  fishlw_sum,
                   by = c("Month", "Reef", "Year"), all = TRUE)

###Halo data ONLY - for analyses looking only at observations with halos
halo_veg <- merge(result, mean_veg)

halo_veg <- halo_veg %>%
  filter(HaloBins1SD > 0 & IsHalo == TRUE) %>%
  subset(NumericPart <= HaloBins1SD) %>%
  group_by(Year, Month, Day, Reef) %>%
  summarize(halo_veg_cover = mean(veg_cover))

full_data_halos_only <- merge(full_data, halo_veg, by = c("Year", "Month", "Day", "Reef"), all = FALSE)

#numeric reef random effect for GAM
full_data$Reef2 <- ifelse(full_data$Reef == "R1", 1,
                          ifelse(full_data$Reef == "R2", 2,
                                 ifelse(full_data$Reef == "R3", 3, full_data$Reef)))
full_data$Reef2 <- as.numeric(full_data$Reef2)

#same for halo only data
full_data_halos_only$Reef2 <- ifelse(full_data_halos_only$Reef == "R1", 1,
                                     ifelse(full_data_halos_only$Reef == "R2", 2,
                                            ifelse(full_data_halos_only$Reef == "R3", 3, full_data_halos_only$Reef)))
full_data_halos_only$Reef2 <- as.numeric(full_data_halos_only$Reef2)

####### Looking at non herbivores
full_data$total_mass_nonherbivore <- full_data$total_mass - full_data$total_mass_herbivore
full_data_halos_only$total_mass_nonherbivore <- full_data_halos_only$total_mass - full_data_halos_only$total_mass_herbivore



##General trends over time
p1<-ggplot(full_data, aes(x = Date, y = halosize, color = Reef)) + 
  geom_line(linewidth = 1) +
  labs(x = NULL, y = "Halo width (cm)") +
  theme_minimal() +
  theme(
    text = element_text(size = 20),  # Set text size to 20
    axis.title = element_text(size = 20),  # Set axis title size to 20
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)# Set axis text size to 20
  )

p2<-ggplot(full_data, aes(x = Date, y = mean_cover, color = Reef)) + 
  geom_line(linewidth = 1) +
  labs(x = NULL, y = "Mean veg cover (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 20),  # Set text size to 20
    axis.title = element_text(size = 20),  # Set axis title size to 20
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)# Set axis text size to 20
  )

p3<-ggplot(full_data_halos_only, aes(x = Date, y = halo_veg_cover, color = Reef)) + 
  geom_point(size = 3) +
  labs(x = NULL, y = "Halo veg cover (%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 20),  # Set text size to 20
    axis.title = element_text(size = 20),  # Set axis title size to 20
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)# Set axis text size to 20
  )

p4<-ggplot(full_data, aes(x = Date, y = rolling_average_14)) + 
  geom_line(linewidth = 1) +
  labs(x = NULL, y = "SST (°C)") +
  theme_minimal() +
  theme(
    text = element_text(size = 20),  # Set text size to 20
    axis.title = element_text(size = 20),  # Set axis title size to 20
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)# Set axis text size to 20
  )

p5<-ggplot(full_data, aes(x = Date, y = total_mass_herbivore, color = Reef)) + 
  geom_line(linewidth = 1) +
  labs(x = NULL, y = expression('Herbivore biomass (g/m'^2*')')) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),  # Set text size to 20
    axis.title = element_text(size = 20),  # Set axis title size to 20
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)# Set axis text size to 20
  )

p6<-ggplot(full_data, aes(x = Date, y = total_mass_nonherbivore, color = Reef)) + 
  geom_line(linewidth = 1) +
  labs(x = NULL, y = expression('Non herbivore biomass (g/m'^2*')')) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),  # Set text size to 20
    axis.title = element_text(size = 20),  # Set axis title size to 20
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)# Set axis text size to 20
  )

p7<-ggplot(full_data, aes(x = Date, y = N.N, color = Reef)) + 
  geom_line(linewidth = 1) +
  labs(x = NULL, y = "N+N (µmol/L)") +
  theme_minimal() +
  theme(
    text = element_text(size = 20),  # Set text size to 20
    axis.title = element_text(size = 20),  # Set axis title size to 20
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)# Set axis text size to 20
  )

p8<-ggplot(full_data, aes(x = Date, y = Phosphate, color = Reef)) + 
  geom_line(linewidth = 1) +
  labs(x = NULL, y = "Phosphate (µmol/L)") +
  theme_minimal() +
  theme(
    text = element_text(size = 20),  # Set text size to 20
    axis.title = element_text(size = 20),  # Set axis title size to 20
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)# Set axis text size to 20
  )

p9<-ggplot(full_data, aes(x = Date, y = Sillicate, color = Reef)) + 
  geom_line(linewidth = 1) +
  labs(x = NULL, y = "Silicate (µmol/L)") +
  theme_minimal() +
  theme(
    text = element_text(size = 20),  # Set text size to 20
    axis.title = element_text(size = 20),  # Set axis title size to 20
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)# Set axis text size to 20
  )

p10<-ggplot(full_data, aes(x = Date, y = Ammonia, color = Reef)) + 
  geom_line(linewidth = 1) +
  labs(x = NULL, y = "Ammonia (µmol/L)") +
  theme_minimal() +
  theme(
    text = element_text(size = 20),  # Set text size to 20
    axis.title = element_text(size = 20),  # Set axis title size to 20
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)# Set axis text size to 20
  )

p11<-ggplot(full_data, aes(x = Date, y = Chlorophyll, color = Reef)) + 
  geom_line(linewidth = 1) +
  labs(x = NULL, y = "Chl (µg/L)") +
  theme_minimal() +
  theme(
    text = element_text(size = 20),  # Set text size to 20
    axis.title = element_text(size = 20),  # Set axis title size to 20
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)# Set axis text size to 20
  )


############ FIGURE S2 ####################
figureS2 <- ggarrange(p1, p2, p3, p5, p6, p7, p8, p4, nrow=2, ncol=4, common.legend = TRUE)
 
ggsave("~/Desktop/FigureS2.png", plot = figureS2, width = 14, height = 8, bg = "transparent")

#can save as csvs if needed
#write.csv(full_data, "full_data.csv")
#write.csv(full_data_halos_only, "full_data_halos_only.csv")





############################# STATISTICAL MODELS #################################################

############################## HALO SIZE #############################################

#
#halo only data

m1_gam <- gam(halosize ~ mean_cover +
                s(rolling_average_14, k = 6) +
                total_mass_herbivore +
                total_mass_nonherbivore +
                s(rolling_average_14, by = total_mass_herbivore) +
                s(Month, Reef2, bs="re"),
              data = full_data_halos_only)

summary(m1_gam)

# Extract fixed effects
fixed_terms <- tidy(m1_gam, parametric = TRUE)

# Extract smooth terms
smooth_terms <- tidy(m1_gam, parametric = FALSE)

# Format the fixed effects table
fixed_table <- fixed_terms %>%
  kbl(digits = 3, col.names = c("Term", "Estimate", "Std. Error", "t-value", "p-value")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F) %>%
  add_header_above(c("Halo Width - Fixed Effects" = 5))

# Format the smooth terms table
smooth_table <- smooth_terms %>%
  kbl(digits = 3, col.names = c("Smooth Term", "Estimated df", "Reference df", "F Statistic", "p-value")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F) %>%
  add_header_above(c("Halo Width - Smooth Terms" = 5))

# Display the tables
fixed_table
smooth_table

# Plot of significant predictor
p1<-ggplot(full_data_halos_only, aes(x = rolling_average_14, y = halosize)) + 
  geom_point(position="jitter") +
  geom_smooth(method = "gam") +
  labs(x="SST (°C)", y = "Halo width (cm)") +
  ggtitle("C.")+
  theme_minimal() +
  theme(
    text = element_text(size = 20), 
    axis.title = element_text(size = 20),  
    axis.text = element_text(size = 20)  
  )







################################# HALO VEG DENSITY #############################################

# GAM
m1_gam <- gam(cbind(I(round(halo_veg_cover*100)), 100)  ~ s(mean_cover, k=3) +
                rolling_average_14 +
                total_mass_herbivore +
                total_mass_nonherbivore +
                rolling_average_14*total_mass_herbivore +
                s(Month, Reef2, bs="re"),
              family = binomial,
              data = full_data_halos_only)


summary(m1_gam)

# Extract fixed effects
fixed_terms <- tidy(m1_gam, parametric = TRUE)

# Extract smooth terms
smooth_terms <- tidy(m1_gam, parametric = FALSE)

# Format the fixed effects table
fixed_table <- fixed_terms %>%
  kbl(digits = 3, col.names = c("Term", "Estimate", "Std. Error", "t-value", "p-value")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F) %>%
  add_header_above(c("Halo Veg Density - Fixed Effects" = 5))

# Format the smooth terms table
smooth_table <- smooth_terms %>%
  kbl(digits = 3, col.names = c("Smooth Term", "Estimated DF", "Reference df", "F Statistic", "p-value")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F) %>%
  add_header_above(c("Halo Veg Density - Smooth Terms" = 5))

######## Supplemental Tables 
fixed_table
smooth_table


#plots of significant predictors
p2<-ggplot(full_data_halos_only, aes(x = mean_cover, y = halo_veg_cover)) + 
  geom_point() +
  geom_smooth(method = "gam") +
  labs(x="Mean vegetation cover (%)", y = "Halo vegetation cover (%)") +
  ggtitle("D.")+
  theme_minimal() +
  theme(
    text = element_text(size = 20), 
    axis.title = element_text(size = 20),  
    axis.text = element_text(size = 20)  
  )

#Making a heatmap for the interaction
rolling_range <- seq(min(full_data_halos_only$rolling_average_14), max(full_data_halos_only$rolling_average_14), length.out = 100)
mass_range <- seq(min(full_data_halos_only$total_mass_herbivore), max(full_data_halos_only$total_mass_herbivore), length.out = 100)

# Create a grid of values for the interaction plot
interaction_grid <- expand.grid(rolling_average_14 = rolling_range,
                                total_mass_herbivore = mass_range,
                                total_mass_nonherbivore = mean(full_data_halos_only$total_mass_nonherbivore), # Set to the mean of nonherbivores
                                mean_cover = mean(full_data_halos_only$mean_cover),  # Set to the mean of mean_cover
                                Reef2 = unique(full_data_halos_only$Reef2),  # Unique levels of Reef
                                Month = unique(full_data_halos_only$Month))  # Unique levels of Month

# Predict the outcome (predicted probabilities) for each combination of values
interaction_grid$predicted_prob <- predict(m1_gam, newdata = interaction_grid, type = "response")

# Plot the interaction
p3<-ggplot(interaction_grid, aes(x = total_mass_herbivore, y = rolling_average_14)) +
  geom_tile(aes(fill = predicted_prob)) +
  scale_fill_viridis_c(name = "Halo vegetation cover (%)", direction=-1, 
                       breaks = c(0.0, 0.1, 0.2),
                       labels = c(0, 0.1, 0.2)) + 
  labs(x = expression('Herbivorous fish biomass (g/m'^2*')'), y = "SST (°C)") +
  ggtitle("E.")+
  theme_minimal() +
  theme(legend.position="top",
        text = element_text(size = 20),
        axis.title = element_text(size = 20), 
        axis.text = element_text(size = 20)  
  )






############################## MEAN VEGETATION #############################################

# one massive outlier for phosphate, filtered it out here
filtered_full_data <- full_data %>%
  filter(Phosphate < 0.2)

# mixed effects model
m1_lme <- glmer(
  cbind(I(round(mean_cover * 100)), 100) ~
    rolling_average_14 +
    total_mass_herbivore +
    N.N + 
    Phosphate +
    rolling_average_14 * total_mass_herbivore +
    (1 | Reef) + (1 | Month), 
  family = binomial,
  data = filtered_full_data
)
summary(m1_lme)

# table of results
model_summary <- tidy(m1_lme)

########## Supplemental Table
model_summary %>%
  kbl(digits = 3, col.names = c("Effect", "Group", "Term", "Estimate", "Std. Error", "Z value", "p-value")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F) %>%
  add_header_above(c("Mean Vegetation Model Output" = 7)) 


#checking residuals
plot(m1_lme, which = 1) 
residuals <- residuals(m1_lme)
# Q-Q plot
qqnorm(residuals)
qqline(residuals)

## plots of significant predictor
p4<-ggplot(filtered_full_data, aes(x = rolling_average_14, y = mean_cover)) + 
  geom_point() +
  geom_smooth(method = "glm") +
  labs(x = "SST (°C)", y = "Mean vegetation cover (%)") +
  ggtitle("A.")+
  theme_minimal() +
  theme(
    text = element_text(size = 20),  # Set text size to 20
    axis.title = element_text(size = 20),  # Set axis title size to 20
    axis.text = element_text(size = 20)  # Set axis text size to 20
  )




################### HALO PRESENCE ###################

#creating a binomial dataframe for halo presence
binomial_df <- full_data
binomial_df$binom_halo1 <- ifelse(binomial_df$halosize > 0, 1, 0)
binomial_df$binom_halo2 <- ifelse(binomial_df$halosize2 > 0, 1, 0)

#using halosize1 (1 SD difference) as the response - not many observations qualified as halos if using the 2 SD cutoff

binom_m1 <- glmer(binom_halo1 ~ rolling_average_14 + total_mass_herbivore + total_mass_nonherbivore + mean_cover +
                    rolling_average_14 * total_mass_herbivore +
                    (1 | Reef) + (1 | Month), data = binomial_df, family = binomial)
summary(binom_m1)
hist(binomial_df$binom_halo1)

#table of results
model_summary <- tidy(binom_m1)

########## Supplemental table
model_summary %>%
  kbl(digits = 3, col.names = c("Effect", "Group", "Term", "Estimate", "Std. Error", "Z value", "p-value")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F) %>%
  add_header_above(c("Halo Presence Model Output" = 7))


#plots of significant predictors
ggplot(binomial_df, aes(x = total_mass_herbivore, y = binom_halo1)) +
  geom_point() +
  labs(y="Halo presence", x=expression('Herbivorous fish biomass (g/m'^2*')')) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),  
    axis.title = element_text(size = 20),  
    axis.text = element_text(size = 20)  
  )

#Interaction plot
rolling_range <- seq(min(binomial_df$rolling_average_14), max(binomial_df$rolling_average_14), length.out = 100)
mass_range <- seq(min(binomial_df$total_mass_herbivore), max(binomial_df$total_mass_herbivore), length.out = 100)

# Create a grid of values for the interaction plot
interaction_grid <- expand.grid(rolling_average_14 = rolling_range,
                                total_mass_herbivore = mass_range,
                                total_mass_nonherbivore = mean(binomial_df$total_mass_nonherbivore), # Set to the mean of nonherivore
                                mean_cover = mean(binomial_df$mean_cover),  # Set to the mean of mean_cover
                                Reef = unique(binomial_df$Reef),  # Unique levels of Reef
                                Month = unique(binomial_df$Month))  # Unique levels of Month

# Predict the outcome (predicted probabilities) for each combination of values
interaction_grid$predicted_prob <- predict(binom_m1, newdata = interaction_grid, type = "response")

# Plot the interaction
p7<-ggplot(interaction_grid, aes(x = total_mass_herbivore, y = rolling_average_14)) +
  geom_tile(aes(fill = predicted_prob)) +
  scale_fill_viridis_c(name = "Halo probability", option = "magma", 
                       breaks = c(0.25, 0.5, 0.75), 
                       labels = c(0.25, 0.5, 0.75)) + 
  labs(x = expression('Herbivorous fish biomass (g/m'^2*')'), y = "SST (°C)") +
  theme_minimal() +
  ggtitle("B.")+
  theme(legend.position="top",
        text = element_text(size = 20), 
        axis.title = element_text(size = 20), 
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 12)
  )



########## FIGURE 4
figure4_final <- ggarrange(p4, p7, p1, p2, p3, nrow=2, ncol=3) 

ggsave("~/Desktop/Figure4.png", plot = figure4_final, width = 12, height =8, bg = "transparent")


########### Making FIGURE S1

#For the mean veg model (hypothesis 1)
predictors <- filtered_full_data[, c("rolling_average_14", "N.N", "Phosphate",
                                     "total_mass_herbivore")]

colnames(predictors) <- c("SST", "Nitrate + nitrite", "Phosphate", "Herbivore biomass")


# Calculate pairwise correlation coefficients
correlation_matrix <- cor(predictors) 
print(correlation_matrix)

# Create heatmap of correlation matrix
p1<-ggplot(data = melt(correlation_matrix), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(name = "Correlation", low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  ggtitle("A.")+
  theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1), text = element_text(size = 20))

#hypotheses 2-4 (mean cover instead of nutrient predictors)
predictors <- filtered_full_data[, c("rolling_average_14", "mean_cover", "total_mass_herbivore", "total_mass_nonherbivore")]
colnames(predictors) <- c("SST", "Mean vegetation cover", "Herbivore biomass", "Non herbivore biomass")

# Calculate pairwise correlation coefficients
correlation_matrix <- cor(predictors) 
print(correlation_matrix)
# Create heatmap of correlation matrix
p2<-ggplot(data = melt(correlation_matrix), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(name = "Correlation", low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  ggtitle("B.")+
  theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1), text = element_text(size = 20))



###FIGURE S1
figureS1 <- ggarrange(p1, p2)

ggsave("~/Desktop/FigureS1.png", plot = figureS1, width = 12, height = 6, bg = "transparent")

