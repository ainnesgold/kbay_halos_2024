#Reading in environmental data

########################################### TEMPERATURE #####################################################
temp_raw <- read.csv("Data/temp_combined2.csv")
temp_raw$Date <- mdy(temp_raw$Date)

# Extract day, month, and year into separate columns
temp_raw$Day <- day(temp_raw$Date)
temp_raw$Month <- month(temp_raw$Date)
temp_raw$Year <- year(temp_raw$Date)

temp_raw$avg_temp <- rowMeans(temp_raw[, c("Temp1", "Temp2", "Temp3")], na.rm = TRUE)
temp_raw$avg_light <- rowMeans(temp_raw[, c("Intensity1", "Intensity2", "Intensity3")], na.rm = TRUE)

#create a mean temp for each day
temp_mean <- temp_raw %>%
  group_by(Date, Day, Month, Year) %>%
  summarise(daily_temp = mean(avg_temp, na.rm=TRUE))

# Create a rolling average and rolling max column using RcppRoll
#temp_mean$rolling_average_30 <- roll_mean(temp_mean$daily_temp, n = 30, align = "right", fill = NA)
temp_mean$rolling_average_14 <- roll_mean(temp_mean$daily_temp, n = 14, align = "right", fill = NA)
temp_mean$rolling_average_7 <- roll_mean(temp_mean$daily_temp, n = 7, align = "right", fill = NA)

#temp_mean$rolling_max_30 <- roll_max(temp_mean$daily_temp, n = 30, align = "right", fill = NA)
temp_mean$rolling_max_14 <- roll_max(temp_mean$daily_temp, n = 14, align = "right", fill = NA)
temp_mean$rolling_max_7 <- roll_max(temp_mean$daily_temp, n = 7, align = "right", fill = NA)




########################################### NUTRIENTS #####################################################

nutrient_raw <- read.csv("Data/nutrients_combined.csv")
colnames(nutrient_raw)[colnames(nutrient_raw) == "Reef_Number"] <- "Reef"
nutrient_raw$Reef <- paste0("R", nutrient_raw$Reef)

nutrient_sum <- nutrient_raw %>%
  group_by(Month, Year, Reef) %>%
  summarize(avg_total_N = mean(Total_N),
            avg_total_P = mean(Total_P),
            avg_phosphate = mean(Phosphate),
            avg_sillicate = mean(Sillicate),
            avg_NN = mean(N.N),
            avg_ammonia = mean(Ammonia),
            avg_chl = mean(Chlorophyll))





########################################### FISH #####################################################
fish <- read.csv("Data/fish surveys.csv")

#l-w data
lw <- read.csv("Data/length-weight.csv")
lw <- lw[,-8] #delete notes column
fishlw <- merge(fish, lw, by = c("Scientific_Name", "Common_Name", "Family"))
fishlw$mass <- fishlw$Intercept * fishlw$Size_cm ^ fishlw$Slope

fishlw_sum <- fishlw %>%
  group_by(Month, Year, Site) %>%
  summarize(total_number = sum(Number), total_mass = sum(mass))

fishlw_sum_herbivore <- fishlw %>%
  filter(Herbivore == "Y") %>%
  group_by(Month, Year, Site) %>%
  summarize(total_number_herbivore = sum(Number), total_mass_herbivore = sum(mass))

fishlw_sum <- merge(fishlw_sum, fishlw_sum_herbivore, by = c("Month", "Year", "Site"))

colnames(fishlw_sum)[3] ="Reef"