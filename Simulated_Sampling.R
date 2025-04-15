#Simulated Sampling - Figure 6
library(reshape2)
library(scales)

source('Model_Sensitivity_SST.R')

# Define initial conditions
initial_conditions <- c(A = 0.5, H = 0.434, q=0.8)
times <- seq(0, length(sst_values2), by=1)
q_values <- numeric(length(times))

#Define parameters - stable
parameters <- c(A0 = 0.8, R = 0.13, rc = 2, g = 2, s = 6, r = 8, k = 5, m = 0.03)


sensitivity_results <- perform_sensitivity_analysis(T_opt_values)

# Initialize an empty list to store dataframes for each T_opt
dataframes <- list()

# Create a dataframe for each T_opt and store it in the list
for (T_opt in names(sensitivity_results)) {
  out <- sensitivity_results[[T_opt]]
  df <- data.frame(Time = out[, 1], A = out[, 2], H = out[, 3], q = out[, 4])
  df$T_opt <- as.numeric(T_opt)
  dataframes[[T_opt]] <- df
}

# Combine all dataframes into one
combined_df <- do.call(rbind, dataframes)

# Organize combined dataframe by T_opt
combined_df <- combined_df[order(combined_df$T_opt, combined_df$Time), ]

q_values <- numeric(length(combined_df$Time))
q_values[1] = combined_df$q[1]
for (i in 2:length(combined_df$Time)) {
  new_value <- combined_df$q[i] - combined_df$q[i-1]
  q_values[i] <- new_value
}
combined_df$q <- q_values


######################## LOW OPTIMAL TEMP - 24 ######################################
# Define the limits
herbivore_limits <- c(0.94, 0.95)
seagrass_limits <- c(0.45, 0.55)
A0 = 0.8

# Ensure no NA values and filter data for T_opt == 24
filtered_data <- combined_df %>%
  filter(T_opt == 24) %>%
  filter(!is.na(A) & !is.na(H))

filtered_data$halosize <- A0 - filtered_data$A
filtered_data <- filtered_data %>%
  mutate(A_scaled = scales::rescale(halosize, to = herbivore_limits))

# Create a sequence of time points at every 30 units
time_points <- seq(200, max(combined_df$Time), by = 30)

# Calculate intersection points
intersection_points <- filtered_data %>%
  filter(Time %in% time_points) %>%
  select(Time, A_scaled, H) #%>%

# Combine the intersection points
intersection_data <- melt(intersection_points, id.vars = "Time", variable.name = "Line")

# Change the entries in column Line
intersection_data <- intersection_data %>%
  mutate(Line = ifelse(Line == "A_scaled", "Seagrass", ifelse(Line == "H", "Herbivore", Line)))

# Create the plot
p1<-ggplot(filtered_data %>% filter(Time >=200), aes(x = Time-200)) +
  geom_line(aes(y = H, color = "Herbivore"), size = 1.2) + 
  geom_line(aes(y = A_scaled, color = "Seagrass"), size = 1.2) + 
  scale_color_manual(values = c("Seagrass" = "#009E73", "Herbivore" = "#0072B2")) + 
  scale_y_continuous(
    name = "Herbivore density", 
    limits = herbivore_limits, 
    sec.axis = sec_axis(~ scales::rescale(., to = seagrass_limits), name = "Relative halo width", 
                        breaks = seq(seagrass_limits[1], seagrass_limits[2], by = 0.025), 
                        labels = number_format(accuracy = 0.01))
  ) +
  geom_vline(xintercept = time_points-200, linetype = "dashed", size = 0.5, color = "black") +
  geom_point(data = intersection_data, aes(x = Time-200, y = value, color = Line), size = 3) +
  labs(x = "Time (days)", y = "Density", color = "Legend") +
  ggtitle("A. Fixed; Optimal SST = 24°C") +
  theme_minimal() +
  theme(
    axis.title.y.right = element_text(color = "#009E73"),
    axis.text.y.right = element_text(color = "#009E73"),
    axis.title.y.left = element_text(color = "#0072B2"),
    axis.text.y.left = element_text(color = "#0072B2"),
    legend.position = "none", 
    text = element_text(size = 20)
  )


# Filter data for T_opt == 24 and Time in time_points
overlay_data <- combined_df %>%
  filter(T_opt == 24)

overlay_data <- merge(overlay_data, sst_df, by = "Time")

#filter only for where it intersects with survey points
overlay_data <- overlay_data %>%
  filter(Time %in% time_points)

p2<-ggplot(overlay_data, aes(x=H, y=SST, color = A0 - A)) + 
  geom_point(size=3) +
  labs(x = "Herbivore density", y = "SST (°C)", color = "Relative halo width") +
  ggtitle("") +
  scale_color_viridis_c(option="viridis") +
  theme_minimal() +
  theme(text = element_text(size = 20))



######################### alternate time point sampling #####################################
filtered_data <- combined_df %>%
  filter(T_opt == 24) %>%
  filter(!is.na(A) & !is.na(H))

filtered_data$halosize <- A0 - filtered_data$A
filtered_data <- filtered_data %>%
  mutate(A_scaled = scales::rescale(halosize, to = herbivore_limits))


# Create a sequence of time points at every 30 units
time_points2 <- seq(215, max(combined_df$Time), by = 30)

# Calculate intersection points
intersection_points <- filtered_data %>%
  filter(Time %in% time_points2) %>%
  select(Time, A_scaled, H) #%>%

# Combine the intersection points
intersection_data <- melt(intersection_points, id.vars = "Time", variable.name = "Line")

# Change the entries in column Line
intersection_data <- intersection_data %>%
  mutate(Line = ifelse(Line == "A_scaled", "Seagrass", ifelse(Line == "H", "Herbivore", Line)))

# Create the plot
p1_extra<-ggplot(filtered_data %>% filter(Time >=200), aes(x = Time-200)) +
  geom_line(aes(y = H, color = "Herbivore"), size = 1.2) + 
  geom_line(aes(y = A_scaled, color = "Seagrass"), size = 1.2) + 
  scale_color_manual(values = c("Seagrass" = "#009E73", "Herbivore" = "#0072B2")) + 
  scale_y_continuous(
    name = "Herbivore density", 
    limits = herbivore_limits, 
    sec.axis = sec_axis(~ scales::rescale(., to = seagrass_limits), name = "Relative halo width", 
                        breaks = seq(seagrass_limits[1], seagrass_limits[2], by = 0.025), 
                        labels = number_format(accuracy = 0.01))
  ) +
  geom_vline(xintercept = time_points2-200, linetype = "dashed", size = 0.5, color = "black") +
  geom_point(data = intersection_data, aes(x = Time-200, y = value, color = Line), size = 3) +
  labs(x = "Time (days)", y = "Density", color = "Legend") +
  ggtitle("A. Fixed, Optimal SST = 24°C") +
  theme_minimal() +
  theme(
    axis.title.y.right = element_text(color = "#009E73"),
    axis.text.y.right = element_text(color = "#009E73"),
    axis.title.y.left = element_text(color = "#0072B2"),
    axis.text.y.left = element_text(color = "#0072B2"),
    legend.position = "none", 
    text = element_text(size = 20)
  )


# Filter data for T_opt == 24 and Time in time_points
overlay_data <- combined_df %>%
  filter(T_opt == 24)

overlay_data <- merge(overlay_data, sst_df, by = "Time")

#filter only for where it intersects with survey points
overlay_data <- overlay_data %>%
  filter(Time %in% time_points2)

p2_extra<-ggplot(overlay_data, aes(x=H, y=SST, color = A0 - A)) + 
  geom_point(size=3) +
  labs(x = "Herbivore density", y = "SST (°C)", color = "Relative halo width") +
  ggtitle("") +
  scale_color_viridis_c(option="viridis") +
  theme_minimal() +
  theme(text = element_text(size = 20))




########################### HIGH 28 OPTIMAL TEMP 
filtered_data <- combined_df %>%
  filter(T_opt == 28) %>%
  filter(!is.na(A) & !is.na(H))

filtered_data$halosize <- A0 - filtered_data$A
filtered_data <- filtered_data %>%
  mutate(A_scaled = scales::rescale(halosize, to = herbivore_limits))

# Create a sequence of time points at every 30 units
time_points <- seq(200, max(combined_df$Time), by = 30)
#time_points <- seq(min(filtered_data$Time), max(filtered_data$Time), by = 30)

# Calculate intersection points
intersection_points <- filtered_data %>%
  filter(Time %in% time_points) %>%
  select(Time, A_scaled, H) #%>%

# Combine the intersection points
intersection_data <- melt(intersection_points, id.vars = "Time", variable.name = "Line")

# Change the entries in column Line
intersection_data <- intersection_data %>%
  mutate(Line = ifelse(Line == "A_scaled", "Seagrass", ifelse(Line == "H", "Herbivore", Line)))

# Create the plot
p3<-ggplot(filtered_data %>% filter(Time >=200), aes(x = Time-200)) +
  geom_line(aes(y = H, color = "Herbivore"), size = 1.2) + 
  geom_line(aes(y = A_scaled, color = "Seagrass"), size = 1.2) + 
  scale_color_manual(values = c("Seagrass" = "#009E73", "Herbivore" = "#0072B2")) + 
  scale_y_continuous(
    name = "Herbivore density", 
    limits = herbivore_limits, 
    sec.axis = sec_axis(~ scales::rescale(., to = seagrass_limits), name = "Relative halo width", 
                        breaks = seq(seagrass_limits[1], seagrass_limits[2], by = 0.025), 
                        labels = number_format(accuracy = 0.01))
  ) +
  geom_vline(xintercept = time_points-200, linetype = "dashed", size = 0.5, color = "black") +
  geom_point(data = intersection_data, aes(x = Time-200, y = value, color = Line), size = 3) +
  #geom_point(data = intersections, aes(x = Time, y = ifelse(Line == "Seagrass", scales::rescale(Value, to = herbivore_limits), Value), color = Line), size = 3) + 
  labs(x = "Time (days)", y = "Density", color = "Legend") +
  ggtitle("B. Fixed; Optimal SST = 28°C") +
  theme_minimal() +
  theme(
    axis.title.y.right = element_text(color = "#009E73"),
    axis.text.y.right = element_text(color = "#009E73"),
    axis.title.y.left = element_text(color = "#0072B2"),
    axis.text.y.left = element_text(color = "#0072B2"),
    legend.position = "none", 
    text = element_text(size = 20)
  )


# Filter data for T_opt == 24 and Time in time_points
overlay_data <- combined_df %>%
  filter(T_opt == 28)

overlay_data <- merge(overlay_data, sst_df, by = "Time")

#filter only for where it intersects with survey points
overlay_data <- overlay_data %>%
  filter(Time %in% time_points)

p4<-ggplot(overlay_data, aes(x=H, y=SST, color = A0 - A)) + 
  geom_point(size=3) +
  labs(x = "Herbivore density", y = "SST (°C)", color = "Relative halo width") +
  ggtitle("") +
  scale_color_viridis_c(option="viridis") +
  theme_minimal() +
  theme(text = element_text(size = 20))



######################### alternate time point sampling #####################################

filtered_data <- combined_df %>%
  filter(T_opt == 28) %>%
  filter(!is.na(A) & !is.na(H))

filtered_data$halosize <- A0 - filtered_data$A
filtered_data <- filtered_data %>%
  mutate(A_scaled = scales::rescale(halosize, to = herbivore_limits))

# Create a sequence of time points at every 30 units
#time points 2
time_points2 <- seq(215, max(combined_df$Time), by = 30)

# Calculate intersection points
intersection_points <- filtered_data %>%
  filter(Time %in% time_points2) %>%
  select(Time, A_scaled, H) #%>%

# Combine the intersection points
intersection_data <- melt(intersection_points, id.vars = "Time", variable.name = "Line")

# Change the entries in column Line
intersection_data <- intersection_data %>%
  mutate(Line = ifelse(Line == "A_scaled", "Seagrass", ifelse(Line == "H", "Herbivore", Line)))

# Create the plot
p3_extra<-ggplot(filtered_data %>% filter(Time >=200), aes(x = Time-200)) +
  geom_line(aes(y = H, color = "Herbivore"), size = 1.2) + 
  geom_line(aes(y = A_scaled, color = "Seagrass"), size = 1.2) + 
  scale_color_manual(values = c("Seagrass" = "#009E73", "Herbivore" = "#0072B2")) + 
  scale_y_continuous(
    name = "Herbivore density", 
    limits = herbivore_limits, 
    sec.axis = sec_axis(~ scales::rescale(., to = seagrass_limits), name = "Relative halo width", 
                        breaks = seq(seagrass_limits[1], seagrass_limits[2], by = 0.025), 
                        labels = number_format(accuracy = 0.01))
  ) +
  geom_vline(xintercept = time_points2-200, linetype = "dashed", size = 0.5, color = "black") +
  geom_point(data = intersection_data, aes(x = Time-200, y = value, color = Line), size = 3) +
  labs(x = "Time (days)", y = "Density", color = "Legend") +
  ggtitle("B. Fixed; Optimal SST = 28°C") +
  theme_minimal() +
  theme(
    axis.title.y.right = element_text(color = "#009E73"),
    axis.text.y.right = element_text(color = "#009E73"),
    axis.title.y.left = element_text(color = "#0072B2"),
    axis.text.y.left = element_text(color = "#0072B2"),
    legend.position = "none", 
    text = element_text(size = 20)
  )


# Filter data for T_opt == 24 and Time in time_points
overlay_data <- combined_df %>%
  filter(T_opt == 28)

overlay_data <- merge(overlay_data, sst_df, by = "Time")

#filter only for where it intersects with survey points
overlay_data <- overlay_data %>%
  filter(Time %in% time_points2)

p4_extra<-ggplot(overlay_data, aes(x=H, y=SST, color = A0 - A)) + 
  geom_point(size=3) +
  labs(x = "Herbivore density", y = "SST (°C)", color = "Relative halo width") +
  ggtitle("") +
  scale_color_viridis_c(option="viridis") +
  theme_minimal() +
  theme(text = element_text(size = 20))




##################################### CYCLES ######################################
# Define initial conditions
initial_conditions <- c(A = 0.5, H = 0.434, q=0.8)
q_values <- numeric(length(times))

#Define parameters
parameters <- c(A0 = 0.8, R = 1.28, rc = 2, g = 2, s = 6, r = 8, k = 5, m = 0.03)

sensitivity_results <- perform_sensitivity_analysis(T_opt_values)

# Initialize an empty list to store dataframes for each T_opt
dataframes <- list()

# Create a dataframe for each T_opt and store it in the list
for (T_opt in names(sensitivity_results)) {
  out <- sensitivity_results[[T_opt]]
  df <- data.frame(Time = out[, 1], A = out[, 2], H = out[, 3], q = out[, 4])
  df$T_opt <- as.numeric(T_opt)
  dataframes[[T_opt]] <- df
}

# Combine all dataframes into one
combined_df <- do.call(rbind, dataframes)

# Organize combined dataframe by T_opt
combined_df <- combined_df[order(combined_df$T_opt, combined_df$Time), ]

q_values <- numeric(length(combined_df$Time))
q_values[1] = combined_df$q[1]
for (i in 2:length(combined_df$Time)) {
  new_value <- combined_df$q[i] - combined_df$q[i-1]
  q_values[i] <- new_value
}
combined_df$q <- q_values




######################## PLOTS ######################################

######################## OPTIMAL TEMP 24 ######################################

# Assuming your combined_df has a column named 'Time'
# Define the limits
herbivore_limits2 <- c(0, 1)
seagrass_limits2 <- c(0, 0.5)
# Ensure no NA values and filter data for T_opt == 24
filtered_data <- combined_df %>%
  filter(T_opt == 24) %>%
  filter(!is.na(A) & !is.na(H))

filtered_data$halosize <- A0 - filtered_data$A
filtered_data <- filtered_data %>%
  mutate(A_scaled = scales::rescale(halosize, to = herbivore_limits2))

# Create a sequence of time points at every 30 units
time_points <- seq(200, max(combined_df$Time), by = 30)

# Calculate intersection points
intersection_points <- filtered_data %>%
  filter(Time %in% time_points) %>%
  select(Time, A_scaled, H) #%>%

# Combine the intersection points
intersection_data <- melt(intersection_points, id.vars = "Time", variable.name = "Line")

# Change the entries in column Line
intersection_data <- intersection_data %>%
  mutate(Line = ifelse(Line == "A_scaled", "Seagrass", ifelse(Line == "H", "Herbivore", Line)))

# Create the plot
p5<-ggplot(filtered_data %>% filter(Time >=200), aes(x = Time-200)) +
  geom_line(aes(y = H, color = "Herbivore"), size = 1.2) + 
  geom_line(aes(y = A_scaled, color = "Seagrass"), size = 1.2) + 
  scale_color_manual(values = c("Seagrass" = "#009E73", "Herbivore" = "#0072B2")) + 
  scale_y_continuous(
    name = "Herbivore density", 
    limits = herbivore_limits2, 
    sec.axis = sec_axis(~ scales::rescale(., to = seagrass_limits2), name = "Relative halo width", 
                        breaks = seq(seagrass_limits2[1], seagrass_limits2[2], by = 0.1), 
                        labels = number_format(accuracy = 0.1))
  ) +
  geom_vline(xintercept = time_points-200, linetype = "dashed", size = 0.5, color = "black") +
  geom_point(data = intersection_data, aes(x = Time-200, y = value, color = Line), size = 3) +
  labs(x = "Time (days)", y = "Density", color = "Legend") +
  ggtitle("C. Cyclic; Optimal SST = 24°C") +
  theme_minimal() +
  theme(
    axis.title.y.right = element_text(color = "#009E73"),
    axis.text.y.right = element_text(color = "#009E73"),
    axis.title.y.left = element_text(color = "#0072B2"),
    axis.text.y.left = element_text(color = "#0072B2"),
    legend.position = "none", 
    text = element_text(size = 20)
  )


# Filter data for T_opt == 24 and Time in time_points
overlay_data <- combined_df %>%
  filter(T_opt == 24)

overlay_data <- merge(overlay_data, sst_df, by = "Time")

#filter only for where it intersects with survey points
overlay_data <- overlay_data %>%
  filter(Time %in% time_points)

p6<-ggplot(overlay_data, aes(x=H, y=SST, color = A0 - A)) + 
  geom_point(size=3) +
  labs(x = "Herbivore density", y = "SST (°C)", color = "Relative halo width") +
  ggtitle("") +
  scale_color_viridis_c(option="viridis") +
  theme_minimal() +
  theme(text = element_text(size = 20))



######################### alternate time point sampling #####################################

filtered_data <- combined_df %>%
  filter(T_opt == 24) %>%
  filter(!is.na(A) & !is.na(H))

filtered_data$halosize <- A0 - filtered_data$A
filtered_data <- filtered_data %>%
  mutate(A_scaled = scales::rescale(halosize, to = herbivore_limits2))

# Create a sequence of time points at every 30 units
time_points2 <- seq(215, max(combined_df$Time), by = 30)

# Calculate intersection points
intersection_points <- filtered_data %>%
  filter(Time %in% time_points2) %>%
  select(Time, A_scaled, H) #%>%

# Combine the intersection points
intersection_data <- melt(intersection_points, id.vars = "Time", variable.name = "Line")

# Change the entries in column Line
intersection_data <- intersection_data %>%
  mutate(Line = ifelse(Line == "A_scaled", "Seagrass", ifelse(Line == "H", "Herbivore", Line)))

# Create the plot
p5_extra<-ggplot(filtered_data %>% filter(Time >=200), aes(x = Time-200)) +
  geom_line(aes(y = H, color = "Herbivore"), size = 1.2) + 
  geom_line(aes(y = A_scaled, color = "Seagrass"), size = 1.2) + 
  scale_color_manual(values = c("Seagrass" = "#009E73", "Herbivore" = "#0072B2")) + 
  scale_y_continuous(
    name = "Herbivore density", 
    limits = herbivore_limits2, 
    sec.axis = sec_axis(~ scales::rescale(., to = seagrass_limits2), name = "Relative halo width", 
                        breaks = seq(seagrass_limits2[1], seagrass_limits2[2], by = 0.1), 
                        labels = number_format(accuracy = 0.1))
  ) +
  geom_vline(xintercept = time_points2-200, linetype = "dashed", size = 0.5, color = "black") +
  geom_point(data = intersection_data, aes(x = Time-200, y = value, color = Line), size = 3) +
  labs(x = "Time (days)", y = "Density", color = "Legend") +
  ggtitle("C. Cyclic; Optimal SST = 24°C") +
  theme_minimal() +
  theme(
    axis.title.y.right = element_text(color = "#009E73"),
    axis.text.y.right = element_text(color = "#009E73"),
    axis.title.y.left = element_text(color = "#0072B2"),
    axis.text.y.left = element_text(color = "#0072B2"),
    legend.position = "none", 
    text = element_text(size = 20)
  )


# Filter data for T_opt == 24 and Time in time_points
overlay_data <- combined_df %>%
  filter(T_opt == 24)

overlay_data <- merge(overlay_data, sst_df, by = "Time")

#filter only for where it intersects with survey points
overlay_data <- overlay_data %>%
  filter(Time %in% time_points2)

p6_extra<-ggplot(overlay_data, aes(x=H, y=SST, color = A0 - A)) + 
  geom_point(size=3) +
  labs(x = "Herbivore density", y = "SST (°C)", color = "Relative halo width") +
  ggtitle("") +
  scale_color_viridis_c(option="viridis") +
  theme_minimal() +
  theme(text = element_text(size = 20))


########################### HIGH 28 OPTIMAL TEMP #
filtered_data <- combined_df %>%
  filter(T_opt == 28) %>%
  filter(!is.na(A) & !is.na(H))

filtered_data$halosize <- A0 - filtered_data$A
filtered_data <- filtered_data %>%
  mutate(A_scaled = scales::rescale(halosize, to = herbivore_limits2))

# Create a sequence of time points at every 30 units
time_points <- seq(200, max(combined_df$Time), by = 30)

# Calculate intersection points
intersection_points <- filtered_data %>%
  filter(Time %in% time_points) %>%
  select(Time, A_scaled, H) #%>%

# Combine the intersection points
intersection_data <- melt(intersection_points, id.vars = "Time", variable.name = "Line")

# Change the entries in column Line
intersection_data <- intersection_data %>%
  mutate(Line = ifelse(Line == "A_scaled", "Seagrass", ifelse(Line == "H", "Herbivore", Line)))

# Create the plot
p7<-ggplot(filtered_data %>% filter(Time >=200), aes(x = Time-200)) +
  geom_line(aes(y = H, color = "Herbivore"), size = 1.2) + 
  geom_line(aes(y = A_scaled, color = "Seagrass"), size = 1.2) + 
  scale_color_manual(values = c("Seagrass" = "#009E73", "Herbivore" = "#0072B2")) + 
  scale_y_continuous(
    name = "Herbivore density", 
    limits = herbivore_limits2, 
    sec.axis = sec_axis(~ scales::rescale(., to = seagrass_limits2), name = "Relative halo width", 
                        breaks = seq(seagrass_limits2[1], seagrass_limits2[2], by = 0.1), 
                        labels = number_format(accuracy = 0.1))
  ) +
  geom_vline(xintercept = time_points-200, linetype = "dashed", size = 0.5, color = "black") +
  geom_point(data = intersection_data, aes(x = Time-200, y = value, color = Line), size = 3) +
  labs(x = "Time (days)", y = "Density", color = "Legend") +
  ggtitle("D. Cyclic; Optimal SST = 28°C") +
  theme_minimal() +
  theme(
    axis.title.y.right = element_text(color = "#009E73"),
    axis.text.y.right = element_text(color = "#009E73"),
    axis.title.y.left = element_text(color = "#0072B2"),
    axis.text.y.left = element_text(color = "#0072B2"),
    legend.position = "none", 
    text = element_text(size = 20)
  )


overlay_data <- combined_df %>%
  filter(T_opt == 28)

overlay_data <- merge(overlay_data, sst_df, by = "Time")

#filter only for where it intersects with survey points
overlay_data <- overlay_data %>%
  filter(Time %in% time_points)

p8<-ggplot(overlay_data, aes(x=H, y=SST, color = A0 - A)) + 
  geom_point(size=3) +
  labs(x = "Herbivore density", y = "SST (°C)", color = "Relative halo width") +
  ggtitle("") +
  scale_color_viridis_c(option="viridis") +
  theme_minimal() +
  theme(text = element_text(size = 20))



######################### alternate time point sampling #####################################

filtered_data <- combined_df %>%
  filter(T_opt == 28) %>%
  filter(!is.na(A) & !is.na(H))

filtered_data$halosize <- A0 - filtered_data$A
filtered_data <- filtered_data %>%
  mutate(A_scaled = scales::rescale(halosize, to = herbivore_limits2))

# Create a sequence of time points at every 30 units
time_points2 <- seq(215, max(combined_df$Time), by = 30)

# Calculate intersection points
intersection_points <- filtered_data %>%
  filter(Time %in% time_points2) %>%
  select(Time, A_scaled, H) 

# Combine the intersection points
intersection_data <- melt(intersection_points, id.vars = "Time", variable.name = "Line")

# Change the entries in column Line
intersection_data <- intersection_data %>%
  mutate(Line = ifelse(Line == "A_scaled", "Seagrass", ifelse(Line == "H", "Herbivore", Line)))

# Create the plot
p7_extra<-ggplot(filtered_data %>% filter(Time >=200), aes(x = Time-200)) +
  geom_line(aes(y = H, color = "Herbivore"), size = 1.2) + 
  geom_line(aes(y = A_scaled, color = "Seagrass"), size = 1.2) + 
  scale_color_manual(values = c("Seagrass" = "#009E73", "Herbivore" = "#0072B2")) + 
  scale_y_continuous(
    name = "Herbivore density", 
    limits = herbivore_limits2, 
    sec.axis = sec_axis(~ scales::rescale(., to = seagrass_limits2), name = "Relative halo width", 
                        breaks = seq(seagrass_limits2[1], seagrass_limits2[2], by = 0.1), 
                        labels = number_format(accuracy = 0.1))
  ) +
  geom_vline(xintercept = time_points2-200, linetype = "dashed", size = 0.5, color = "black") +
  geom_point(data = intersection_data, aes(x = Time-200, y = value, color = Line), size = 3) +
  labs(x = "Time (days)", y = "Density", color = "Legend") +
  ggtitle("D. Cyclic; Optimal SST = 28°C") +
  theme_minimal() +
  theme(
    axis.title.y.right = element_text(color = "#009E73"),
    axis.text.y.right = element_text(color = "#009E73"),
    axis.title.y.left = element_text(color = "#0072B2"),
    axis.text.y.left = element_text(color = "#0072B2"),
    legend.position = "none", 
    text = element_text(size = 20)
  )


# Filter data for T_opt == 24 and Time in time_points
overlay_data <- combined_df %>%
  filter(T_opt == 28)

overlay_data <- merge(overlay_data, sst_df, by = "Time")

#filter only for where it intersects with survey points
overlay_data <- overlay_data %>%
  filter(Time %in% time_points2)

p8_extra<-ggplot(overlay_data, aes(x=H, y=SST, color = A0 - A)) + 
  geom_point(size=3) +
  labs(x = "Herbivore density", y = "SST (°C)", color = "Relative halo width") +
  ggtitle("") +
  scale_color_viridis_c(option="viridis") +
  theme_minimal() +
  theme(text = element_text(size = 20))


############### Figure 6
Figure6_total <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow=4, ncol=2)

ggsave("~/Desktop/Figure6.png", plot = Figure6_total, width = 16, height = 12, bg = "transparent")


#individual panels of figure 6
ggarrange(p1+rremove("xlab"),p3 +rremove("xlab"),p5 +rremove("xlab"),p7, nrow=4, ncol=1)

ggarrange(p2+rremove("xlab"),p4+rremove("xlab"),p6+rremove("xlab"),p8, nrow=4, ncol=1)


################# Figure S4, alternate time points
FigureS4 <- ggarrange(p1_extra, p2_extra, p3_extra, p4_extra, p5_extra, p6_extra, p7_extra, p8_extra, nrow=4, ncol=2)

ggsave("~/Desktop/FigureS4.png", plot = FigureS4, width = 16, height = 12, bg = "transparent")


