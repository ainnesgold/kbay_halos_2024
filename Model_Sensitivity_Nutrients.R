source('Stable_Cycle_Baseline.R')
#################################### NUTRIENT DATA ####################################
nutrient_raw <- read.csv("Data/nutrients_combined.csv")

nutrient_sum <- nutrient_raw %>%
  group_by(Month, Year) %>%
  summarize(avg_total_N = mean(Total_N),
            avg_total_P = mean(Total_P),
            avg_phosphate = mean(Phosphate),
            avg_sillicate = mean(Sillicate),
            avg_NN = mean(N.N),
            avg_ammonia = mean(Ammonia),
            avg_chl = mean(Chlorophyll))

nutrient_repeated <- nutrient_sum[rep(seq_len(nrow(nutrient_sum)), each = 30), ]
nutrient_sorted <- nutrient_repeated[order(nutrient_repeated$Year, nutrient_repeated$Month),]
nutrient_sorted$avg_total_N_ug <- nutrient_sorted$avg_total_N * 28.02 #molar mass of Nitrogen, to convert total N from umol to ug
nutrient_sorted$avg_NN_ug <- nutrient_sorted$avg_NN * 54.00495 #avg of nitrate and nitrite molar mass, to convert N+N from umol to ug
nutrient_sorted$avg_ammonia_ug <- nutrient_sorted$avg_ammonia * 17.031 #ammonia molar mass, to convert ammonia from umol to ug
nutrient_sorted$inorganic_N_ug <- nutrient_sorted$avg_ammonia_ug + nutrient_sorted$avg_NN_ug

#Selecting totalN
nutrient_values <- nutrient_sorted$avg_total_N_ug

#################################### NUTRIENT DEPENDENT SEAGRASS GROWTH FUNCTION ####################################
sg_growth_nutrients <- function(slope, nutrients) {
  q <- 1 / (1.25 + exp(slope * nutrients))
  q <- ifelse(q < 0, 0, q) 
  return(q)
}

slope_values <- seq(-0.0025, -0.05, by = -0.005)

x <- seq(0, 2500, length.out = 1000)

# Calculate y values for each T_opt value
data_sg_list <- lapply(slope_values, function(slope) {
  y <- sg_growth_nutrients(x, slope)
  data.frame(x = x, y = y, slope = as.factor(slope))
})

# Combine all dataframes into one
data_sg <- do.call(rbind, data_sg_list)

# Create the plot using ggplot2
nutrient_function_plot <- ggplot(data_sg, aes(x = x, y = y, color = slope)) +
  geom_line(size = 1.5) +
  scale_color_viridis_d()+
  labs( title = "E.",
        x = "Total N (ug/L)",
        y = "Seagrass growth rate") +
  theme_minimal() +
  guides(color = guide_legend(label.wrap = 20, title="Slope parameter")) +
  theme(text = element_text(size = 20))

# Print the plot
#print(nutrient_function_plot)

###############################

#################################### MODEL SET UP ####################################


# Define the ODE system
ode_system_nutrients <- function(t, y, parms, nutrient_values2, slope) {
  nutrients <- nutrient_values2[which.min(abs(times - t))]
  with(as.list(c(y, parms)), {
    if (t <= 200) {
      q <- 0.8
    } else {
      q <- 1 / (1.25 + exp(slope * nutrients))
      q <- ifelse(q < 0, 0, q)
    }
    dA <- q * A * (1 - A / (A0 * exp(-(R * (1 / rc) * (A0 - A))))) - ((g * A * H) / (1 + g * s * A))
    dH <- (r * A * H) / (1 + g * s * A) * (1 - H / ((1 - A0) * k)) - (m * H)
    return(list(c(dA, dH, q)))
  })
}


# Define initial conditions
initial_conditions <- c(A = 0.5, H = 0.434, q=0.8)

#repeating first 200 because of the 200 timestep burn in period
nutrient_values2 <- c(nutrient_values[1:200], nutrient_values) 
times <- seq(0, length(nutrient_values2), by=1)
q_values <- numeric(length(times))


#################################### FIXED MODEL ####################################

#Define parameters
parameters <- c(A0 = 0.8, R = 0.13, rc = 2, g = 2, s = 6, r = 8, k = 5, m = 0.03)

#Sensitivity code
# Define function for sensitivity analysis
perform_sensitivity_analysis <- function(slope_values) {
  results <- list()
  for (slope in slope_values) {
    # Solve ODE system for the given T_opt
    out <- ode(y = initial_conditions, times = times, func = ode_system_nutrients, parms = parameters, 
               nutrient_values = nutrient_values2, slope = slope)
    # Store results
    results[[as.character(slope)]] <- out
  }
  return(results)
}

# Test sensitivity analysis 
sensitivity_results <- perform_sensitivity_analysis(slope_values)

# Initialize an empty list to store dataframes for each T_opt
dataframes <- list()

# Create a dataframe for each T_opt and store it in the list
for (slope in names(sensitivity_results)) {
  out <- sensitivity_results[[slope]]
  df <- data.frame(Time = out[, 1], A = out[, 2], H = out[, 3], q = out[, 4])
  df$slope <- as.numeric(slope)
  dataframes[[slope]] <- df
}

# Combine all dataframes into one
combined_df <- do.call(rbind, dataframes)

# Organize combined dataframe by T_opt
combined_df <- combined_df[order(combined_df$slope, combined_df$Time), ]

#calculate the q (seagrass growth rate) at each time step
q_values <- numeric(length(combined_df$Time))
q_values[1] = combined_df$q[1]
for (i in 2:length(combined_df$Time)) {
  new_value <- combined_df$q[i] - combined_df$q[i-1]
  q_values[i] <- new_value
}
combined_df$q <- q_values

#seagrass plot
a_gg <- ggplot(combined_df %>% filter(Time >= 200), aes(x=Time-200, y=A, color=factor(slope, levels = rev(levels(factor(slope)))))) +
  geom_line(size=1) +
  geom_line(data = solution_stable %>% filter(time >=200 & time < 620), aes(x = time-200, y = A_stable_baseline), 
            color = "black", linetype = "dashed", size = 1) +
  scale_color_viridis_d() +
  theme_minimal() +
  labs(x="Time (days)", y = "Seagrass density") +
  theme(legend.position = "right", 
        text = element_text(size = 20),
        legend.key.height = unit(1.5, "lines")) +  # Adjust legend key height for more space
  guides(color = guide_legend(label.wrap = 20, title="Nutrient function slope")) +
  ggtitle("A. Fixed baseline (R = 0.13)") +
  scale_y_continuous(limits = c(0, 0.6))

#herbivore plot
h_gg<-ggplot(combined_df %>% filter(Time >= 200), aes(x=Time-200, y=H, color=factor(slope, levels = rev(levels(factor(slope)))))) +
  geom_line(size=1) +
  geom_line(data = solution_stable %>% filter(time >= 200 & time < 620), aes(x = time-200, y = H_stable_baseline), 
            color = "black", linetype = "dashed", size = 1) +
  scale_color_viridis_d() +
  theme_minimal() +
  labs(x="Time (days)", y = "Herbivore density") +
  theme(legend.position = "right", 
        text = element_text(size = 20),
        legend.key.height = unit(1.5, "lines")) +  # Adjust legend key height for more space
  guides(color = guide_legend(label.wrap = 20, title="Nutrient function slope")) +
  ggtitle("") +
  scale_y_continuous(limits = c(0, 1))

#plot of seagrass growth rates
q_gg_nutrients<-ggplot(combined_df %>% filter(Time >= 200), aes(x=Time-200, y=q, color=factor(slope, levels = rev(levels(factor(slope)))))) +
  geom_line(size=1) +
  scale_color_viridis_d() +
  theme_minimal() +
  labs(x="Time (days)", y = "Seagrass growth rate") +
  theme(legend.position = "right", 
        text = element_text(size = 20),
        legend.key.height = unit(1.5, "lines")) +  # Adjust legend key height for more space
  guides(color = guide_legend(label.wrap = 20, title="Nutrient function slope")) +
  ggtitle("F.") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1)



#################################### CYCLIC BASELINE MODEL ####################################

# Define initial conditions
initial_conditions <- c(A = 0.5, H = 0.434, q=0.8)
times <- seq(0, length(nutrient_values2), by=1)
q_values <- numeric(length(times))

#Define parameters
parameters <- c(A0 = 0.8, R = 1.28, rc = 2, g = 2, s = 6, r = 8, k = 5, m = 0.03)
sensitivity_results <- perform_sensitivity_analysis(slope_values)
# Initialize an empty list to store dataframes for each T_opt
dataframes <- list()

# Create a dataframe for each T_opt and store it in the list
for (slope in names(sensitivity_results)) {
  out <- sensitivity_results[[slope]]
  df <- data.frame(Time = out[, 1], A = out[, 2], H = out[, 3], q = out[, 4])
  df$slope <- as.numeric(slope)
  dataframes[[slope]] <- df
}

# Combine all dataframes into one
combined_df <- do.call(rbind, dataframes)

# Organize combined dataframe by T_opt
combined_df <- combined_df[order(combined_df$slope, combined_df$Time), ]

q_values <- numeric(length(combined_df$Time))
q_values[1] = combined_df$q[1]
for (i in 2:length(combined_df$Time)) {
  new_value <- combined_df$q[i] - combined_df$q[i-1]
  q_values[i] <- new_value
}
combined_df$q <- q_values


a_gg2 <-ggplot(combined_df %>% filter(Time > 200), aes(x=Time-200, y=A, color=factor(slope, levels = rev(levels(factor(slope)))))) +
  geom_line(size=1) +
  geom_line(data = solution_cycles %>% filter(time > 200 & time < 620), 
            aes(x = time-200, y = A_cycles_baseline), color = "black", linetype = "dashed", size = 1) +
  scale_color_viridis_d() +
  theme_minimal() +
  labs(x="Time (days)", y = "Seagrass density") +
  theme(legend.position = "right", 
        text = element_text(size = 20),
        legend.key.height = unit(1.5, "lines")) +  # Adjust legend key height for more space
  guides(color = guide_legend(label.wrap = 20, title="Nutrient function slope")) +
  ggtitle("B. Cyclic baseline (R = 1.28)") +
  scale_y_continuous(limits = c(0, 0.6))



h_gg2<-ggplot(combined_df %>% filter(Time > 200), aes(x=Time-200, y=H, color=factor(slope, levels = rev(levels(factor(slope)))))) +
  geom_line(size=1) +
  geom_line(data = solution_cycles %>% filter(time > 200 & time < 620), 
            aes(x = time-200, y = H_cycles_baseline), color = "black", linetype = "dashed", size = 1) +
  scale_color_viridis_d() +
  theme_minimal() +
  labs(x="Time (days)", y = "Herbivore density") +
  theme(legend.position = "right", 
        text = element_text(size = 20),
        legend.key.height = unit(1.5, "lines")) +  # Adjust legend key height for more space
  guides(color = guide_legend(label.wrap = 20, title="Nutrient function slope")) +
  ggtitle("") +
  scale_y_continuous(limits = c(0, 1))

#Figure S3
#ggarrange(a_gg, h_gg, a_gg2, h_gg2, nrow=2, ncol=2, common.legend = TRUE, legend="right") #12 x 8


#Plot the nutrient time series
# Create a data frame with your data
nutrient_df <- data.frame(Time = times[-1], Nutrient = nutrient_values2)

# Create the ggplot object
nutrient_plot<-ggplot(nutrient_df %>% filter(Time >200), aes(x = Time-200, y = Nutrient)) +
  geom_line() +
  labs(x = "Time (days)", y = "Total N (ug/L)") +
  theme_minimal() +
  ggtitle("D.")+
  theme(text = element_text(size = 20))

#ggarrange(nutrient_plot, nutrient_function_plot, q_gg_nutrients, nrow=1, ncol=3, common.legend = TRUE, legend="right")










