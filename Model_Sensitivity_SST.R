source('Stable_Cycle_Baseline.R')
#################################### TEMPERATURE DATA ####################################
temp_raw <- read.csv("Data/Temp_combined.csv")
temp_raw$avg_temp <- rowMeans(temp_raw[, c("Temp_S1", "Temp_S2", "Temp_S3")], na.rm = TRUE)
temp_raw$avg_light <- rowMeans(temp_raw[, c("Intensity_S1", "Intensity_S2", "Intensity_S3")], na.rm = TRUE)

#temp logs every 6 hours. take the average of every 4 times to get a daily avg
# Function to take the average of every 4 rows for multiple columns
average_of_four_rows <- function(data) {
  n <- nrow(data)
  num_groups <- ceiling(n / 4)
  result <- matrix(NA, nrow = num_groups, ncol = ncol(data))
  for (i in 1:num_groups) {
    start_row <- (i - 1) * 4 + 1
    end_row <- min(i * 4, n)
    result[i, ] <- colMeans(data[start_row:end_row, ], na.rm = TRUE)
  }
  result_df <- as.data.frame(result)
  return(result_df)
}

# Applying the function to the dataframe
averages <- average_of_four_rows(temp_raw %>% select(Month, Year, avg_temp))

#convert to Celsius
sst_values <- (5/9) * (averages$V3 - 32)

#################################### TEMP DEPENDENT SEAGRASS GROWTH FUNCTION ####################################
# Define seagrass_growth function with T_opt as an argument
seagrass_growth <- function(sst, T_opt) {
  q <- 0.8 + 0 * sst + -0.0037 * (sst - T_opt)^2
  q <- ifelse(q < 0, 0, q) 
  return(q)
}

x <- seq(5, 45, by = 0.1)
T_opt_values <- 20:30

data_sg_list <- lapply(T_opt_values, function(T_opt) {
  y <- seagrass_growth(x, T_opt)
  data.frame(x = x, y = y, T_opt = as.factor(T_opt))
})

data_sg <- do.call(rbind, data_sg_list)

sst_function_plot <- ggplot(data_sg %>% filter(T_opt == 20 | T_opt == 24 | T_opt == 28), aes(x = x, y = y, color = T_opt)) +
  geom_line(size = 1.5) +
  scale_color_viridis_d()+
  labs(title = "B.",
       x = "SST (°C)",
       y = "Vegetation growth rate") +
  theme_minimal() +
  guides(color = guide_legend(label.wrap = 20, title="Optimal SST")) +
  theme(text = element_text(size = 20))



#################################### MODEL SET UP ####################################

ode_system_sst <- function(t, y, parms, sst_values2, T_opt) {
  sst <- sst_values2[which.min(abs(times - t))]  
  with(as.list(c(y, parms)), {
    if (t <= 200) {
      q <- 0.8
    } else {
      q <- 0.8 - 0.0037 * (sst - T_opt)^2
      q <- ifelse(q < 0, 0, q)
    }
    dA <- q * A * (1 - A / (A0 * exp(-(R * (1 / rc) * (A0 - A))))) - ((g * A * H) / (1 + g * s * A))
    dH <- (r * A * H) / (1 + g * s * A) * (1 - H / ((1 - A0) * k)) - (m * H)
    return(list(c(dA, dH, q)))
  })
}


initial_conditions <- c(A = 0.5, H = 0.434, q=0.8)
sst_values2 <- c(sst_values[1:200], sst_values) #repeating first 30 since it won't use them the first time when there's no temp dependent q
times <- seq(0, length(sst_values2), by=1)
q_values <- numeric(length(times))


#################################### STABLE MODEL ####################################

parameters <- c(A0 = 0.8, R = 0.13, rc = 2, g = 2, s = 6, r = 8, k = 5, m = 0.03)

perform_sensitivity_analysis <- function(T_opt_values) {
  results <- list()
  for (T_opt in T_opt_values) {
    # Solve ODE system for the given T_opt
    out <- ode(y = initial_conditions, times = times, func = ode_system_sst, parms = parameters, sst_values = sst_values2, T_opt = T_opt)
    # Store results
    results[[as.character(T_opt)]] <- out
  }
  return(results)
}

T_opt_values <- 20:30
sensitivity_results <- perform_sensitivity_analysis(T_opt_values)

dataframes <- list()

for (T_opt in names(sensitivity_results)) {
  out <- sensitivity_results[[T_opt]]
  df <- data.frame(Time = out[, 1], A = out[, 2], H = out[, 3], q = out[, 4])
  df$T_opt <- as.numeric(T_opt)
  dataframes[[T_opt]] <- df
}

combined_df <- do.call(rbind, dataframes)

combined_df <- combined_df[order(combined_df$T_opt, combined_df$Time), ]

q_values <- numeric(length(combined_df$Time))
q_values[1] = combined_df$q[1]
for (i in 2:length(combined_df$Time)) {
  new_value <- combined_df$q[i] - combined_df$q[i-1]
  q_values[i] <- new_value
}
combined_df$q <- q_values


A0 <- 0.8

a_gg <-ggplot(combined_df %>% filter(Time >=200) %>% filter(T_opt == 20 | T_opt == 24 | T_opt == 28),
              aes(x=Time-200, y=A0 - A, color=as.factor(T_opt))) +
  geom_line(size=1) +
  geom_line(data = solution_stable %>% 
              filter(time >=200), aes(x = time-200, y = A0 - A_stable_baseline), 
            color = "black", linetype = "dashed", size = 1) +
  scale_color_viridis_d() +
  theme_minimal() +
  labs(x="Time (days)", y = "Relative halo width") +
  theme(legend.position = "right", 
        text = element_text(size = 20),
        legend.key.height = unit(1.5, "lines")) +  # Adjust legend key height for more space
  guides(color = guide_legend(label.wrap = 20, title="Optimal SST")) +
  ggtitle("A. Fixed baseline (R = 0.13)") +
scale_y_continuous(limits = c(0.2, 0.8))



h_gg<-ggplot(combined_df %>% filter(Time >=200) %>% filter(T_opt == 20 | T_opt == 24 | T_opt == 28), 
             aes(x=Time-200, y=H, color=as.factor(T_opt))) +
  geom_line(size=1) +
  geom_line(data = solution_stable %>% filter(time >=200), aes(x = time-200, y = H_stable_baseline), color = "black", linetype = "dashed", size = 1) +
  scale_color_viridis_d() +
  theme_minimal() +
  labs(x="Time (days)", y = "Herbivore density") +
  theme(legend.position = "right", 
        text = element_text(size = 20),
        legend.key.height = unit(1.5, "lines")) +  # Adjust legend key height for more space
  guides(color = guide_legend(label.wrap = 20, title="Optimal SST")) +
  ggtitle("") +
  scale_y_continuous(limits = c(0, 1))



q_gg_sst<-ggplot(combined_df %>% filter(Time >=200) %>% filter(T_opt == 20 | T_opt == 24 | T_opt == 28),
                 aes(x=Time-200, y=q, color=as.factor(T_opt))) +
  geom_line(size=1) +
  scale_color_viridis_d() +
  theme_minimal() +
  labs(x="Time (days)", y = "Vegetation growth rate") +
  theme(legend.position = "right", 
        text = element_text(size = 20),
        legend.key.height = unit(1.5, "lines")) +  # Adjust legend key height for more space
  guides(color = guide_legend(label.wrap = 20, title="Optimal SST")) +
  ggtitle("C.") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1)


#################################### CYCLIC MODEL ####################################

initial_conditions <- c(A = 0.5, H = 0.434, q=0.8)
times <- seq(0, length(sst_values2), by=1)
q_values <- numeric(length(times))

parameters <- c(A0 = 0.8, R = 1.28, rc = 2, g = 2, s = 6, r = 8, k = 5, m = 0.03)

sensitivity_results <- perform_sensitivity_analysis(T_opt_values)

dataframes <- list()

for (T_opt in names(sensitivity_results)) {
  out <- sensitivity_results[[T_opt]]
  df <- data.frame(Time = out[, 1], A = out[, 2], H = out[, 3], q = out[, 4])
  df$T_opt <- as.numeric(T_opt)
  dataframes[[T_opt]] <- df
}

combined_df <- do.call(rbind, dataframes)

combined_df <- combined_df[order(combined_df$T_opt, combined_df$Time), ]

q_values <- numeric(length(combined_df$Time))
q_values[1] = combined_df$q[1]
for (i in 2:length(combined_df$Time)) {
  new_value <- combined_df$q[i] - combined_df$q[i-1]
  q_values[i] <- new_value
}
combined_df$q <- q_values

a_gg2 <-ggplot(combined_df %>% filter(Time >=200) %>% filter(T_opt == 20 | T_opt == 24 | T_opt == 28),
               aes(x=Time-200, y=A0 - A, color=as.factor(T_opt))) +
  geom_line(size=1) +
  geom_line(data = solution_cycles %>% filter(time >=200), aes(x = time-200, y = A0 - A_cycles_baseline), color = "black", linetype = "dashed", size = 1) +
  scale_color_viridis_d() +
  theme_minimal() +
  labs(x="Time (days)", y = "Relative halo width") +
  theme(legend.position = "right", 
        text = element_text(size = 20),
        legend.key.height = unit(1.5, "lines")) +  # Adjust legend key height for more space
  guides(color = guide_legend(label.wrap = 20, title="Optimal SST")) +
  ggtitle("B. Cyclic baseline (R = 1.28)") +
  scale_y_continuous(limits = c(0.2, 0.8))


h_gg2<-ggplot(combined_df %>% filter(Time >=200) %>% filter(T_opt == 20 | T_opt == 24 | T_opt == 28),
              aes(x=Time-200, y=H, color=as.factor(T_opt))) +
  geom_line(size=1) +
  geom_line(data = solution_cycles %>% filter(time >=200), aes(x = time-200, y = H_cycles_baseline), color = "black", linetype = "dashed", size = 1) +
  scale_color_viridis_d() +
  theme_minimal() +
  labs(x="Time (days)", y = "Herbivore density") +
  theme(legend.position = "right", 
        text = element_text(size = 20),
        legend.key.height = unit(1.5, "lines")) +  # Adjust legend key height for more space
  guides(color = guide_legend(label.wrap = 20, title="Optimal SST")) +
  ggtitle("") +
  scale_y_continuous(limits = c(0, 1))


# PLOT SST DATA
sst_df <- data.frame(Time = times[-1], SST = sst_values2)

sst_plot<-ggplot(sst_df %>% filter(Time > 200), aes(x = Time-200, y = SST)) +
  geom_line() +
  labs(x = "Time (days)", y = "SST (°C)") +
  theme_minimal() +
  ggtitle("A.")+
  theme(text = element_text(size = 20))

