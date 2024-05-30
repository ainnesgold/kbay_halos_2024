##Baseline (no temp, nutrients) model from Ong et al.

library(deSolve)
library(tidyverse)
library(ggpubr)

#Halo model
ode_system <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    dA <- q * A * (1 - A / (A0 * exp(-(R * (1 / rc) * (A0 - A))))) - ((g * A * H) / (1 + g * s * A))
    dH <- (r * A * H) / (1 + g * s * A) * (1 - H / ((1 - A0) * k)) - (m * H)
    return(list(c(dA, dH)))
  })
}

# Set the time points where you want to evaluate the solution
times <- seq(0, 641, by = 1) #same length as SST model (200 burn in + 441 SST time series points)

# Define initial conditions
initial_conditions <- c(A = 0.5, H = 0.434)


#STABLE
# Define parameter values
parameters <- c(q = 0.8, A0 = 0.8, R = 0.13, rc = 2, g = 2, s = 6, r = 8, k = 5, m = 0.03)  #original stable

# Solve the system of differential equations
solution_stable <- ode(y = initial_conditions, times = times, func = ode_system, parms = parameters)

p1<-ggplot(data=as.data.frame(solution_stable), aes(x=time, y=A)) +
  geom_line(col="dark green") +
  theme_minimal() +
  labs(x="Time", y = "Seagrass Density")

p2<-ggplot(data=as.data.frame(solution_stable), aes(x=time, y=H)) +
  geom_line(col="dark blue") +
  theme_minimal() +
  labs(x="Time", y = "Herbivore Density")
#ggarrange(p1, p2, nrow=2)


##CYCLES
parameters <- c(q = 0.8, A0 = 0.8, R = 1.28, rc = 2, g = 2, s = 6, r = 8, k = 5, m = 0.03)  #original cycle

# Solve the system of differential equations
solution_cycles <- ode(y = initial_conditions, times = times, func = ode_system, parms = parameters)

p3<-ggplot(data=as.data.frame(solution_cycles), aes(x=time, y=A)) +
  geom_line(col="dark green") +
  theme_minimal() +
  labs(x="Time", y = "Seagrass Density")

p4<-ggplot(data=as.data.frame(solution_cycles), aes(x=time, y=H)) +
  geom_line(col="dark blue") +
  theme_minimal() +
  labs(x="Time", y = "Herbivore Density")

#ggarrange(p3, p4, nrow=2)

solution_stable <- as.data.frame(solution_stable)
solution_cycles <- as.data.frame(solution_cycles)

colnames(solution_stable)[colnames(solution_stable) == "A"] <- "A_stable_baseline"
colnames(solution_stable)[colnames(solution_stable) == "H"] <- "H_stable_baseline"

colnames(solution_cycles)[colnames(solution_cycles) == "A"] <- "A_cycles_baseline"
colnames(solution_cycles)[colnames(solution_cycles) == "H"] <- "H_cycles_baseline"

solution_baseline <- merge(solution_stable, solution_cycles, by = "time")
solution_baseline <- pivot_longer(solution_baseline, A_stable_baseline:H_cycles_baseline, names_to = "type", values_to = "value")

