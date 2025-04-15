# Define parameters
parameters <- c(A0 = 0.8, rc = 2, g = 2, s = 6, r = 8, k = 5, m = 0.03)

# Function to compute the equation
f <- function(A, A0, rc, R) {
  1 - A / (A0 * exp(-(R * (1 / rc) * (A0 - A))))
}

# Range of A values
A_values <- seq(0, parameters["A0"], length.out = 100)

# R values to plot
R_values <- c(0.13, 1.28)

# Plot setup
plot(A_values, f(A_values, parameters["A0"], parameters["rc"], R_values[1]), 
     type = "l", col = "blue", ylim = c(0, 1), 
     ylab = "Term", xlab = "A", lwd = 2)

# Add second curve
lines(A_values, f(A_values, parameters["A0"], parameters["rc"], R_values[2]), col = "red", lwd = 2)

# Add legend
legend("topright", legend = paste("R =", R_values), col = c("blue", "red"), lwd = 2)
