# Iterates the discrete-time predator-prey 
# model from Zhang et al 2018 Discret Dyn Nat Soc

# Arguments:
# t_max: number of time steps
# x: vector of initial abundances
# p: vector of parameters 

# Output:
# results_df: data frame with abundances over time

predator_prey_2sp_map <- function(t_max, x, p) {
  # initial condition
  x1 <- x[1]
  x2 <- x[2]
  # parameters
  a = p[1]
  b = p[2] 
  c = p[3]
  d = p[4]
  alpha = p[5]
  beta = p[6]
  tau = p[7]
  # iterate model
  for (t in 1:(t_max-1)) {
    x1[t+1] <- x1[t] + tau * x1[t] * (a - x1[t] - (b * x2[t]) / ((1 + alpha * x1[t]) * (1 + beta * x2[t])))
    x2[t+1] <- x2[t] + tau * x2[t] * (-c + (d * x1[t]) / ((1 + alpha * x1[t]) * (1 + beta * x2[t])))
  }
  # return results data frame
  results_df <- data.frame(time = 0:(t_max-1), x1 = x1, x2 = x2)
  return(results_df)
}
