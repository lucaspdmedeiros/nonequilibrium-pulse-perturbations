# Iterates the discrete-time Larvae-Pupae-Adults (LPA) model
# (model from Costantino et al 1997 Science)

# Arguments:
# t_max: number of time steps
# x: vector of initial abundances
# p: vector of parameters 

# Output:
# results_df: data frame with abundances over time

larvae_pupae_adults_3sp_map <- function(t_max, x, p) {
  # initial condition
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  # parameters
  b = p[1]
  c_el = p[2] 
  c_ea = p[3]
  u_l = p[4]
  c_pa = p[5]
  u_a = p[6]
  # iterate map
  for (t in 1:(t_max-1)) {
    x1[t+1] <- b * x3[t] * exp(-c_el * x1[t] - c_ea * x3[t])
    x2[t+1] <- (1 - u_l) * x1[t]
    x3[t+1] <- x2[t] * exp(-c_pa * x3[t]) + x3[t] * (1 - u_a)
  }
  # return results data frame
  results_df <- data.frame(time = 0:(t_max-1), x1 = x1, x2 = x2, x3 = x3)
  return(results_df)
}
