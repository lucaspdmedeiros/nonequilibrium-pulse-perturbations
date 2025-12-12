# Returns the rates of change for the Rosenzweig–MacArthur 2-species
# consumer-resource model with time-varying carrying capacity 
# (model modified from Bieg et al 2023 Proc R Soc B)

# Arguments:
# t: vector of time steps
# state: vector of initial abundance values
# parameters: vector of parameters 

# Output:
# dx1: resource rate of change 
# dx2: consumer rate of change 

rosenzweig_macarthur_2sp_ode <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # rate of change
    k_force <- A * sin(p * 2 * pi * t)
    dx1 <- x1 * (r * (1 - (x1 / (k_mean + k_force))) - (a * x2) / (b + x1))
    dx2 <- x2 * ((e * a * x1) / (b + x1) - d)
    # return the rate of change
    list(c(dx1, dx2))
  })
}
