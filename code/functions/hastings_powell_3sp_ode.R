# Returns the rates of change for the Hastings-Powell 3-species
# food chain model with a time-varying predator attack rate
# (model modified from Hastings & Powell 1991 Ecology)

# Arguments:
# t: vector of time steps
# state: vector of initial abundance values
# parameters: vector of parameters 

# Output:
# dx1: resource rate of change 
# dx2: consumer rate of change 
# dx3: predator rate of change 

hastings_powell_3sp_ode <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # rate of change
    a2_force <- r * t 
    dx1 <- x1 * (1 - (x1 / k)) - (a1 * x1 * x2) / (x1 + b1)
    dx2 <- (a1 * e1 * x1 * x2) / (x1 + b1) - d1 * x2 - ((a2_base + a2_force) * x2 * x3) / (x2 + b2)
    dx3 <- ((a2_base + a2_force) * e2 * x2 * x3) / (x2 + b2) - d2 * x3
    # return the rate of change
    list(c(dx1, dx2, dx3))
  })
}