# Numerically integrates a stochastic version of the 
# Rosenzweig–MacArthur 2-species consumer-resource model
# (model modified from Bieg et al 2023 Proc R Soc B)

# Arguments:
# state: vector of initial abundances
# p: vector of parameters, including process noise magnitudes
# times: vector of integration times

# Output:
# out: data frame with species abundances over time

rosenzweig_macarthur_2sp_sde <- function(state, p, times) {
  # set up diffeqr
  de <- diffeqr::diffeq_setup()
  # define deterministic part
  f <- function(u, p, t) {
    du1 <- u[1] * (p[1] * (1 - (u[1] / p[2])) - (p[3] * u[2]) / (p[4] + u[1]))
    du2 <- u[2] * ((p[5] * p[3] * u[1]) / (p[4] + u[1]) - p[6])
    deterministic <- c(du1, du2)
    return(deterministic)
  }
  # define stochastic part
  g <- function(u, p, t) {
    du1 <- p[7] * u[1]
    du2 <- p[8] * u[2]
    stochastic <- c(du1, du2)
    return(stochastic)
  }
  # parameters
  p <- as.numeric(p)
  # integration time steps
  saveat <- times[2] - times[1]
  # time range
  tspan = as.numeric(range(times))
  # initial condition
  u0 <- as.numeric(state)
  # numerical integration
  prob <- de$SDEProblem(f, g, u0, tspan, p)
  sol <- de$solve(prob, alg = de$EM(), dt = saveat, saveat = saveat)
  # organize and return data frame
  out <- as.data.frame(cbind(sol$t, as.data.frame(t(sapply(sol$u, identity)))))
  names(out) <- c("time", paste("x", 1:length(state), sep = ""))
  return(out)
}
