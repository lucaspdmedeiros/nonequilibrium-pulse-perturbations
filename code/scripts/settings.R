# Parameter and simulation settings for different population dynamics models

# Loading packages ------------------------------
if (!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
if (!require(rootSolve)) {install.packages("rootSolve"); library(rootSolve)}
if (!require(plyr)) {install.packages("plyr"); library(plyr)}
if (!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if (!require(ggforce)) {install.packages("ggforce"); library(ggforce)}
if (!require(reshape2)) {install.packages("reshape2"); library(reshape2)}
if (!require(viridis)) {install.packages("viridis"); library(viridis)}
if (!require(ggsci)) {install.packages("ggsci"); library(ggsci)}
if (!require(expm)) {install.packages("expm"); library(expm)}
if (!require(scales)) {install.packages("scales"); library(scales)}
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(numDeriv)) {install.packages("numDeriv"); library(numDeriv)}
if (!require(gganimate)) {install.packages("gganimate"); library(gganimate)}
if (!require(gifski)) {install.packages("gifski"); library(gifski)}
if (!require(plotly)) {install.packages("plotly"); library(plotly)}
if (!require(latex2exp)) {install.packages("latex2exp"); library(latex2exp)}
if (!require(RColorBrewer)) {install.packages("RColorBrewer"); library(RColorBrewer)}
if (!require(diffeqr)) {install.packages("diffeqr"); library(diffeqr)}
if (!require(JuliaCall)) {install.packages("JuliaCall"); library(JuliaCall)}
if (!require(pbapply)) {install.packages("pbapply"); library(pbapply)}
if (!require(glmnet)) {install.packages("glmnet"); library(glmnet)}
if (!require(mvtnorm)) {install.packages("mvtnorm"); library(mvtnorm)}

# Rosenzweig–MacArthur model with fixed point ------------------------------
if (func_name == "rosenzweig_macarthur_2sp_fixed_point") {
  time_type <- "continuous"
  n_sp <- 2
  func <- rosenzweig_macarthur_2sp_ode
  A <- 0
  p <- 0.1
  r <- 5
  a <- 1.3
  b <- 1
  e <- 0.7
  d <- 0.2
  k_hopf <- (-b * (d + a * e) / (d - a * e))
  k_mean <- 1
  parms <- c(A = A, p = p, r = r, k_mean = k_mean, a = a, b = b, e = e, d = d)
  equilibrium <- c(x1 = - (b * d) / (d - a * e), x2 = - (b * e * (b * d + d * k_mean - a * e * k_mean) * r) / 
                     (k_mean * (d - a * e)^2))
  state <- equilibrium
  n_points <- 10000
  time_step <- 10 / n_points 
  times_pert <- seq(0, time_step * n_points, by = time_step)
  lim_x1 <- c(0.235, 0.325)
  lim_x2 <- c(3.47, 3.605)
  k <- k_mean + (A * sin(p * 2 * pi * times_pert))
  sample_size <- 20
}

# Rosenzweig–MacArthur model with cycle driven by time-varying carrying capacity ------------------------------
if (func_name == "rosenzweig_macarthur_2sp_forced_fixed_point") {
  time_type <- "continuous"
  n_sp <- 2
  func <- rosenzweig_macarthur_2sp_ode
  A <- 0.5
  p <- 0.1
  r <- 5
  a <- 1.3
  b <- 1
  e <- 0.7
  d <- 0.2
  k_hopf <- (-b * (d + a * e) / (d - a * e))
  k_mean <- 1
  parms <- c(A = A, p = p, r = r, k_mean = k_mean, a = a, b = b, e = e, d = d)
  equilibrium <- c(x1 = - (b * d) / (d - a * e), x2 = - (b * e * (b * d + d * k_mean - a * e * k_mean) * r) / 
                     (k_mean * (d - a * e)^2))
  state <- c(x1 = 0.3418031, x2 = 2.747028)
  n_recurrences <- 2
  n_points <- 10000
  time_step <- 10 / n_points # dominant period divided by n_points
  times <- seq(0, time_step * n_points * n_recurrences, by = time_step)
  time_step_pert <- time_step
  times_pert <- seq(0, time_step_pert * n_points, by = time_step_pert)
  lim_x1 <- c(0.07, 0.69)
  lim_x2 <- c(2.7, 4.6)
  k <- k_mean + (A * sin(p * 2 * pi * times))
  sample_size <- 20
}

# Rosenzweig–MacArthur model with limit cycle ------------------------------
if (func_name == "rosenzweig_macarthur_2sp_limit_cycle") {
  time_type <- "continuous"
  n_sp <- 2
  func <- rosenzweig_macarthur_2sp_ode
  A <- 0
  p <- 0.1
  r <- 5
  a <- 1.3
  b <- 1
  e <- 0.7
  d <- 0.2
  k_mean <- 1.8
  k_hopf <- (-b * (d + a * e) / (d - a * e))
  parms <- c(A = A, p = p, r = r, k_mean = k_mean, a = a, b = b, e = e, d = d)
  equilibrium <- c(x1 = - (b * d) / (d - a * e), x2 = - (b * e * (b * d + d * k_mean - a * e * k_mean) * r) / 
                     (k_mean * (d - a * e)^2))
  state <- c(x1 = 0.51113819, x2 = 4.905690)
  n_recurrences <- 2
  n_points <- 10000
  time_step <- 8.5 / n_points # dominant period divided by n_points
  times <- seq(0, time_step * n_points * n_recurrences, by = time_step)
  time_step_pert <- time_step
  times_pert <- seq(0, time_step_pert * n_points, by = time_step_pert)
  lim_x1 <- c(0, 1)
  lim_x2 <- c(2.8, 5.2)
  k <- k_mean + (A * sin(p * 2 * pi * times))
  sample_size <- 20
}

# Stochastic Rosenzweig–MacArthur model with limit cycle ------------------------------
if (func_name == "rosenzweig_macarthur_2sp_limit_cycle_stochastic") {
  time_type <- "continuous"
  n_sp <- 2
  func <- rosenzweig_macarthur_2sp_sde
  parms <- c(r = 5, k = 1.8, a = 1.3, b = 1, e = 0.7, d = 0.2, s1 = 0.04, s2 = 0.04)
  state <- c(x1 = 0.51113819, x2 = 4.905690)
  n_recurrences <- 10
  n_points <- 100000
  time_step <- 8.5 / n_points # dominant period divided by n_points
  times <- seq(0, time_step * n_points * n_recurrences, by = time_step)
  lim_x1 <- c(0, 1.15)
  lim_x2 <- c(2.4, 5.7)
}

# Hastings-Powell model with limit cycle and time-varying attack rate ------------------------------
if (func_name == "hastings_powell_3sp_forced_cycle") {
  time_type <- "continuous"
  n_sp <- 3
  func <- hastings_powell_3sp_ode
  parms <- c(r = 0.000048, k = 0.99, a1 = 0.8036, a2_base = 0.1984, 
             e1 = 1, e2 = 1, b1 = 0.16129, b2 = 0.5, d1 = 0.4, d2 = 0.08)
  state <- c(x1 = 0.6538001, x2 = 0.3563991, x3 = 0.9953945)
  n_recurrences <- 8
  n_points <- 10000
  time_step <- 50 / n_points # dominant period divided by n_points
  times <- seq(0, time_step * n_points * n_recurrences, by = time_step)
  time_step_pert <- time_step
  times_pert <- seq(0, time_step_pert * n_points, by = time_step_pert)
  a2 <- parms[4] + parms[1] * times
  lim_x1 <- c(0.3, 0.9)
  lim_x2 <- c(0.15, 0.5)
  lim_x3 <- c(0.75, 1.1)
  sample_size <- 20
}

# Hastings-Powell model with chaotic attractor ------------------------------
if (func_name == "hastings_powell_3sp_chaos") {
  time_type <- "continuous"
  n_sp <- 3
  func <- hastings_powell_3sp_ode
  parms <- c(r = 0, k = 0.99, a1 = 0.8036, a2_base = 0.23008, 
             e1 = 1, e2 = 1, b1 = 0.16129, b2 = 0.5, d1 = 0.4, d2 = 0.08)
  state <- c(x1 = 0.2445642, x2 = 0.2500377, x3 = 0.8477652)
  n_recurrences <- 10
  n_points <- 10000
  time_step <- 42 / n_points # dominant period divided by n_points
  times <- seq(0, time_step * n_points * n_recurrences, by = time_step)
  time_step_pert <- time_step
  times_pert <- seq(0, time_step_pert * n_points, by = time_step_pert)
  a2 <- parms[6] + parms[1] * times
  lim_x1 <- c(0.1, 0.95)
  lim_x2 <- c(0.1, 0.6)
  lim_x3 <- c(0.5, 1.1)
  sample_size <- 20
}

# Stochastic Hastings-Powell model with chaotic attractor ------------------------------
if (func_name == "hastings_powell_3sp_chaos_stochastic") {
  time_type <- "continuous"
  n_sp <- 3
  func <- hastings_powell_3sp_sde
  parms <- c(k = 0.99, a1 = 0.8036, b1 = 0.16129,
             e1 = 1, d1 = 0.4, a2 = 0.23008, 
             b2 = 0.5, e2 = 1, d2 = 0.08,
             s1 = 0.01, s2 = 0.01, s3 = 0.01)
  state <- c(x1 = 0.2445642, x2 = 0.2500377, x3 = 0.8477652)
  n_recurrences <- 10
  n_points <- 100000
  time_step <- 42 / n_points # dominant period divided by n_points
  times <- seq(0, time_step * n_points * n_recurrences, by = time_step)
  lim_x1 <- c(0.15, 0.95)
  lim_x2 <- c(0.1, 0.6)
  lim_x3 <- c(0.55, 1.15)
}

# Discrete-time predator-prey model with limit cycle ------------------------------
if (func_name == "predator_prey_2sp_limit_cycle") {
  time_type <- "discrete"
  n_sp <- 2
  func <- predator_prey_2sp_map
  a = 2
  b = 2
  c = 2
  d = 1.85
  alpha = 0.1
  beta = 0.1
  tau = 1.1
  parms <- c(a = a, b = b, c = c, d = d, 
             alpha = alpha, beta = beta, tau = tau)
  state <- c(x1 = 1.128124, x2 = 0.2990047)
  n_recurrences <- 4
  n_points <- 5
  lim_x1 <- c(1, 1.7)
  lim_x2 <- c(0.2, 0.6)
}

# Discrete-time Larvae-Pupae-Adults (LPA) model with chaotic attractor ------------------------------
if (func_name == "larvae_pupae_adults_3sp_chaos") {
  time_type <- "discrete"
  n_sp <- 3
  func <- larvae_pupae_adults_3sp_map
  b <- 6.598
  c_el <- 0.01209
  c_ea <- 0.01155
  u_l <- 0.2055
  c_pa <- 0.35
  u_a <- 0.96
  parms <- c(b = b, c_el = c_el, c_ea = c_ea, u_l = u_l, 
             c_pa = c_pa, u_a = u_a)
  state <- c(x1 = 8.067926, x2 = 1.234136, x3 = 57.151993)
  n_recurrences <- 6
  n_points <- 3
  lim_x1 <- c(1, 192)
  lim_x2 <- c(1, 153)
  lim_x3 <- c(0, 76)
}
