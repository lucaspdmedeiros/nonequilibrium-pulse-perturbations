# Computes the power spectrum for a population dynamics model 
# from time series generate from the model

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/rosenzweig_macarthur_2sp_ode.R")
source("code/functions/hastings_powell_3sp_ode.R")
if (!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
if (!require(rootSolve)) {install.packages("rootSolve"); library(rootSolve)}
if (!require(plyr)) {install.packages("plyr"); library(plyr)}
if (!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if (!require(ggforce)) {install.packages("ggforce"); library(ggforce)}
if (!require(latex2exp)) {install.packages("latex2exp"); library(latex2exp)}
if (!require(reshape2)) {install.packages("reshape2"); library(reshape2)}
if (!require(viridis)) {install.packages("viridis"); library(viridis)}
if (!require(ggsci)) {install.packages("ggsci"); library(ggsci)}
if (!require(expm)) {install.packages("expm"); library(expm)}
if (!require(scales)) {install.packages("scales"); library(scales)}
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(numDeriv)) {install.packages("numDeriv"); library(numDeriv)}
if (!require(gganimate)) {install.packages("gganimate"); library(gganimate)}
if (!require(gifski)) {install.packages("gifski"); library(gifski)}

# settings ------------------------------
# to reproduce results 
set.seed(42)
# whether to save plots
save_plots <- FALSE
# species to use
species <- "x1"
# sampling frequency
sampling_freq <- 1000
# model to use
func_name <- "rosenzweig_macarthur_2sp_limit_cycle"
# Rosenzweig–MacArthur model with limit cycle
if (func_name == "rosenzweig_macarthur_2sp_limit_cycle") {
  time_type <- "continuous"
  n_sp <- 2
  func <- rosenzweig_macarthur_2sp_ode
  parms <- c(A = 0, p = 0.1, r = 5, k_mean = 1.8, a = 1.3, b = 1, e = 0.7, d = 0.2)
  state <- c(x1 = 0.51113819, x2 = 4.905690)
  n_points <- 500000
  time_step <- 0.001
  times <- seq(0, time_step * n_points, by = time_step)
}
# Hastings-Powell model with chaotic attractor
if (func_name == "hastings_powell_3sp_chaos") {
  time_type <- "continuous"
  n_sp <- 3
  func <- hastings_powell_3sp_ode
  parms <- c(r = 0, k = 0.99, a1 = 0.8036, a2_base = 0.23008, 
             e1 = 1, e2 = 1, b1 = 0.16129, b2 = 0.5, d1 = 0.4, d2 = 0.08)
  state <- c(x1 = 0.2445642, x2 = 0.2500377, x3 = 0.8477652)
  n_points <- 500000
  time_step <- 0.001
  times <- seq(0, time_step * n_points, by = time_step)
}

# generate time series ------------------------------
# species labels
sp_names <- paste("x", 1:n_sp, sep = "")
# generate abundance trajectory
ts <- as.data.frame(ode(y = state, times = times, func = func, parms = parms, method = "ode45"))
ts <- ts[seq(1, nrow(ts), by = sampling_freq), ]
plot(ts$time, ts$x1, type = "l", ylim = c(min(ts[ , -1]), max(ts[ , -1])), col = "red")
lines(ts$time, ts$x2, type = "l", ylim = c(min(ts[ , -1]), max(ts[ , -1])), col = "blue")
if (n_sp == 3) {
  lines(ts$time, ts$x3, type = "l", ylim = c(min(ts[ , -1]), max(ts[ , -1])), col = "green")
}
# select species time series
xts <- ts[ , species]

# compute power spectrum (code adapted from Rogers et al 2023 Ecol Lett) ------------------------------
# settings
lam <- 0.01
scale <- TRUE
trim <- TRUE
# standardize data
if(scale == TRUE) {
  xtss <- (xts - mean(xts, na.rm = TRUE)) / sd(xts, na.rm = TRUE)
} else {
  xtss <- xts - mean(xts, na.rm = TRUE)
}
# get spectrum
if (trim) {  # ts length, exclude NAs on ends
  tsl <- length(which(!is.na(xtss))[1]:which(!is.na(xtss))[length(which(!is.na(xtss)))])
} else { # do not trim NAs on ends
  tsl <- length(xtss)
}
tsle <- floor(tsl/2) * 2 #ts length rounded down to even number (if odd)
fi <- (1:(tsle/2)) / tsle # frequencies (cycles per timestep)
per <- 1 / fi # periods (timesteps per cycle)
wi <- 2 * pi * fi # frequencies in radians per timestep
times <- 1:length(xtss)
cosbf <- cos(outer(times, wi))
sinbf <- sin(outer(times,wi))
allbf <- cbind(cosbf,sinbf) # all basis functions
y <- xtss[complete.cases(xtss)] # remove missing timepoints
X <- allbf[complete.cases(xtss), ] # remove missing timepoints
coefsdir <- solve(t(X) %*% X + lam*diag(ncol(X))) %*% t(X) %*% y
lmr <- sqrt(coefsdir[1:(length(coefsdir)/2), ]^2 + coefsdir[(length(coefsdir)/2+1):length(coefsdir), ]^2)
results_df <- data.frame(frequency = fi, period = per, power = lmr)
results_df$frequency[which.max(results_df$power)] # frequency of maximum power
results_df$period[which.max(results_df$power)] # period of maximum power
# plot power as a function of period
plot(results_df$frequency, results_df$power, type = "l")

# compute power spectrum using spectrum function from stats package ------------------------------
mspect <- spectrum(xts, log = "no", spans = c(2, 2), plot = FALSE)
specx <- mspect$freq
specy <- 2 * mspect$spec
specx[which.max(specy)] # frequency of maximum power
(1 / specx[which.max(specy)]) # period of maximum power
plot(specx, specy, xlab = "Frequency", ylab = "Spectral Density", type = "l")
