# Code for Figure 5

# cleaning environment and loading functions ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/rosenzweig_macarthur_2sp_ode.R")
source("code/functions/power_spectrum.R")

# settings ------------------------------
# to reproduce results 
set.seed(42)
# whether to save plots
save_plots <- TRUE
# model to use
func_name <- "rosenzweig_macarthur_2sp_limit_cycle"
# load model settings
source("code/scripts/settings.R")
# whether to run analysis or just load saved results
run_analysis <- FALSE
# modify parameter and simulation settings
time_type <- "continuous"
n_sp <- 2
func <- rosenzweig_macarthur_2sp_ode
A <- 0
p <- 0.1
a <- 1.3
b <- 1
e <- 0.7
d <- 0.2
k_hopf <- (-b * (d + a * e) / (d - a * e))
state <- c(x1 = 0.5, x2 = 4.0)
lim_r <- c(1, 21)
lim_k <- c(0.25, 2.95)
# species labels
sp_names <- paste("x", 1:n_sp, sep = "")
# number of points to use for analysis
sample_size <- 200
# color palette
palette <- colorRampPalette(brewer.pal(11, "RdYlBu"))(sample_size)
# set of parameters to perform analysis
k_list <- seq(0.4, 2.8, by = 0.2)
r_list <- seq(2, 20, by = 1.5)
param_list <- expand.grid(k_list, r_list)
# covariance matrix to use
Sigma <- diag(rep(1, n_sp))

# compute equilibrium and non-equilibrium metrics for a set of parameter values ------------------------------
if (run_analysis) {
  # to store stability metrics for equilibrium point
  stability_eq_df <- data.frame()
  # to store stability metrics for non-equilibrium trajectory
  stability_noneq_df <- data.frame()
  # loop over parameter values
  for (i in 1:nrow(param_list)) {
    print(i)
    # current parameters and equilibrium point
    k_mean <- param_list[i, 1]
    r <- param_list[i, 2]
    parms <- c(A = A, p = p, r = r, k_mean = k_mean, a = a, b = b, e = e, d = d)
    equilibrium <- c(x1 = - (b * d) / (d - a * e), x2 = - (b * e * (b * d + d * k_mean - a * e * k_mean) * r) / 
                       (k_mean * (d - a * e)^2))
    if (any(equilibrium < 0)) {
      print("Unfeasible equilibrium")
    }
    # compute stability metrics at the equilibrium point
    curr_state <- as.numeric(equilibrium)
    names(curr_state) <- sp_names
    # define time to evolve perturbations
    tau_eq <- 1.5
    # compute Jacobian and stability metrics
    J <- jacobian.full(y = curr_state, func = func, parms = parms)
    Phi <- expm(tau_eq * J)
    max_eigen <- max(Re(eigen(J)$values))
    min_singular <- min(Re(eigen(Phi %*% t(Phi))$values))
    max_singular <- max(Re(eigen(Phi %*% t(Phi))$values))
    min_growth_rate <- log(min_singular) / (2 * tau_eq)
    max_growth_rate <- log(max_singular) / (2 * tau_eq)
    expect_avg_growth_rate <- (log(sum(diag(Phi %*% Sigma %*% t(Phi)))) - log(sum(diag(Sigma)))) / (2 * tau_eq)
    # save results
    curr_stability_eq_df <- data.frame(k = k_mean,
                                       r = r,
                                       x1 = as.numeric(equilibrium[1]),
                                       x2 = as.numeric(equilibrium[2]),
                                       max_eigen = max_eigen,
                                       min_growth_rate = min_growth_rate,
                                       max_growth_rate = max_growth_rate,
                                       expect_avg_growth_rate = expect_avg_growth_rate)
    stability_eq_df <- rbind(stability_eq_df, curr_stability_eq_df)
    # compute stability metrics for the non-equilibrium trajectory
    if (k_mean > k_hopf) {
      # define time steps
      n_points <- 500000
      time_step <- 0.001
      times <- seq(0, time_step * n_points, by = time_step)
      # generate time series
      ts <- as.data.frame(ode(y = state, times = times, func = func, parms = parms, method = "ode45"))
      ts <- ts[seq(1, nrow(ts), by = 1000), ]
      # dominant period
      dom_period <- round(power_spectrum(ts[ , c(1, 2)], lam = 0.01, scale = TRUE, trim = TRUE, plot = FALSE)[[2]][2], 1)
      # redefine time points
      n_recurrences <- 2
      n_points <- 10000
      time_step <- dom_period / n_points
      times <- seq(0, time_step * n_points * n_recurrences, by = time_step)
      # extract last point
      curr_state <- as.numeric(ts[nrow(ts), -1])
      names(curr_state) <- sp_names
      # generate time series with initial condition on limit cycle
      ts <- as.data.frame(ode(y = curr_state, times = times, func = func, parms = parms, method = "ode45"))
      ts$time <- ts$time - ts$time[1]
      # subset of points to store results
      n_points_keep <- floor(nrow(ts) * ((n_recurrences - 1) / n_recurrences))
      ts_sub <- head(ts, n_points_keep)[floor(seq(1, n_points_keep, length = sample_size)), ]
      # compute Jacobian matrix along trajectory
      J <- dlply(ts, "time", function(x) jacobian.full(y = unlist(c(x[2:(n_sp + 1)])), 
                                                       func = func,
                                                       parms = parms))
      # rename Jacobian matrices
      names(J) <- sprintf("%.5f", ts$time)
      tau <- ts$time[which.min(abs(ts$time - tau_eq))]
      # compute stability metrics along trajectory
      max_eigen <- c()
      min_growth_rate <- c() 
      max_growth_rate <- c() 
      expect_avg_growth_rate <- c()
      for (j in 1:nrow(ts_sub)) {
        curr_time <- ts_sub$time[j]
        curr_state <- as.numeric(ts_sub[j, sp_names])
        names(curr_state) <- sp_names
        J_curr <- J[sprintf("%.5f", seq(curr_time, (curr_time + tau), by = time_step))]
        Phi <- Reduce("%*%", rev(lapply(J_curr, function(A) expm(time_step * A))))
        max_eigen[j] <- max(Re(eigen(Phi)$values))
        min_singular <- min(Re(eigen(Phi %*% t(Phi))$values))
        max_singular <- max(Re(eigen(Phi %*% t(Phi))$values))
        min_growth_rate[j] <- log(min_singular) / (2 * tau)
        max_growth_rate[j] <- log(max_singular) / (2 * tau)
        expect_avg_growth_rate[j] <- (log(sum(diag(Phi %*% Sigma %*% t(Phi)))) - log(sum(diag(Sigma)))) / (2 * tau)
      }
      # save results
      curr_stability_noneq_df <- data.frame(k = k_mean,
                                            r = r,
                                            x1 = ts_sub$x1,
                                            x2 = ts_sub$x2,
                                            max_eigen = max_eigen,
                                            min_growth_rate = min_growth_rate,
                                            max_growth_rate = max_growth_rate,
                                            expect_avg_growth_rate = expect_avg_growth_rate)
      stability_noneq_df <- rbind(stability_noneq_df, curr_stability_noneq_df, rep(NA, ncol(curr_stability_noneq_df)))
    }
  }
  # save results
  write.csv(x = stability_eq_df, file = "tables/stability_eq_df.csv", row.names = FALSE)
  write.csv(x = stability_noneq_df, file = "tables/stability_noneq_df.csv", row.names = FALSE)
} else {
  # load results from previous analysis
  stability_eq_df <- read.csv(file = "tables/stability_eq_df.csv")
  stability_noneq_df <- read.csv(file = "tables/stability_noneq_df.csv")
}

# plot median perturbation growth rate for each parameter combination  ------------------------------
# compute maximum median perturbation growth rate along nonequilibrium trajectory
summ_stability_noneq_df <- ddply(stability_noneq_df, c("k", "r"), summarise,
                                 max_expect_avg_growth_rate = max(expect_avg_growth_rate))
# define axis limits for plots
axis_edge_eq <- max(abs(c(min(subset(stability_eq_df, k < k_hopf)$expect_avg_growth_rate, na.rm = TRUE), 
                          max(subset(stability_eq_df, k < k_hopf)$expect_avg_growth_rate, na.rm = TRUE))))
axis_edge_noneq <- max(abs(c(min(subset(summ_stability_noneq_df, k > k_hopf)$max_expect_avg_growth_rate, na.rm = TRUE), 
                             max(subset(summ_stability_noneq_df, k > k_hopf)$max_expect_avg_growth_rate, na.rm = TRUE))))
# plot for equilibrium case
fig <- ggplot(data = subset(stability_eq_df, k < k_hopf), aes(x = k, y = r, fill = expect_avg_growth_rate)) +
  geom_tile() +
  scale_fill_gradientn(colors = brewer.pal(11, "RdYlBu"), 
                       limits = c(-axis_edge_eq, axis_edge_eq)) +
  geom_vline(xintercept = k_hopf, size = 1, linetype = "dashed") +
  scale_x_continuous(limits = lim_k, expand = c(0, 0)) +
  scale_y_continuous(limits = lim_r, expand = c(0, 0)) +
  xlab("K") +
  ylab("r") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.width = unit(1.2, "cm"),
        legend.key.height = unit(0.8, "cm"))
if (save_plots) {
  ggsave(paste("figs/", func_name, "_median_growth_rate_eq_many_k_and_r", ".pdf", sep = ""), 
         fig, width = 15, height = 15, units = "cm")
}
# plot for nonequilibrium case
fig <- ggplot(data = subset(summ_stability_noneq_df, k > k_hopf), aes(x = k, y = r, fill = max_expect_avg_growth_rate)) +
  geom_tile() +
  scale_fill_gradientn(colors = brewer.pal(11, "RdYlBu"), 
                        limits = c(-axis_edge_noneq, axis_edge_noneq)) +
  geom_vline(xintercept = k_hopf, size = 1, linetype = "dashed") +
  scale_x_continuous(limits = lim_k, expand = c(0, 0)) +
  scale_y_continuous(limits = lim_r, expand = c(0, 0)) +
  xlab("K") +
  ylab("r") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.width = unit(1.2, "cm"),
        legend.key.height = unit(0.8, "cm"))
if (save_plots) {
  ggsave(paste("figs/", func_name, "_median_growth_rate_noneq_many_k_and_r", ".pdf", sep = ""), 
         fig, width = 15, height = 15, units = "cm")
}

# plot nonequilibrium trajectories colored by median perturbation growth rate ------------------------------
# unstable equilibrium point for k = 1.8
equilibrium_df <- data.frame(x1 = - (b * d) / (d - a * e), 
                             x2 = - (b * e * (b * d + d * 1.8 - a * e * 1.8) * 5) / 
                               (1.8 * (d - a * e)^2))
# plot for k = 1.8
fig <- ggplot() +
  geom_point(data = equilibrium_df, aes(x = x1, y = x2), fill = "gray70", size = 5, shape = 21) +
  geom_point(data = subset(stability_noneq_df, round(k, 1) == 1.8 & round(r, 1) == 5.0), 
             aes(x = x1, y = x2, color = expect_avg_growth_rate), size = 5) +
  scale_color_gradientn(colors = brewer.pal(11, "RdYlBu"), 
                        limits = c(-axis_edge_noneq, axis_edge_noneq)) +
  xlab(latex2exp::TeX("Resource abundance")) +
  ylab(latex2exp::TeX("Consumer abundance")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 17),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 17),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position = "none")
if (save_plots) {
  ggsave(paste("figs/", func_name, "_trajectory_median_growth_rate_noneq_k_1.8", ".pdf", sep = ""), 
         fig, width = 9, height = 9, units = "cm")
}
# unstable equilibrium point for k = 2.6
equilibrium_df <- data.frame(x1 = - (b * d) / (d - a * e), 
                             x2 = - (b * e * (b * d + d * 2.6 - a * e * 2.6) * 17) / 
                               (2.6 * (d - a * e)^2))
# plot for k = 2.6
fig <- ggplot() +
  geom_point(data = equilibrium_df, aes(x = x1, y = x2), fill = "gray70", size = 5, shape = 21) +
  geom_point(data = subset(stability_noneq_df, round(k, 1) == 2.6 & round(r, 1) == 17.0), 
             aes(x = x1, y = x2, color = expect_avg_growth_rate), size = 5) +
  scale_color_gradientn(colors = brewer.pal(11, "RdYlBu"), 
                        limits = c(-axis_edge_noneq, axis_edge_noneq)) +
  xlab(latex2exp::TeX("Resource abundance")) +
  ylab(latex2exp::TeX("Consumer abundance")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 17),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 17),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position = "none")
if (save_plots) {
  ggsave(paste("figs/", func_name, "_trajectory_median_growth_rate_noneq_k_2.6", ".pdf", sep = ""), 
         fig, width = 9, height = 9, units = "cm")
}
