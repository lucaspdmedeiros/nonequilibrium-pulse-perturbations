# Similar to the code for Figure 2, but using different structures for the 
# covariance matrix of perturbations (Sigma)

# cleaning environment and loading functions ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/rosenzweig_macarthur_2sp_ode.R")
source("code/functions/hypersphere_sampling.R")

# settings ------------------------------
# to reproduce results 
set.seed(42)
# whether to save plots
save_plots <- TRUE
# scenario to use (this analysis is performed only for this scenario)
func_name <- "rosenzweig_macarthur_2sp_limit_cycle"
# load scenario settings
source("code/scripts/settings.R")
# species labels
sp_names <- paste("x", 1:n_sp, sep = "")
# species labels
sp_names <- paste("x", 1:n_sp, sep = "")
# limits for plot axis
lim_x1 <- c(0, 1.2)
lim_x2 <- c(2.7, 5.5)
# perturbation magnitude
pert_magn <- 0.15
# magnitude of covariance in Sigma (number between -1 and 1)
cov_magn <- 0
# which covariance matrix to use to compute analytical stability metrics (true or identity)
cov_mat <- "true"

# generate non-perturbed abundance trajectory ------------------------------
# generate time series
ts <- as.data.frame(ode(y = state, times = times, func = func, parms = parms, method = "ode45"))
n_points_keep <- floor(nrow(ts) * ((n_recurrences - 1) / n_recurrences))
# subset of points to store results
ts_sub <- head(ts, n_points_keep)[floor(seq(1, n_points_keep, length = sample_size)), ]
# times to store analytical results
select_times <- ts_sub$time
select_times[1] <- times_pert[11]
# rescale time
ts_sub$time_plot <- (ts_sub$time / max(ts_sub$time)) * 100
# plot
fig <- ggplot() +
  geom_path(data = ts, aes(x = x1, y = x2), size = 1, color = "gray70") +
  geom_point(data = ts_sub, aes(x = x1, y = x2), size = 5) +
  xlab(latex2exp::TeX("Resource abundance ($N_1$)")) +
  ylab(latex2exp::TeX("Consumer abundance ($N_2$)")) +
  scale_x_continuous(limits = lim_x1) +
  scale_y_continuous(limits = lim_x2) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.size = unit(0.5, "cm"))

# compute Jacobian matrix along trajectory ------------------------------
J <- dlply(ts, "time", function(x) jacobian.full(y = unlist(c(x[2:(n_sp + 1)])), 
                                                 func = func,
                                                 parms = parms))
names(J) <- sprintf("%.5f", ts$time)

# sample perturbation directions from multivariate normal distribution ------------------------------
# variance and covariance will be proportional to standard deviation of abundances
pert_size_init <- pert_magn * mean(apply(ts[ , sp_names], 2, sd))
# number of perturbations
n_pert <- 50
# mean vector
mu <- rep(0, n_sp)
# covariance matrix
Sigma <- matrix(c(1, cov_magn, cov_magn, 1), nrow = n_sp, ncol = n_sp, byrow = TRUE)
Sigma <- (pert_size_init^2) * Sigma
# sample n_pert perturbations
pert_list <- replicate(n_pert, rmvnorm(n = 1, mean = mu, sigma = Sigma), simplify = FALSE)
# transform into data frame
pert_df <- data.frame(matrix(unlist(pert_list), nrow = length(pert_list), byrow = TRUE))
names(pert_df) <- sp_names

# generate and plot time series of perturbed abundances ------------------------------
sim_results_df <- data.frame()
i <- 1
indeces_j <- 1:nrow(pert_df)
curr_state <- as.numeric(ts_sub[i, sp_names])
names(curr_state) <- sp_names
ts_unpert <- as.data.frame(ode(y = curr_state, times = times_pert, func = func, parms = parms, method = "ode45"))
names(ts_unpert) <- c("time", sp_names)
sim_results_df <- data.frame()
for (j in indeces_j) { 
  curr_state <- as.numeric(ts_sub[i, sp_names] + pert_df[j, ])
  names(curr_state) <- sp_names
  curr_state[curr_state <= 0] <- 0
  ts_pert <- as.data.frame(ode(y = curr_state, times = times_pert, func = func, parms = parms, method = "ode45"))
  names(ts_pert) <- c("time", sp_names)
  ts_pert$perturbation <- j
  ts_pert$time_ref <- ts_sub[i, "time"]
  ts_pert$x1_ref <- ts_unpert$x1
  ts_pert$x2_ref <- ts_unpert$x2
  # rescale time
  ts_pert$time_plot <- (ts_pert$time / max(ts_pert$time)) * 100
  # add initial perturbation size
  ts_pert$pert_size_init <- sqrt(sum((curr_state - ts_sub[i, sp_names])^2))
  # concatenate results
  sim_results_df <- rbind(sim_results_df, ts_pert)
}
# plot perturbation trajectories in state space
plot_df_pert <- sim_results_df[!is.na(match(sim_results_df$time, ts_sub$time)), c("time", "time_plot", "x1", "x2", "perturbation")]
plot_df_pert$perturbation <- as.factor(plot_df_pert$perturbation)
plot_df_unpert <- ts_sub
plot_df_unpert$perturbation <- 0
plot_df_unpert$perturbation <- as.factor(plot_df_unpert$perturbation)
sim_results_df$Perturbation <- as.factor(sim_results_df$perturbation)
fig <- ggplot() +
  geom_path(data = ts, aes(x = x1, y = x2), size = 1, color = "gray70") +
  geom_path(data = subset(sim_results_df, time < 8), aes(x = x1, y = x2, color = Perturbation), size = 1) +
  geom_point(data = subset(sim_results_df, time == 0), aes(x = x1, y = x2), fill = "gray70", size = 2, shape = 21) +
  xlab(latex2exp::TeX("Resource abundance ($N_1$)")) +
  ylab(latex2exp::TeX("Consumer abundance ($N_2$)")) +
  scale_x_continuous(limits = lim_x1) +
  scale_y_continuous(limits = lim_x2) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position = "none")
if (save_plots) {
  ggsave(paste("figs/", func_name, "_perturbation_trajectories_time_", ts_sub$time[i],
               "_multiple_gaussian_covariance_", cov_magn, ".pdf", sep = ""), 
         fig, width = 14, height = 14, units = "cm")
}

# compute perturbation growth rate over time from simulations ------------------------------
sim_results_df$size <- apply(sim_results_df, 1, function(x) sum((as.numeric(x[sp_names]) - 
                                                                   as.numeric(x[paste("x", 1:n_sp, "_ref", sep = "")]))^2))
sim_results_df$growth_rate <- (log(sim_results_df$size) - log(sim_results_df$pert_size_init^2)) / (2 * sim_results_df$time)
mean_sim_results_df <- ddply(sim_results_df, c("time_ref", "time", "time_plot"), summarise,
                             mean_size = mean(size),
                             mean_growth_rate = mean(growth_rate),
                             median_size = median(size),
                             median_growth_rate = median(growth_rate),
                             mean_pert_size_init = mean(pert_size_init))
mean_sim_results_df$growth_rate_of_mean_size <- (log(mean_sim_results_df$mean_size) - log(mean_sim_results_df$mean_pert_size_init^2)) / (2 * mean_sim_results_df$time)

# compute and plot analytical perturbation growth rate over time ------------------------------
if (cov_mat == "identity") {
  Sigma <- diag(1, n_sp)
}
mean_ana_results_df <- data.frame()
curr_time <- ts_sub$time[i]
curr_state <- as.numeric(ts_sub[i, sp_names])
names(curr_state) <- sp_names
max_eigen <- c()
min_singular <- c() 
max_singular <- c() 
min_growth_rate <- c() 
max_growth_rate <- c() 
expect_inst_growth_rate <- c()
expect_avg_size <- c()
expect_avg_growth_rate <- c()
for (j in 1:length(select_times)) {
  J_curr <- J[sprintf("%.5f", seq(curr_time, (curr_time + select_times[j]), by = time_step))]
  Phi <- Reduce("%*%", rev(lapply(J_curr, function(A) expm(time_step * A))))
  max_eigen[j] <- max(Re(eigen(Phi)$values))
  min_singular[j] <- min(Re(eigen(Phi %*% t(Phi))$values))
  max_singular[j] <- max(Re(eigen(Phi %*% t(Phi))$values))
  min_growth_rate[j] <- log(min_singular[j]) / (2 * select_times[j])
  max_growth_rate[j] <- log(max_singular[j]) / (2 * select_times[j])
  expect_inst_growth_rate[j] <- log(det(Phi)) / (n_sp * select_times[j])
  expect_avg_size[j] <- sqrt(sum(diag(Phi %*% Sigma %*% t(Phi))))
  expect_avg_growth_rate[j] <- (log(sum(diag(Phi %*% Sigma %*% t(Phi)))) - log(sum(diag(Sigma)))) / (2 * select_times[j])
}
J_curr <- J[names(J) == sprintf("%.5f", curr_time)][[1]]
max_eigen_continuous <- max(Re(eigen(J_curr)$values))
max_growth_rate_continuous <- max(Re(eigen((J_curr + t(J_curr)) / 2)$values))
expect_inst_growth_rate_continuous <- 2 * sum(diag(J_curr)) / n_sp
curr_mean_ana_results_df <- data.frame(time = select_times, 
                                       time_ref = ts_sub$time[i],
                                       max_eigen = max_eigen,
                                       min_singular = min_singular,
                                       max_singular = max_singular,
                                       min_growth_rate = min_growth_rate,
                                       max_growth_rate = max_growth_rate,
                                       expect_inst_growth_rate = expect_inst_growth_rate,
                                       expect_avg_size = expect_avg_size,
                                       expect_avg_growth_rate = expect_avg_growth_rate,
                                       max_eigen_continuous = max_eigen_continuous,
                                       max_growth_rate_continuous = max_growth_rate_continuous,
                                       expect_inst_growth_rate_continuous = expect_inst_growth_rate_continuous)
mean_ana_results_df <- rbind(mean_ana_results_df, curr_mean_ana_results_df)
# rescale time
mean_ana_results_df$time_plot <- (mean_ana_results_df$time / max(mean_ana_results_df$time)) * 100
# plot perturbation growth rate over time
sim_results_df$perturbation <- as.factor(sim_results_df$perturbation)
fig <- ggplot() +
  geom_hline(yintercept = 0, size = 0.7) +
  geom_line(data = subset(sim_results_df, time >= select_times[1] & time_ref == ts_sub$time[i]), 
            aes(x = time, y = growth_rate, group = perturbation, color = perturbation), size = 0.7) +
  geom_line(data = subset(mean_sim_results_df, time >= select_times[1] & time_ref == ts_sub$time[i]), 
            aes(x = time, y = median_growth_rate), size = 1.5) +
  geom_point(data = mean_ana_results_df, aes(x = time, y = expect_avg_growth_rate),
             size = 3.5, shape = 22, fill = "gray70") +
  geom_point(data = mean_ana_results_df, aes(x = time, y = max_growth_rate),
             size = 3.5, shape = 23, fill = "gray70") +
  geom_point(data = mean_ana_results_df, aes(x = time, y = min_growth_rate),
             size = 3.5, shape = 24, fill = "gray70") +
  scale_x_continuous(limits = c(min(sim_results_df$time) - time_step * 10, max(sim_results_df$time) + time_step * 10)) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  xlab(latex2exp::TeX(r"(Time ($\tau$))")) +
  ylab(latex2exp::TeX(r"(Perturbation growth rate ($r_{\tau}$))")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position = "none",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.size = unit(0.5, "cm"))
if (save_plots) {
  ggsave(paste("figs/", func_name, "_perturbation_growth_rate_time_", ts_sub$time[i], 
               "_multiple_gaussian_covariance_", cov_magn, "_cov_matrix_", cov_mat, ".pdf", sep = ""), 
         fig, width = 12, height = 6, units = "cm")
}