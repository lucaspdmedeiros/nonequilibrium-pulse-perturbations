# Code for Figure 3

# cleaning environment and loading functions ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/rosenzweig_macarthur_2sp_ode.R")
source("code/functions/hastings_powell_3sp_ode.R")
source("code/functions/hypersphere_sampling.R")

# settings ------------------------------
# to reproduce results 
set.seed(42)
# whether to save plots
save_plots <- TRUE
# model to use
func_name <- "rosenzweig_macarthur_2sp_forced_fixed_point"
# load model settings
source("code/scripts/settings.R")
# species labels
sp_names <- paste("x", 1:n_sp, sep = "")
# how to sample perturbations (equidistant or uniform)
pert_sampling <- "uniform"
# times to store analytical results
select_times <- times_pert[floor(seq(11, length(times_pert), length = sample_size))]
# color palettes
palette_1 <- c(brewer.pal(9, "Blues")[8], brewer.pal(9, "YlOrRd")[4], 
               brewer.pal(9, "Reds")[8])
palette_2 <- c(brewer.pal(9, "Blues")[5], brewer.pal(9, "YlOrRd")[2], 
               brewer.pal(9, "Reds")[5])
# number of points to use for analysis
sample_size <- 50 * 2^n_sp
# perturbation magnitude
pert_magn <- 0.05

# generate non-perturbed abundance trajectory ------------------------------
# generate time series
ts <- as.data.frame(ode(y = state, times = times, func = func, parms = parms, method = "ode45"))
n_points_keep <- floor(nrow(ts) * ((n_recurrences - 1) / n_recurrences))
# subset of points to store results
ts_sub <- head(ts, n_points_keep)[floor(seq(1, n_points_keep, length = sample_size)), ]
# value of forcing parameter over time
if (func_name == "rosenzweig_macarthur_2sp_forced_fixed_point") {
  k_sub <- head(k, n_points_keep)[floor(seq(1, n_points_keep, length = sample_size))]
}
if (func_name == "hastings_powell_3sp_forced_cycle") {
  a2_sub <- head(a2, n_points_keep)[floor(seq(1, n_points_keep, length = sample_size))]
}

# compute Jacobian matrix along trajectory ------------------------------
J <- dlply(ts, "time", function(x) jacobian.full(y = unlist(c(x[2:(n_sp + 1)])), 
                                                 func = func,
                                                 parms = parms))
if (func_name == "rosenzweig_macarthur_2sp_forced_fixed_point") {
  J <- list()
  x1 <- ts$x1
  x2 <- ts$x2
  for (i in 1:nrow(ts)) {
    J[[i]] <- matrix(c(r - ((2 * r * x1[i]) / k[i]) - ((a * x2[i] * (b + x1[i]) - a * x1[i] * x2[i]) / (b + x1[i])^2), 
                       - (a * x1[i]) / (b + x1[i]),
                       (e * a * x2[i] * (b + x1[i]) - e * a * x1[i] * x2[i]) / (b + x1[i])^2, 
                       ((e * a * x1[i]) / (b + x1[i])) - d), nrow = 2, byrow = TRUE)
  }
}
names(J) <- sprintf("%.5f", ts$time)

# sample perturbation directions ------------------------------
pert_size_init <- pert_magn * mean(apply(ts[ , sp_names], 2, sd))
if (pert_sampling == "equidistant" & n_sp == 2) {
  n_pert <- 20
  angles <- seq(0.1, 2*pi + 0.1, length = n_pert)
  sin_angles <- sin(angles) * pert_size_init
  cos_angles <- cos(angles) * pert_size_init
  pert_df <- data.frame(x1 = cos_angles, x2 = sin_angles)
}
if (pert_sampling == "uniform") {
  n_pert <- 500
  pert_list <- replicate(n_pert, hypersphere_sampling(n_sp, positive = FALSE, within = FALSE), simplify = FALSE)
  pert_df <- matrix(unlist(pert_list), nrow = length(pert_list), byrow = TRUE)
  pert_df <- data.frame(t(t(pert_df) * pert_size_init))
  names(pert_df) <- sp_names
}

# determine points to perform perturbation analysis ------------------------------
if (func_name == "rosenzweig_macarthur_2sp_forced_fixed_point") {
  id_1 <- which.min(apply(ts_sub[ , -1], 1, function(x) sum((x - c(0.11, 2.99))^2)))
  id_2 <- which.min(apply(ts_sub[ , -1], 1, function(x) sum((x - c(0.19, 4.43))^2)))
  id_3 <- which.min(apply(ts_sub[ , -1], 1, function(x) sum((x - c(0.67, 3.45))^2)))
}
if (func_name == "rosenzweig_macarthur_2sp_limit_cycle") {
  id_1 <- which.min(apply(ts_sub[ , -1], 1, function(x) sum((x - c(0.15, 3.08))^2)))
  id_2 <- which.min(apply(ts_sub[ , -1], 1, function(x) sum((x - c(0.51, 4.90))^2)))
  id_3 <- which.min(apply(ts_sub[ , -1], 1, function(x) sum((x - c(0.82, 3.79))^2)))
}
if (func_name == "hastings_powell_3sp_forced_cycle") {
  id_1 <- which.min(apply(ts_sub[ , -1], 1, function(x) sum((x - c(0.45, 0.42, 0.97))^2)))
  id_2 <- which.min(apply(ts_sub[ , -1], 1, function(x) sum((x - c(0.65, 0.36, 0.90))^2)))
  id_3 <- which.min(apply(ts_sub[ , -1], 1, function(x) sum((x - c(0.80, 0.21, 0.92))^2)))
}
if (func_name == "hastings_powell_3sp_chaos") {
  id_1 <- which.min(apply(ts_sub[ , -1], 1, function(x) sum((x - c(0.20, 0.22, 0.81))^2)))
  id_2 <- which.min(apply(ts_sub[ , -1], 1, function(x) sum((x - c(0.57, 0.44, 0.74))^2)))
  id_3 <- which.min(apply(ts_sub[ , -1], 1, function(x) sum((x - c(0.88, 0.13, 0.83))^2)))
}
time_ids <- c(id_1, id_2, id_3)
ts_sub_sub <- ts_sub[time_ids, ]
ts_sub_sub$time_plot <- as.factor(round(ts_sub_sub$time, 0))
# plot abundance trajectories with chosen points
if (n_sp == 2) {
  fig <- ggplot() +
    geom_path(data = ts, aes(x = x1, y = x2), size = 1.5, color = "gray70") +
    geom_point(data = ts_sub_sub, aes(x = x1, y = x2), 
               size = 15, shape = 21, fill = palette_2) +
    xlab(latex2exp::TeX("Resource abundance ($N_1$)")) +
    ylab(latex2exp::TeX("Consumer abundance ($N_2$)")) +
    scale_x_continuous(limits = lim_x1) +
    scale_y_continuous(limits = lim_x2) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 18),
          axis.title = element_text(size = 24),
          axis.text.x = element_text(size = 18),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          legend.title = element_blank(),
          legend.text = element_text(size = 16),
          legend.key.size = unit(0.5, "cm"))
  if (save_plots) {
    ggsave(paste("figs/", func_name, "_sample_points_trajectory.pdf", sep = ""), 
           fig, width = 12, height = 12, units = "cm")
  }
}
if (n_sp == 3) {
  if (func_name == "hastings_powell_3sp_chaos") {
    palette_ordered <- c(palette_2[2], palette_2[1], palette_2[3])
  }
  if (func_name == "hastings_powell_3sp_forced_cycle") {
    palette_ordered <- c(palette_2[2], palette_2[3], palette_2[1])
  }
  fig <- plot_ly(x = ~x1, y = ~x2, z = ~x3, showlegend = FALSE,
                 colors = palette_ordered) %>% 
    add_trace(data = ts, type = 'scatter3d', mode = "lines",
              color = I('gray70'), line = list(width = 7)) %>% 
    add_markers(data = ts_sub_sub, type = 'scatter3d', marker = list(size = 26, line = list(color = "black", width = 1)),
                color = ~time_plot) %>% 
    layout(scene = list(xaxis = list(title = "",
                                     titlefont = list(size = 26, 
                                                      family = "Arial, sans-serif"), 
                                     tickfont = list(size = 16,
                                                     family = "Arial, sans-serif"),
                                     range = lim_x1,
                                     ticklen = 6,
                                     gridwidth = 1.2,
                                     zerolinewidth = 0,
                                     showgrid = FALSE, 
                                     showline = TRUE),
                        yaxis = list(title = "",
                                     titlefont = list(size = 26, 
                                                      family = "Arial, sans-serif"),
                                     tickfont = list(size = 16,
                                                     family = "Arial, sans-serif"),
                                     range = lim_x2,
                                     ticklen = 6,
                                     gridwidth = 1.2,
                                     zerolinewidth = 0,
                                     showgrid = FALSE, 
                                     showline = TRUE),
                        zaxis = list(title = "",
                                     titlefont = list(size = 26, 
                                                      family = "Arial, sans-serif"),
                                     tickfont = list(size = 16,
                                                     family = "Arial, sans-serif"),
                                     range = lim_x3,
                                     ticklen = 6,
                                     gridwidth = 1.2,
                                     zerolinewidth = 0,
                                     showgrid = FALSE, 
                                     showline = TRUE),
                        camera = list(eye = list(x = 1.8, y = 0.8, z = 1)),
                        aspectratio = list(x = 0.8, y = 0.8, z = 0.8)),
           margin = list(l = 0,
                         r = 0,
                         b = 0,
                         t = 0))
  if (save_plots) {
    orca(fig, paste("figs/", func_name, "_sample_points_trajectory.pdf", sep = ""), format = "pdf",
         width = 800, height = 600)
  }
}

# apply perturbations and compute perturbation growth rate for selected points ------------------------------
sim_results_df <- data.frame()
for (i in 1:length(time_ids)) {
  curr_state <- as.numeric(ts_sub[time_ids[i], sp_names])
  names(curr_state) <- sp_names
  time_step_pert <- time_step
  times_pert <- seq(ts_sub$time[time_ids[i]], ts_sub$time[time_ids[i]] + time_step_pert * n_points, 
                    by = time_step_pert)
  ts_unpert <- as.data.frame(ode(y = curr_state, times = times_pert, func = func, parms = parms, method = "ode45"))
  names(ts_unpert) <- c("time", sp_names)
  sim_results_df <- data.frame()
  for (j in 1:nrow(pert_df)) { 
    print(j)
    curr_state <- as.numeric(ts_sub[time_ids[i], sp_names] + pert_df[j, ])
    names(curr_state) <- sp_names
    curr_state[curr_state <= 0] <- 0
    ts_pert <- as.data.frame(ode(y = curr_state, times = times_pert, func = func, parms = parms, method = "ode45"))
    names(ts_pert) <- c("time", sp_names)
    ts_pert$perturbation <- j
    ts_pert$time_ref <- ts_sub[time_ids[i], "time"]
    ts_pert$x1_ref <- ts_unpert$x1
    ts_pert$x2_ref <- ts_unpert$x2
    if (n_sp == 3) {
      ts_pert$x3_ref <- ts_unpert$x3
    }
    # update time
    ts_pert$time <- seq(0, time_step_pert * n_points, by = time_step_pert)
    # add initial perturbation size
    ts_pert$pert_size_init <- sqrt(sum((curr_state - ts_sub[time_ids[i], sp_names])^2))
    # concatenate results
    sim_results_df <- rbind(sim_results_df, ts_pert)
  }
  
  # compute perturbation growth rate over time from simulations ------------------------------
  sim_results_df$size <- apply(sim_results_df, 1, function(x) sum((as.numeric(x[sp_names]) - 
                                                                     as.numeric(x[paste("x", 1:n_sp, "_ref", sep = "")]))^2))
  sim_results_df$growth_rate <- (log(sim_results_df$size) - log(sim_results_df$pert_size_init^2)) / (2 * sim_results_df$time)
  mean_sim_results_df <- ddply(sim_results_df, c("time_ref", "time"), summarise,
                               mean_size = mean(size),
                               mean_growth_rate = mean(growth_rate),
                               median_size = median(size),
                               median_growth_rate = median(growth_rate),
                               min_growth_rate = min(growth_rate),
                               Q_min_growth_rate = quantile(growth_rate, probs = 0.25, na.rm = TRUE),
                               Q_max_growth_rate = quantile(growth_rate, probs = 0.75, na.rm = TRUE),
                               max_growth_rate = max(growth_rate),
                               mean_pert_size_init = mean(pert_size_init))
  mean_sim_results_df$growth_rate_of_mean_size <- (log(mean_sim_results_df$mean_size) - log(mean_sim_results_df$mean_pert_size_init^2)) / (2 * mean_sim_results_df$time)
  
  # compute perturbation growth rate analytically ------------------------------
  mean_ana_results_df <- data.frame()
  Sigma <- diag(apply(pert_df, 2, var))
  curr_time <- ts_sub$time[time_ids[i]]
  curr_state <- as.numeric(ts_sub[time_ids[i], sp_names])
  names(curr_state) <- sp_names
  min_growth_rate_eigen <- c()
  max_growth_rate_eigen <- c()
  min_growth_rate <- c() 
  max_growth_rate <- c() 
  expect_inst_growth_rate <- c()
  expect_avg_size <- c()
  expect_avg_growth_rate <- c()
  for (j in 1:length(select_times)) {
    J_curr <- J[sprintf("%.5f", seq(curr_time, (curr_time + select_times[j]), by = time_step))]
    Phi <- Reduce("%*%", rev(lapply(J_curr, function(A) expm(time_step * A))))
    singular_phi <- svd(Phi)$d
    min_growth_rate[j] <- log(min(singular_phi)) / select_times[j]
    max_growth_rate[j] <- log(max(singular_phi)) / select_times[j]
    eigen_phi <- eigen(Phi)$values
    abs_eigen_phi <- sqrt(Re(eigen_phi)^2 + Im(eigen_phi)^2)
    min_growth_rate_eigen[j] <- log(min(abs_eigen_phi)) / select_times[j]
    max_growth_rate_eigen[j] <- log(max(abs_eigen_phi)) / select_times[j]
    expect_inst_growth_rate[j] <- log(det(Phi)) / (n_sp * select_times[j])
    expect_avg_size[j] <- sqrt(sum(diag(Phi %*% Sigma %*% t(Phi))))
    expect_avg_growth_rate[j] <- (log(sum(diag(Phi %*% Sigma %*% t(Phi)))) - log(sum(diag(Sigma)))) / (2 * select_times[j])
  }
  min_growth_rate[!is.finite(min_growth_rate)] <- min(min_growth_rate[is.finite(min_growth_rate)])
  J_curr <- J[names(J) == sprintf("%.5f", curr_time)][[1]]
  curr_mean_ana_results_df <- data.frame(time = select_times, 
                                         time_ref = ts_sub$time[time_ids[i]],
                                         min_growth_rate_eigen = min_growth_rate_eigen,
                                         max_growth_rate_eigen = max_growth_rate_eigen,
                                         min_growth_rate = min_growth_rate,
                                         max_growth_rate = max_growth_rate,
                                         expect_inst_growth_rate = expect_inst_growth_rate,
                                         expect_avg_size = expect_avg_size,
                                         expect_avg_growth_rate = expect_avg_growth_rate)
  mean_ana_results_df <- rbind(mean_ana_results_df, curr_mean_ana_results_df)
  # rescale time
  mean_ana_results_df$time_plot <- (mean_ana_results_df$time / max(mean_ana_results_df$time)) * 100
  # plot perturbation growth rate over time
  fig <- ggplot() +
    geom_hline(yintercept = 0, size = 0.7) +
    geom_ribbon(data = subset(mean_sim_results_df, time >= select_times[1]),
                aes(x = time, ymin = Q_min_growth_rate, ymax = Q_max_growth_rate),
                fill = palette_1[i], alpha = 0.8) + 
    geom_ribbon(data = subset(mean_sim_results_df, time >= select_times[1]),
                aes(x = time, ymin = min_growth_rate, ymax = Q_min_growth_rate),
                fill = palette_2[i], alpha = 0.8) + 
    geom_ribbon(data = subset(mean_sim_results_df, time >= select_times[1]),
                aes(x = time, ymin = Q_max_growth_rate, ymax = max_growth_rate),
                fill = palette_2[i], alpha = 0.8) + 
    geom_line(data = subset(mean_sim_results_df, time >= select_times[1]), 
              aes(x = time, y = median_growth_rate), size = 1.5) +
    geom_point(data = mean_ana_results_df, aes(x = time, y = expect_avg_growth_rate),
               size = 4, shape = 22, fill = "gray70") +
    geom_point(data = mean_ana_results_df, aes(x = time, y = max_growth_rate),
               size = 4, shape = 23, fill = "gray70") +
    geom_point(data = mean_ana_results_df, aes(x = time, y = min_growth_rate),
               size = 4, shape = 24, fill = "gray70") +
    scale_x_continuous(limits = c(min(sim_results_df$time) - time_step * 10, max(sim_results_df$time) + time_step * 10)) +
    scale_y_continuous(labels = label_number(accuracy = 0.1)) +
    xlab(latex2exp::TeX(r"(Time ($\tau$))")) +
    ylab(latex2exp::TeX(r"(Perturbation growth rate ($r_{\tau}$))")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.5),
          axis.text.y = element_text(size = 18),
          axis.title = element_text(size = 24),
          axis.text.x = element_text(size = 18),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          legend.position = "none",
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16),
          legend.key.size = unit(0.5, "cm"))
  if (save_plots) {
    ggsave(paste("figs/", func_name, "_perturbation_growth_rate_time_", ts_sub$time[time_ids[i]], 
                 ".pdf", sep = ""), 
           fig, width = 17, height = 9, units = "cm")
  }
}
