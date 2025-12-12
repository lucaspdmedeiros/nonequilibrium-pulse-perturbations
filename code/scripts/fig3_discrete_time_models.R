# Similar to the code for Figure 3, but using discrete-time models

# cleaning environment and loading functions ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/larvae_pupae_adults_3sp_map.R")
source("code/functions/predator_prey_2sp_map.R")
source("code/functions/hypersphere_sampling.R")

# settings ------------------------------
# to reproduce results 
set.seed(42)
# whether to save plots
save_plots <- TRUE
# there are 2 possible scenarios to use for func_name: 
# 1) predator_prey_2sp_limit_cycle
# 2) larvae_pupae_adults_3sp_chaos 
func_name <- "predator_prey_2sp_limit_cycle"
# load scenario settings
source("code/scripts/settings.R")
# species labels
sp_names <- paste("x", 1:n_sp, sep = "")
# how to sample perturbations (equidistant or uniform)
pert_sampling <- "uniform"
# color palettes
palette_1 <- c(brewer.pal(9, "Blues")[8], brewer.pal(9, "YlOrRd")[4], 
               brewer.pal(9, "Reds")[8])
palette_2 <- c(brewer.pal(9, "Blues")[5], brewer.pal(9, "YlOrRd")[2], 
               brewer.pal(9, "Reds")[5])
# perturbation magnitude
pert_magn <- 0.05

# generate non-perturbed abundance trajectory ------------------------------
# generate time series
ts <- func(t_max = n_points * n_recurrences, x = state, p = parms)
n_points_keep <- floor(nrow(ts) * ((n_recurrences - 1) / n_recurrences))
# subset of points to store results
ts_sub <- head(ts, n_points_keep)
# times to store analytical results
select_times <- 0:(n_points_keep - 1)

# compute Jacobian matrix along trajectory ------------------------------
J <- list()
if (func_name == "predator_prey_2sp_limit_cycle") {
  x1 <- ts$x1
  x2 <- ts$x2 
  exp_f_1 <- expression(x1 + tau * x1 * (a - x1 - (b * x2) / ((1 + alpha * x1) * (1 + beta * x2))))
  exp_j_11 <- D(exp_f_1, 'x1')
  j_11 <- eval(exp_j_11)
  exp_j_12 <- D(exp_f_1, 'x2')
  j_12 <- eval(exp_j_12)
  exp_f_2 <- expression(x2 + tau * x2 * (-c + (d * x1) / ((1 + alpha * x1) * (1 + beta * x2))))
  exp_j_21 <- D(exp_f_2, 'x1')
  j_21 <- eval(exp_j_21)
  exp_j_22 <- D(exp_f_2, 'x2')
  j_22 <- eval(exp_j_22)
  J_list <- list()
  for (t in 1:nrow(ts)) {
    J[[t]] <- matrix(c(j_11[t], j_12[t], j_21[t], j_22[t]), 
                     nrow = n_sp, ncol = n_sp, byrow = TRUE)
  }
}
if (func_name == "larvae_pupae_adults_3sp_chaos") {
  x1 <- ts$x1
  x2 <- ts$x2
  x3 <- ts$x3
  for (i in 1:nrow(ts)) {
    J[[i]] <- matrix(c(b * x3[i] * (-c_el) * exp(-c_el * x1[i] - c_ea * x3[i]), 0,
                       b * exp(-c_el * x1[i] - c_ea * x3[i]) + b * x3[i] * (-c_ea) * exp(-c_el * x1[i] - c_ea * x3[i]),
                       (1 - u_l), 0, 0, 0, exp(-c_pa * x3[i]),
                       x2[i] * (-c_pa) * exp(-c_pa * x3[i]) + (1 - u_a)), nrow = 3, byrow = TRUE)
  }
}
names(J) <- ts$time

# sample perturbation directions ------------------------------
pert_size_init <- pert_magn * mean(apply(ts[ , sp_names], 2, sd))
if (pert_sampling == "uniform") {
  n_pert <- 500
  pert_list <- replicate(n_pert, hypersphere_sampling(n_sp, positive = FALSE, within = FALSE), simplify = FALSE)
  pert_df <- matrix(unlist(pert_list), nrow = length(pert_list), byrow = TRUE)
  pert_df <- data.frame(t(t(pert_df) * pert_size_init))
  names(pert_df) <- sp_names
}

# determine points to perform perturbation analysis ------------------------------
if (func_name == "predator_prey_2sp_limit_cycle") {
  id_1 <- 2
  id_2 <- 1
  id_3 <- 5
}
if (func_name == "larvae_pupae_adults_3sp_chaos") {
  id_1 <- 3
  id_2 <- 1
  id_3 <- 2
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
  palette_ordered <- c(palette_2[2], palette_2[1], palette_2[3])
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
  ts_unpert <- func(t_max = max(select_times) + 2, x = curr_state, p = parms)
  names(ts_unpert) <- c("time", sp_names)
  sim_results_df <- data.frame()
  for (j in 1:nrow(pert_df)) { 
    print(j)
    curr_state <- as.numeric(ts_sub[time_ids[i], sp_names] + pert_df[j, ])
    names(curr_state) <- sp_names
    curr_state[curr_state <= 0] <- 0
    ts_pert <- func(t_max = max(select_times) + 2, x = curr_state, p = parms)
    names(ts_pert) <- c("time", sp_names)
    ts_pert$perturbation <- j
    ts_pert$time_ref <- ts_sub[time_ids[i], "time"]
    ts_pert$x1_ref <- ts_unpert$x1
    ts_pert$x2_ref <- ts_unpert$x2
    if (n_sp == 3) {
      ts_pert$x3_ref <- ts_unpert$x3
    }
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
    J_curr <- J[as.character(seq(curr_time, (curr_time + select_times[j]), by = 1))]
    Phi <- Reduce("%*%", rev(lapply(J_curr, function(A) A)))
    singular_phi <- svd(Phi)$d
    min_growth_rate[j] <- log(min(singular_phi)) / (select_times[j] + 1)
    max_growth_rate[j] <- log(max(singular_phi)) / (select_times[j] + 1)
    eigen_phi <- eigen(Phi)$values
    abs_eigen_phi <- sqrt(Re(eigen_phi)^2 + Im(eigen_phi)^2)
    min_growth_rate_eigen[j] <- log(min(abs_eigen_phi)) / (select_times[j] + 1)
    max_growth_rate_eigen[j] <- log(max(abs_eigen_phi)) / (select_times[j] + 1)
    expect_inst_growth_rate[j] <- log(det(Phi)) / (n_sp * (select_times[j] + 1))
    expect_avg_size[j] <- sqrt(sum(diag(Phi %*% Sigma %*% t(Phi))))
    expect_avg_growth_rate[j] <- (log(sum(diag(Phi %*% Sigma %*% t(Phi)))) - log(sum(diag(Sigma)))) / (2 * (select_times[j] + 1))
  }
  curr_mean_ana_results_df <- data.frame(time = select_times + 1, 
                                         time_ref = ts_sub$time[time_ids[i]],
                                         min_growth_rate_eigen = min_growth_rate_eigen,
                                         max_growth_rate_eigen = max_growth_rate_eigen,
                                         min_growth_rate = min_growth_rate,
                                         max_growth_rate = max_growth_rate,
                                         expect_inst_growth_rate = expect_inst_growth_rate,
                                         expect_avg_size = expect_avg_size,
                                         expect_avg_growth_rate = expect_avg_growth_rate)
  mean_ana_results_df <- rbind(mean_ana_results_df, curr_mean_ana_results_df)
  # plot perturbation growth rate over time
  if (func_name == "predator_prey_2sp_limit_cycle") {
    y_scale <- scale_y_continuous(labels = label_number(accuracy = 0.1))
  }
  if (func_name == "larvae_pupae_adults_3sp_chaos") {
    y_scale <- scale_y_continuous(limits = c(-4, NA), labels = label_number(accuracy = 0.1)) 
  }
  fig <- ggplot() +
    geom_hline(yintercept = 0, size = 0.7) +
    geom_ribbon(data = subset(mean_sim_results_df, time > select_times[1]),
                aes(x = time, ymin = Q_min_growth_rate, ymax = Q_max_growth_rate),
                fill = palette_1[i], alpha = 0.8) + 
    geom_ribbon(data = subset(mean_sim_results_df, time > select_times[1]),
                aes(x = time, ymin = min_growth_rate, ymax = Q_min_growth_rate),
                fill = palette_2[i], alpha = 0.8) + 
    geom_ribbon(data = subset(mean_sim_results_df, time > select_times[1]),
                aes(x = time, ymin = Q_max_growth_rate, ymax = max_growth_rate),
                fill = palette_2[i], alpha = 0.8) + 
    geom_line(data = subset(mean_sim_results_df, time > select_times[1]), 
              aes(x = time, y = median_growth_rate), size = 1.5) +
    geom_point(data = mean_ana_results_df, aes(x = time, y = expect_avg_growth_rate),
               size = 4, shape = 22, fill = "gray70") +
    geom_point(data = mean_ana_results_df, aes(x = time, y = max_growth_rate),
               size = 4, shape = 23, fill = "gray70") +
    geom_point(data = mean_ana_results_df, aes(x = time, y = min_growth_rate),
               size = 4, shape = 24, fill = "gray70") +
    y_scale +
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
