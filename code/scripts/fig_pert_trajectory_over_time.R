# Code that plots perturbed abundances over time using a given population
# dynamics model

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
func_name <- "rosenzweig_macarthur_2sp_limit_cycle"
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
pert_magn <- 0.1
# percentage of recurrence time to evolve perturbations
perc_rec_time <- 0.8
tau <- round(time_step * n_points * perc_rec_time, 1)
# define axis range
if (func_name == "hastings_powell_3sp_forced_cycle") {
  lim_x1 <- c(0.4, 0.85)
  lim_x2 <- c(0.1, 0.55)
  lim_x3 <- c(0.7, 1.15)
}
if (func_name == "hastings_powell_3sp_chaos") {
  lim_x1 <- c(0.2, 0.9)
  lim_x2 <- c(0, 0.7)
  lim_x3 <- c(0.45, 1.15)
}

# generate non-perturbed abundance trajectory ------------------------------
# generate time series
ts <- as.data.frame(ode(y = state, times = times, func = func, parms = parms, method = "ode45"))
n_points_keep <- floor(nrow(ts) * ((n_recurrences - 1) / n_recurrences))
# subset of points to store results
ts_sub <- head(ts, n_points_keep)[seq(1, n_points_keep, length = sample_size), ]
# value of forcing parameter over time
if (func_name == "rosenzweig_macarthur_2sp_forced_fixed_point") {
  k_sub <- head(k, n_points_keep)[seq(1, n_points_keep, length = sample_size)]
}
if (func_name == "hastings_powell_3sp_forced_cycle") {
  a2_sub <- head(a2, n_points_keep)[seq(1, n_points_keep, length = sample_size)]
}

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
  n_pert <- 20
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
indeces_j <- 1:nrow(pert_df)
ts_sub_sub <- ts_sub[time_ids, ]
ts_sub_sub$time_plot <- as.factor(round(ts_sub_sub$time, 0))

# apply perturbations and compute perturbation growth rate for a selected point ------------------------------
# to store results
sim_results_df <- data.frame()
# which point to use
i <- 1
# evolve perturbations over time
curr_state <- as.numeric(ts_sub[time_ids[i], sp_names])
names(curr_state) <- sp_names
time_step_pert <- time_step
times_pert <- seq(ts_sub$time[time_ids[i]], ts_sub$time[time_ids[i]] + tau, 
                  by = time_step_pert)
ts_unpert <- as.data.frame(ode(y = curr_state, times = times_pert, func = func, parms = parms, method = "ode45"))
names(ts_unpert) <- c("time", sp_names)
sim_results_df <- data.frame()
for (j in indeces_j) { 
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
  ts_pert$time <- seq(0, tau, by = time_step_pert)
  # add initial perturbation size
  ts_pert$pert_size_init <- sqrt(sum((curr_state - ts_sub[time_ids[i], sp_names])^2))
  # concatenate results
  sim_results_df <- rbind(sim_results_df, ts_pert)
}
# plot perturbations over time
sim_results_df$Perturbation <- as.factor(sim_results_df$perturbation)
if (n_sp == 2) {
  fig <- ggplot() +
    geom_path(data = sim_results_df, aes(x = x1, y = x2, color = Perturbation), size = 0.8) +
    geom_path(data = ts, aes(x = x1, y = x2), size = 0.6, color = "gray70") +
    geom_point(data = subset(sim_results_df, time == 0), aes(x = x1, y = x2), fill = "gray100", size = 1.4, shape = 21) +
    geom_point(data = subset(sim_results_df, time == tau), aes(x = x1, y = x2), fill = "gray40", size = 1.4, shape = 21) +
    xlab(latex2exp::TeX("Resource abundance ($N_1$)")) +
    ylab(latex2exp::TeX("Consumer abundance ($N_2$)")) +
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
    ggsave(paste("figs/", func_name, "_trajectory_perturbations.pdf", sep = ""), 
           fig, width = 24, height = 18, units = "cm")
  }
}
if (n_sp == 3) {
  fig <- plot_ly(x = ~x1, y = ~x2, z = ~x3, showlegend = FALSE) %>% 
    add_trace(data = sim_results_df, type = 'scatter3d', mode = "lines",
              color = ~Perturbation, line = list(width = 6)) %>% 
    add_trace(data = ts, type = 'scatter3d', mode = "lines",
              color = I('gray70'), line = list(width = 4)) %>% 
    add_markers(data = subset(sim_results_df, time == 0), type = 'scatter3d', 
                marker = list(size = 3, line = list(color = "black", width = 1)),
                color = I('gray100')) %>% 
    add_markers(data = subset(sim_results_df, time == tau), type = 'scatter3d', 
                marker = list(size = 3, line = list(color = "black", width = 1)),
                color = I('gray40')) %>% 
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
    orca(fig, paste("figs/", func_name, "_trajectory_perturbations.pdf", sep = ""), format = "pdf",
         width = 800, height = 600)
  }
}
