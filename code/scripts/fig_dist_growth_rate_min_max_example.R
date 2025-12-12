# Code that plots the distribution of perturbation growth rates (similar to Fig 3),
# but highlights the perturbations that give the maximum and minimum growth rates

# cleaning environment and loading functions ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/rosenzweig_macarthur_2sp_ode.R")
source("code/functions/hastings_powell_3sp_ode.R")
source("code/functions/hypersphere_sampling.R")

# settings ------------------------------
# to reproduce results 
set.seed(42)
# scenario to use (this analysis is performed only for this scenario)
func_name <- "hastings_powell_3sp_chaos"
# load scenario settings
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
time_ids <- id_1
i <- 1
ts_sub_sub <- ts_sub[time_ids, ]
ts_sub_sub$time_plot <- as.factor(round(ts_sub_sub$time, 0))

# apply perturbations and compute perturbation growth rate for selected points ------------------------------
sim_results_df <- data.frame()
# add perturbation that gives the minimum growth rate
if ((func_name == "hastings_powell_3sp_chaos") & i == 1) {
  v <- c(c(-0.00552497841827141994, 
           0.00402508479412213866, 
           -0.00142123758506340582))
  pert_df <- rbind(pert_df, v)
}
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
# verifying which perturbations give the maximum and minimum growth rates at tau = 42
sub_sim_results_df <- subset(sim_results_df, time == tail(select_times, 1))
sub_sim_results_df[which.max(sub_sim_results_df$growth_rate), ]
sub_sim_results_df[which.min(sub_sim_results_df$growth_rate), ]
# removing perturbation that gives the minimum growth rate
sub_sim_results_df <- sub_sim_results_df[-which(sub_sim_results_df$perturbation == 501), ]

# plot distribution of perturbation growth rate at last tau value ------------------------------
dens <- density(sub_sim_results_df$growth_rate)
plot_df <- data.frame(x = dens$x, y = dens$y)
probs <- c(0.25, 0.5, 0.75)
quantiles <- quantile(sub_sim_results_df$growth_rate, prob = probs)
plot_df$quant <- factor(findInterval(plot_df$x, quantiles))
fig <- ggplot(data = plot_df, aes(x = x, y = y)) +
  geom_line(size = 1, color = "white") + 
  geom_ribbon(aes(ymin = 0, ymax = y, fill = quant)) + 
  geom_vline(xintercept = c(0.04120289), 
             size = 1.3, color = "#FF7F00") + 
  geom_vline(xintercept = c(-0.3025808), 
             size = 1.3, color = "#E41A1C") + 
  scale_x_continuous(limits = c(-0.31, 0.06)) +
  scale_fill_manual(values = c(palette_2[i], 
                               palette_1[i], palette_1[i],
                               palette_2[i])) +
  xlab(latex2exp::TeX(r"(Perturbation growth rate ($r_{\tau}$))")) +
  ylab("Density") +
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

# plot distribution of perturbation growth rate for all tau ------------------------------
# extracting perturbation that gives the maximum
max_sub_sim_results_df <- subset(sim_results_df, perturbation == 33)
pert_df[33, ]
# extracting perturbation that gives the minimum
min_sub_sim_results_df <- subset(sim_results_df, perturbation == 501)
pert_df[501, ]
# removing perturbation that gives the minimum growth rate
sim_results_df <- sim_results_df[-which(sim_results_df$perturbation == 501), ]
# compute median and average perturbation growth rate
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
# plot
fig <- ggplot() +
  geom_hline(yintercept = 0, size = 0.7) +
  geom_ribbon(data = subset(mean_sim_results_df, time >= select_times[1]),
              aes(x = time, ymin = Q_min_growth_rate, ymax = Q_max_growth_rate),
              fill = palette_1[i]) + 
  geom_ribbon(data = subset(mean_sim_results_df, time >= select_times[1]),
              aes(x = time, ymin = min_growth_rate, ymax = Q_min_growth_rate),
              fill = palette_2[i]) + 
  geom_ribbon(data = subset(mean_sim_results_df, time >= select_times[1]),
              aes(x = time, ymin = Q_max_growth_rate, ymax = max_growth_rate),
              fill = palette_2[i]) + 
  geom_line(data = subset(max_sub_sim_results_df, time >= select_times[1]), 
            aes(x = time, y = growth_rate), size = 1.3, color = "#FF7F00") +
  geom_line(data = subset(min_sub_sim_results_df, time >= select_times[1]), 
            aes(x = time, y = growth_rate), size = 1.3, color = "#E41A1C") +
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
