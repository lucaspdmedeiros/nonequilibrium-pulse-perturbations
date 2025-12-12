# Code for Figure 4

# cleaning environment and loading functions ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/rosenzweig_macarthur_2sp_ode.R")
source("code/functions/hastings_powell_3sp_ode.R")

# settings ------------------------------
# to reproduce results 
set.seed(42)
# whether to save plots
save_plots <- TRUE
# there are 4 possible scenarios to use for func_name: 
# 1) rosenzweig_macarthur_2sp_forced_fixed_point
# 2) rosenzweig_macarthur_2sp_limit_cycle 
# 3) hastings_powell_3sp_forced_cycle
# 4) hastings_powell_3sp_chaos
func_name <- "rosenzweig_macarthur_2sp_forced_fixed_point"
# load scenario settings
source("code/scripts/settings.R")
# species labels
sp_names <- paste("x", 1:n_sp, sep = "")
# number of points to use for analysis
sample_size <- 50 * 2^n_sp
# color palettes
palette <- colorRampPalette(brewer.pal(11, "RdYlBu"))(sample_size)
# percentage of recurrence time to evolve perturbations (set to 0.2 or 0.8 to produce other plots in Fig 4)
perc_rec_time <- 0.02
tau <- time_step * n_points * perc_rec_time
# which stability metric to use: 
# - expect_avg_growth_rate for Fig 4 
# - max_growth_rate for Fig S9 
# - max_growth_rate_eigen for Fig S10
# - expect_inst_growth_rate for Fig S11
metric <- "expect_avg_growth_rate"

# generate non-perturbed abundance trajectory ------------------------------
# generate time series
ts <- as.data.frame(ode(y = state, times = times, func = func, parms = parms, method = "ode45"))
n_points_keep <- floor(nrow(ts) * ((n_recurrences - 1) / n_recurrences))
# subset of points to store results
ts_sub <- head(ts, n_points_keep)[floor(seq(1, n_points_keep, length = sample_size)), ]
# time series to show time evolution
if (func_name == "rosenzweig_macarthur_2sp_forced_fixed_point") {
  ts_example <- ts[which(round(ts$time, 3) == round(0, 3)):which(round(ts$time, 3) == round(tau, 3)), ]
}
if (func_name == "rosenzweig_macarthur_2sp_limit_cycle") {
  ts_example <- ts[which(round(ts$time, 3) == round(5.865, 3)):which(round(ts$time, 3) == round(5.865+tau, 3)), ]
}
if (func_name == "hastings_powell_3sp_forced_cycle") {
  ts_example <- ts[which(round(ts$time, 3) == round(351, 3)):which(round(ts$time, 3) == round(351+tau, 3)), ]
  ts_example$x3 <- ts_example$x3 + 0.02
}
if (func_name == "hastings_powell_3sp_chaos") {
  ts_example <- ts[which(round(ts$time, 3) == round(369.474, 3)):which(round(ts$time, 3) == round(369.474+tau, 3)), ]
  ts_example$x3 <- ts_example$x3 + 0.015
}

# compute Jacobian matrix along trajectory ------------------------------
if ((func_name == "rosenzweig_macarthur_2sp_forced_fixed_point") | 
    (func_name == "rosenzweig_macarthur_2sp_limit_cycle")) {
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
if ((func_name == "hastings_powell_3sp_forced_cycle") |
    (func_name == "hastings_powell_3sp_chaos")) {
  J <- list()
  x1 <- ts$x1
  x2 <- ts$x2
  x3 <- ts$x3
  for (i in 1:nrow(ts)) {
    J[[i]] <- matrix(c((1 - (x1[i] / k)) - x1[i] * (1 / k) - (a1 * x2[i] / (x1[i] + b1) - (a1 * x1[i] * x2[i]) / (x1[i] + b1)^2),
                       -(a1 * x1[i] / (x1[i] + b1)),
                       0,
                       a1 * e1 * x2[i] / (x1[i] + b1) - (a1 * e1 * x1[i] * x2[i]) / (x1[i] + b1)^2,
                       a1 * e1 * x1[i] / (x1[i] + b1) - d1 - (a2[i] * x3[i] / (x2[i] + b2) - (a2[i] * x2[i] * x3[i]) / (x2[i] + b2)^2),
                       -(a2[i] * x2[i] / (x2[i] + b2)),
                       0,
                       a2[i] * e2 * x3[i] / (x2[i] + b2) - (a2[i] * e2 * x2[i] * x3[i]) / (x2[i] + b2)^2,
                       a2[i] * e2 * x2[i] / (x2[i] + b2) - d2), nrow = 3, byrow = TRUE)
  }
}
names(J) <- sprintf("%.5f", ts$time)

# compute and plot perturbation growth rate along trajectory ------------------------------
mean_ana_results_df <- data.frame()
Sigma <- diag(rep(1, n_sp))
min_growth_rate <- c() 
max_growth_rate <- c() 
min_growth_rate_eigen <- c()
max_growth_rate_eigen <- c()
expect_avg_growth_rate <- c()
expect_inst_growth_rate <- c()
for (j in 1:nrow(ts_sub)) {
  print(j)
  curr_time <- ts_sub$time[j]
  curr_state <- as.numeric(ts_sub[j, sp_names])
  names(curr_state) <- sp_names
  J_curr <- J[sprintf("%.5f", seq(curr_time, (curr_time + tau), by = time_step))]
  Phi <- Reduce("%*%", rev(lapply(J_curr, function(A) expm(time_step * A))))
  singular_squared_phi <- eigen(Phi %*% t(Phi))$values
  min_growth_rate[j] <- log(min(singular_squared_phi)) / (2 * tau)
  max_growth_rate[j] <- log(max(singular_squared_phi)) / (2 * tau)
  eigen_phi <- eigen(Phi)$values
  abs_eigen_phi <- sqrt(Re(eigen_phi)^2 + Im(eigen_phi)^2)
  min_growth_rate_eigen[j] <- log(min(abs_eigen_phi)) / tau
  max_growth_rate_eigen[j] <- log(max(abs_eigen_phi)) / tau
  expect_avg_growth_rate[j] <- (log(sum(diag(Phi %*% Sigma %*% t(Phi)))) - log(sum(diag(Sigma)))) / (2 * tau)
  expect_inst_growth_rate[j] <- log(det(Phi)) / (n_sp * tau)
}
plot_df <- ts_sub
if (metric == "expect_avg_growth_rate") {
  plot_df$metric <- expect_avg_growth_rate
  if (func_name == "rosenzweig_macarthur_2sp_forced_fixed_point") {
    axis_edge <- 0.405
    legend_breaks <- c(-0.4, -0.2, 0, 0.2, 0.4)
  }
  if (func_name == "rosenzweig_macarthur_2sp_limit_cycle") {
    axis_edge <- 1.02
    legend_breaks <- c(-1, -0.5, 0, 0.5, 1)
  }
  if (func_name == "hastings_powell_3sp_forced_cycle") {
    axis_edge <- 0.065
    legend_breaks <- c(-0.06, -0.03, 0, 0.03, 0.06)
  }
  if (func_name == "hastings_powell_3sp_chaos") {
    axis_edge <- 0.13
    legend_breaks <- c(-0.12, -0.06, 0, 0.06, 0.12)
  }
}
if (metric == "max_growth_rate") {
  plot_df$metric <- max_growth_rate
  if (func_name == "rosenzweig_macarthur_2sp_forced_fixed_point") {
    axis_edge <- 1.07
    legend_breaks <- c(-1, -0.5, 0, 0.5, 1)
  }
  if (func_name == "rosenzweig_macarthur_2sp_limit_cycle") {
    axis_edge <- 1.71
    legend_breaks <- c(-1.5, -0.75, 0, 0.75, 1.5)
  }
  if (func_name == "hastings_powell_3sp_forced_cycle") {
    axis_edge <- 0.255
    legend_breaks <- c(-0.2, -0.1, 0, 0.1, 0.2)
  }
  if (func_name == "hastings_powell_3sp_chaos") {
    axis_edge <- 0.4
    legend_breaks <- c(-0.4, -0.2, 0, 0.2, 0.4)
  }
}
if (metric == "max_growth_rate_eigen") {
  plot_df$metric <- max_growth_rate_eigen
  if (func_name == "rosenzweig_macarthur_2sp_forced_fixed_point") {
    axis_edge <- 1.07
    legend_breaks <- c(-1, -0.5, 0, 0.5, 1)
  }
  if (func_name == "rosenzweig_macarthur_2sp_limit_cycle") {
    axis_edge <- 1.71
    legend_breaks <- c(-1.5, -0.75, 0, 0.75, 1.5)
  }
  if (func_name == "hastings_powell_3sp_forced_cycle") {
    axis_edge <- 0.255
    legend_breaks <- c(-0.2, -0.1, 0, 0.1, 0.2)
  }
  if (func_name == "hastings_powell_3sp_chaos") {
    axis_edge <- 0.4
    legend_breaks <- c(-0.4, -0.2, 0, 0.2, 0.4)
  }
}
if (metric == "expect_inst_growth_rate") {
  plot_df$metric <- expect_inst_growth_rate
  if (func_name == "rosenzweig_macarthur_2sp_forced_fixed_point") {
    axis_edge <- 0.51
    legend_breaks <- c(-0.5, -0.25, 0, 0.25, 0.5)
  }
  if (func_name == "rosenzweig_macarthur_2sp_limit_cycle") {
    axis_edge <- 1.02
    legend_breaks <- c(-1, -0.5, 0, 0.5, 1)
  }
  if (func_name == "hastings_powell_3sp_forced_cycle") {
    axis_edge <- 0.21
    legend_breaks <- c(-0.2, -0.1, 0, 0.1, 0.2)
  }
  if (func_name == "hastings_powell_3sp_chaos") {
    axis_edge <- 0.27
    legend_breaks <- c(-0.2, -0.1, 0, 0.1, 0.2)
  }
}
# plot abundance trajectory colored by stability metric (for 2 species)
if (n_sp == 2) {
  fig <- ggplot() +
    geom_point(data = plot_df, aes(x = x1, y = x2, color = metric), 
               size = 5) +
    geom_path(data = ts_example, aes(x = x1, y = x2), size = 1, color = "black", 
              arrow = arrow(ends = "last", length = unit(0.2, "inches"))) +
    scale_color_gradientn(colors = brewer.pal(11, "RdYlBu"), 
                         limits = c(-axis_edge, axis_edge)) +
    xlab(latex2exp::TeX("Resource abundance ($N_1$)")) +
    ylab(latex2exp::TeX("Consumer abundance ($N_2$)")) +
    scale_x_continuous(limits = lim_x1) +
    scale_y_continuous(limits = lim_x2) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 20),
          axis.text.x = element_text(size = 14),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          legend.position = "none")
  if (save_plots) {
    ggsave(paste("figs/", func_name, "_trajectory_color_pert_growth_rate_tau_", perc_rec_time,  ".pdf", sep = ""), 
           fig, width = 11.5, height = 11.5, units = "cm")
  }
  # plot to extract color scale
  if (perc_rec_time == 0.02) {
    fig <- ggplot() +
      geom_point(data = plot_df, aes(x = x1, y = x2, color = metric), 
                 size = 5) +
      geom_path(data = ts_example, aes(x = x1, y = x2), size = 1, color = "black", 
                arrow = arrow(ends = "last", length = unit(0.2, "inches"))) +
      scale_color_gradientn(colors = brewer.pal(11, "RdYlBu"), 
                            limits = c(-axis_edge, axis_edge), breaks = legend_breaks) +
      xlab(latex2exp::TeX("Resource abundance")) +
      ylab(latex2exp::TeX("Consumer abundance")) +
      scale_x_continuous(limits = lim_x1) +
      scale_y_continuous(limits = lim_x2) +
      theme_classic() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = 14),
            axis.title = element_text(size = 20),
            axis.text.x = element_text(size = 14),
            plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
            legend.position = "top",
            legend.title = element_blank(),
            legend.text = element_text(size = 17),
            legend.key.width = unit(1.5, "cm"),
            legend.key.height = unit(0.8, "cm"))
    if (save_plots) {
      ggsave(paste("figs/", func_name, "_legend_tau_", perc_rec_time, ".pdf", sep = ""), 
             fig, width = 11.5, height = 11.5, units = "cm")
    }
  }
}
# plot abundance trajectory colored by stability metric (for 3 species)
if (n_sp == 3) {
  fig <- plot_ly(x = ~x1, y = ~x2, z = ~x3, colors = palette, showlegend = FALSE) %>% 
    add_markers(data = plot_df, type = 'scatter3d', marker = list(size = 8),
                color = ~metric) %>% 
    add_trace(data = ts_example, type = 'scatter3d', mode = "lines",
              color = I('black'), line = list(width = 8)) %>% 
    colorbar(limits = c(-axis_edge, axis_edge), title = "Average\nperturbation\ngrowth\nrate", orientation = "v") %>% 
    hide_colorbar() %>%
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
    orca(fig, paste("figs/", func_name, "_trajectory_color_pert_growth_rate_tau_", perc_rec_time,  ".pdf", sep = ""), 
         format = "pdf", width = 800, height = 600)
  }
  # plot to extract color scale
  if (perc_rec_time == 0.02) { 
    fig <- ggplot() +
      geom_point(data = plot_df, aes(x = x1, y = x2, color = metric), 
                 size = 5) +
      scale_color_gradientn(colors = brewer.pal(11, "RdYlBu"), 
                            limits = c(-axis_edge, axis_edge), breaks = legend_breaks) +
      xlab(latex2exp::TeX("Resource abundance")) +
      ylab(latex2exp::TeX("Consumer abundance")) +
      scale_x_continuous(limits = lim_x1) +
      scale_y_continuous(limits = lim_x2) +
      theme_classic() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = 14),
            axis.title = element_text(size = 20),
            axis.text.x = element_text(size = 14),
            plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
            legend.position = "top",
            legend.title = element_blank(),
            legend.text = element_text(size = 17),
            legend.key.width = unit(1.5, "cm"),
            legend.key.height = unit(0.8, "cm"))
    if (save_plots) {
      ggsave(paste("figs/", func_name, "_legend_tau_", perc_rec_time, ".pdf", sep = ""), 
             fig, width = 11.5, height = 11.5, units = "cm")
    }
  }
}
