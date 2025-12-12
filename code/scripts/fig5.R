# Code for Figure 5

# cleaning environment and loading functions ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/rosenzweig_macarthur_2sp_ode.R")
source("code/functions/rosenzweig_macarthur_2sp_sde.R")
source("code/functions/hastings_powell_3sp_ode.R")
source("code/functions/hastings_powell_3sp_sde.R")
source("code/functions/regularized_smap_fit.R")
source("code/functions/regularized_smap_forecast.R")
source("code/functions/regularized_smap_cv.R")

# settings ------------------------------
# whether to save plots
save_plots <- TRUE
# there are 2 possible scenarios to use for func_name: 
# 1) rosenzweig_macarthur_2sp_limit_cycle_stochastic
# 2) hastings_powell_3sp_chaos_stochastic 
func_name <- "rosenzweig_macarthur_2sp_limit_cycle_stochastic"
# load scenario settings
source("code/scripts/settings.R")
# to reproduce results 
seed <- 10
set.seed(seed)
julia <- julia_setup()
julia_command("using Random")
julia_command(paste("Random.seed!(", seed, ")", sep = ""))
# species labels
sp_names <- paste("x", 1:n_sp, sep = "")
# number of points to use for analysis
sample_size <- 150
# amount of observational noise
obs_noise <- 0.1
# time step used in the S-map
delta_t <- 1
# whether we use abundance or growth rate as response variable in the S-map
response <- "continuous_growth_rate"

# generate non-perturbed abundance trajectory ------------------------------
# generate time series
ts <- func(state = state, p = parms, times = times)
# subset of points to store results
ts_sub <- ts[floor(seq(1, nrow(ts), length = sample_size)), ]

# compute analytical Jacobian matrix using ODE without process noise ------------------------------
if (func_name == "rosenzweig_macarthur_2sp_limit_cycle_stochastic") {
  # parameter values
  parms <- c(A = 0, p = 0, r = 5, k_mean = 1.8, a = 1.3, b = 1, e = 0.7, d = 0.2)
  # compute Jacobian
  J <- dlply(ts_sub, "time", function(x) jacobian.full(y = unlist(c(x[2:(n_sp + 1)])), 
                                                       func = rosenzweig_macarthur_2sp_ode,
                                                       parms = parms))
}
if (func_name == "hastings_powell_3sp_chaos_stochastic") {
  # parameter values
  parms <- c(r = 0, k = 0.99, a1 = 0.8036, a2_base = 0.23008, 
             e1 = 1, e2 = 1, b1 = 0.16129, b2 = 0.5, d1 = 0.4, d2 = 0.08)
  # compute Jacobian
  J <- dlply(ts_sub, "time", function(x) jacobian.full(y = unlist(c(x[2:(n_sp + 1)])), 
                                                       func = hastings_powell_3sp_ode,
                                                       parms = parms))
}
names(J) <- sprintf("%.5f", ts_sub$time)
# remove last Jacobian (cannot be inferred from time-series data)
J <- J[-length(J)]

# add observational noise to time series ------------------------------
ts_noise <- ts_sub
if (obs_noise > 0) {
  for (i in 1:n_sp) {
    for (j in 1:nrow(ts_noise)) {
      noisy_abund <- ts_noise[j, i+1] + rnorm(n = 1, mean = 0, sd = sd(ts_noise[ , i+1]) * obs_noise)
      while (noisy_abund < 0) {
        noisy_abund <- ts_noise[j, i+1] + rnorm(n = 1, mean = 0, sd = sd(ts_noise[ , i+1]) * obs_noise)
      }
      ts_noise[j, i+1] <- noisy_abund
    }
  }
}
ts_noise$time <- 0:(nrow(ts_noise)-1)

# plot abundance time series ------------------------------ 
# define color palette
if (n_sp == 2) {
  palette <- brewer.pal(9, "Set1")[c(3, 4)]
}
if (n_sp == 3) {
  palette <- brewer.pal(9, "Set1")[c(3, 4, 2)]
}
# data frame for plotting
plot_df <- gather(ts_noise, "species", "abundance", -time)
# plot species abundances through time
fig <- ggplot() +
  geom_line(data = plot_df, aes(x = time, y = abundance, color = species), 
            size = 0.5) +
  geom_point(data = plot_df, aes(x = time, y = abundance, fill = species), 
             size = 2.5, shape = 21) +
  geom_line(size = 0.4) +
  geom_point(size = 1.5) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  xlab("Time") +
  ylab("Species\nabundances") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        legend.position = "none")
if (save_plots) {
  ggsave(paste("figs/", func_name, "_time_series_abundances.pdf", sep = ""), 
         fig, width = 17, height = 6, units = "cm")
}

# plot attractor in state space ------------------------------ 
if (n_sp == 2) {
  fig <- ggplot() +
    geom_path(data = ts_noise, aes(x = x1, y = x2), size = 0.5, color = "gray70") +
    geom_point(data = ts_noise, aes(x = x1, y = x2), size = 3, fill = "gray70", shape = 21) +
    xlab(latex2exp::TeX("Resource abundance ($N_1$)")) +
    ylab(latex2exp::TeX("Consumer abundance ($N_2$)")) +
    scale_x_continuous(limits = lim_x1) +
    scale_y_continuous(limits = lim_x2) +
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
  if (save_plots) {
    ggsave(paste("figs/", func_name, "_trajectory_state_space.pdf", sep = ""), 
           fig, width = 9, height = 9, units = "cm")
  }
}
if (n_sp == 3) {
  fig <- plot_ly(x = ~x1, y = ~x2, z = ~x3, showlegend = FALSE) %>% 
    add_markers(data = ts_noise, type = 'scatter3d', color = I('gray70'),
                marker = list(size = 10)) %>% 
    add_trace(data = ts_noise, type = 'scatter3d', mode = "lines",
              color = I('gray70'), line = list(width = 2)) %>% 
    layout(scene = list(xaxis = list(title = "",
                                     titlefont = list(size = 28, 
                                                      family = "Arial, sans-serif"),
                                     tickfont = list(size = 18,
                                                     family = "Arial, sans-serif"),
                                     range = lim_x1,
                                     ticklen = 6,
                                     gridwidth = 1.2,
                                     zerolinewidth = 0,
                                     showgrid = FALSE, 
                                     showline = TRUE),
                        yaxis = list(title = "",
                                     titlefont = list(size = 28, 
                                                      family = "Arial, sans-serif"),
                                     tickfont = list(size = 18,
                                                     family = "Arial, sans-serif"),
                                     range = lim_x2,
                                     ticklen = 6,
                                     gridwidth = 1.2,
                                     zerolinewidth = 0,
                                     showgrid = FALSE, 
                                     showline = TRUE),
                        zaxis = list(title = "",
                                     titlefont = list(size = 28, 
                                                      family = "Arial, sans-serif"),
                                     tickfont = list(size = 18,
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
    orca(fig, paste("figs/", func_name, "_trajectory_state_space.pdf", sep = ""), 
         format = "pdf", width = 800, height = 600)
  }
}

# perform out-of-sample predictions to define S-map settings/hyperparameters ------------------------------
# size of training set
training_frac <- 0.5
# training and test sets
last_train_point <- floor(nrow(ts_noise) * training_frac)
training_ts <- ts_noise[1:last_train_point, ]
test_ts <- ts_noise[(last_train_point+1):nrow(ts_noise), ] 
# hyperparameters
theta <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 3, 4, 5, 6, 7, 8)
alpha <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
lambda <- c(0, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5)
df_hyperparameters <- expand_grid(theta, alpha, lambda)
# fit S-map and perform forecasts for a grid of hyperparameters
mean_R2 <- pbmapply(regularized_smap_cv, theta = df_hyperparameters$theta, 
                    alpha = df_hyperparameters$alpha, lambda = df_hyperparameters$lambda,
                    MoreArgs = list(training_ts = training_ts, test_ts = test_ts, response = response, 
                                    delta_t = delta_t, method = "glmnet", scale = TRUE, output = "mean_R2"))
df_hyperparameters$mean_R2 <- mean_R2
# order data frame according to R2
df_hyperparameters <- df_hyperparameters[order(df_hyperparameters$mean_R2, decreasing = TRUE), ]
# removing cases with negative R2
df_hyperparameters <- df_hyperparameters[df_hyperparameters$mean_R2 > 0, ]
# compute cumulative percent change in R2
R2_diff <- -c(diff(as.numeric(df_hyperparameters$mean_R2)), NA)
R2_perc_diff <- R2_diff / as.numeric(df_hyperparameters$mean_R2)
df_hyperparameters$R2_cum_perc_diff <- cumsum(R2_perc_diff)
# data frame of equivalently optimal R2 values (within 1% of highest R2)
df_hyperparameters_equivalent <- 
  df_hyperparameters[which(df_hyperparameters$R2_cum_perc_diff < 0.01), ]
# order new data frame according to theta, lambda, and alpha
df_hyperparameters_equivalent <- df_hyperparameters_equivalent[order(df_hyperparameters_equivalent$theta,
                                                                     df_hyperparameters_equivalent$lambda,
                                                                     df_hyperparameters_equivalent$alpha), ]
# select best hyperparameters
theta <- df_hyperparameters_equivalent$theta[1]
alpha <- df_hyperparameters_equivalent$alpha[1]
lambda <- df_hyperparameters_equivalent$lambda[1]

# fit S-map to infer Jacobian matrices over time with selected hyperparameters ------------------------------ 
# whether to transform estimated S-map coefficients
if (response == "abundance") {
  structure_transform <- "none"
}
if ((response == "continuous_growth_rate") | (response == "discrete_growth_rate")) {
  structure_transform <- "model_prediction"
}
# fit S-map
smap_results <- regularized_smap_fit(ts = ts_noise, points = 1:(nrow(ts_noise)-1), theta = theta, alpha = alpha, 
                                     lambda = lambda, response = response, delta_t = delta_t, method = "glmnet", 
                                     scale = TRUE, structure_transform = structure_transform)
smap_J <- smap_results[[1]]
# compute correlation between true and inferred Jacobian matrix at all times
cor_list <- c()
for (i in 1:(nrow(ts_noise)-1)) {
  cor_list[i] <- cor(c(J[[i]]), c(smap_J[[i]]), use = "complete.obs")
}

# plot true and inferred Jacobian elements through time ------------------------------ 
# create data frame from list of matrices
smap_J_df <- data.frame(matrix(unlist(smap_J), nrow = length(smap_J), byrow = TRUE))
names(smap_J_df) <- apply(expand.grid(1:n_sp, 1:n_sp), 1, 
                          function(x) paste("j", paste(x, collapse = ""), sep = "_"))
J_df <- data.frame(matrix(unlist(J), nrow = length(J), byrow = TRUE))
names(J_df) <- apply(expand.grid(1:n_sp, 1:n_sp), 1, 
                     function(x) paste("j", paste(x, collapse = ""), sep = "_"))
# add time
smap_J_df$time <- ts_noise$time[-nrow(ts_noise)]
smap_J_df$type <- "Inferred"
J_df$time <- ts_noise$time[-nrow(ts_noise)]
J_df$type <- "True"
# merge data frames to plot
both_J_df <- rbind(smap_J_df, J_df)
# data frame for plotting
plot_df <- gather(both_J_df, "coefficient", "value", -c(time, type))
plot_df$type <- factor(plot_df$type, levels = c("True", "Inferred"))
# plot coefficients through time
fig <- ggplot(data = plot_df, aes(x = time, y = value, color = type)) +
  geom_hline(yintercept = 0, size = 0.5, linetype = "dashed") +
  geom_line(size = 0.7) +
  facet_wrap(~coefficient) +
  scale_color_brewer(palette = "Set2") +
  xlab("Time") +
  ylab("Jacobian element value") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 0.5),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "white", size = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.size = unit(0.6, "cm"))

# compute stability metrics from Jacobian matrices ------------------------------ 
# determine time scale
tau <- 10
# data frames to store results
results_analytical_df <- data.frame(time = ts_noise$time, 
                                    max_growth_rate = rep(NA, nrow(ts_noise)),
                                    median_growth_rate = rep(NA, nrow(ts_noise)),
                                    type = "True")
results_smap_df <- data.frame(time = ts_noise$time, 
                              max_growth_rate = rep(NA, nrow(ts_noise)),
                              median_growth_rate = rep(NA, nrow(ts_noise)),
                              type = "Inferred")
# compute stability metrics at each point in time
for (t in 1:(nrow(ts_noise)-tau)) {
  # compute stability metrics using analytical Jacobian
  J_curr <- J[t:(t+tau-1)]
  Phi <- Reduce("%*%", rev(lapply(J_curr, function(A) expm(delta_t * A))))
  singular_phi <- svd(Phi)$d
  max_growth_rate <- log(max(singular_phi)) / tau
  median_growth_rate <- log(sum(diag(Phi %*% t(Phi))) / n_sp) / (2 * tau)
  results_analytical_df[t, "max_growth_rate"] <- max_growth_rate
  results_analytical_df[t, "median_growth_rate"] <- median_growth_rate
  # compute stability metrics using Jacobian inferred with S-map
  smap_J_curr <- smap_J[t:(t+tau-1)]
  if ((response == "abundance") | (response == "discrete_growth_rate")) {
    smap_Phi <- Reduce("%*%", rev(lapply(smap_J_curr, function(A) A)))
  }
  if (response == "continuous_growth_rate") {
    smap_Phi <- Reduce("%*%", rev(lapply(smap_J_curr, function(A) expm(delta_t * A))))
  }
  smap_singular_phi <- svd(smap_Phi)$d
  smap_max_growth_rate <- log(max(smap_singular_phi)) / tau
  smap_median_growth_rate <- log(sum(diag(smap_Phi %*% t(smap_Phi))) / n_sp) / (2 * tau)
  results_smap_df[t, "max_growth_rate"] <- smap_max_growth_rate
  results_smap_df[t, "median_growth_rate"] <- smap_median_growth_rate
}
# correlation between inferred and true stability metrics
cor(results_smap_df$median_growth_rate, results_analytical_df$median_growth_rate, 
    use = "complete.obs")

# plot time series of true and inferred stability metrics ------------------------------ 
# data frame for plotting
plot_df <- rbind(results_analytical_df, results_smap_df)
plot_df$type <- factor(plot_df$type, levels = c("True", "Inferred"))
# plot median perturbation growth rate through time
fig <- ggplot() +
  geom_hline(yintercept = 0, size = 0.5) +
  geom_line(data = plot_df, aes(x = time, y = median_growth_rate, color = type), 
            size = 0.5) +
  geom_point(data = plot_df, aes(x = time, y = median_growth_rate, fill = type), 
             size = 2.5, shape = 21) +
  scale_color_manual(values = c("gray70", "black")) +
  scale_fill_manual(values = c("gray70", "black")) +
  xlab("Time") +
  ylab("Median perturbation\ngrowth rate") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.6, "cm"))
if (save_plots) {
  ggsave(paste("figs/", func_name, "_time_series_median_pert_growth_rate.pdf", sep = ""), 
         fig, width = 21, height = 6, units = "cm")
}
