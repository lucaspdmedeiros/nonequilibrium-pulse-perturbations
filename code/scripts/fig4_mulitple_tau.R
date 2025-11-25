# Similar to the code for Figure 4 but for multiple values of tau to see how
# expected and maximum growth rate change with tau

# cleaning environment and loading functions ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/rosenzweig_macarthur_2sp_ode.R")
source("code/functions/hastings_powell_3sp_ode.R")

# settings ------------------------------
# models to use
func_names <- c("rosenzweig_macarthur_2sp_forced_fixed_point",
                "rosenzweig_macarthur_2sp_limit_cycle",
                "hastings_powell_3sp_forced_cycle",
                "hastings_powell_3sp_chaos")
# data frame to store results
full_results_df <- data.frame()
# list of tau values to use
tau_list <- seq(0.1, 1, by = 0.1)
# to reproduce results 
set.seed(42)
# whether to run analysis or just load saved results
run_analysis <- TRUE

# perform analysis for each model ------------------------------
if (run_analysis) {
  for (i in 1:length(func_names)) {
    # data frame to store results
    results_df <- data.frame()
    # current model
    func_name <- func_names[i]
    print(func_name)
    # load model settings
    source("code/scripts/settings.R")
    # species labels
    sp_names <- paste("x", 1:n_sp, sep = "")
    # number of points to use for analysis
    sample_size <- 50 * 2^n_sp
    
    # generate non-perturbed abundance trajectory ------------------------------
    # generate time series
    ts <- as.data.frame(ode(y = state, times = times, func = func, parms = parms, method = "ode45"))
    n_points_keep <- floor(nrow(ts) * ((n_recurrences - 1) / n_recurrences))
    # subset of points to store results
    ts_sub <- head(ts, n_points_keep)[floor(seq(1, n_points_keep, length = sample_size)), ]
    
    # compute Jacobian matrix along trajectory ------------------------------
    if ((func_name == "rosenzweig_macarthur_2sp_forced_fixed_point") | 
        (func_name == "rosenzweig_macarthur_2sp_limit_cycle")) {
      J <- list()
      x1 <- ts$x1
      x2 <- ts$x2
      for (j in 1:nrow(ts)) {
        J[[j]] <- matrix(c(r - ((2 * r * x1[j]) / k[j]) - ((a * x2[j] * (b + x1[j]) - a * x1[j] * x2[j]) / (b + x1[j])^2), 
                           - (a * x1[j]) / (b + x1[j]),
                           (e * a * x2[j] * (b + x1[j]) - e * a * x1[j] * x2[j]) / (b + x1[j])^2, 
                           ((e * a * x1[j]) / (b + x1[j])) - d), nrow = 2, byrow = TRUE)
      }
    }
    if ((func_name == "hastings_powell_3sp_forced_cycle") |
        (func_name == "hastings_powell_3sp_chaos")) {
      J <- list()
      x1 <- ts$x1
      x2 <- ts$x2
      x3 <- ts$x3
      for (j in 1:nrow(ts)) {
        J[[j]] <- matrix(c((1 - (x1[j] / k)) - x1[j] * (1 / k) - (a1 * x2[j] / (x1[j] + b1) - (a1 * x1[j] * x2[j]) / (x1[j] + b1)^2),
                           -(a1 * x1[j] / (x1[j] + b1)),
                           0,
                           a1 * e1 * x2[j] / (x1[j] + b1) - (a1 * e1 * x1[j] * x2[j]) / (x1[j] + b1)^2,
                           a1 * e1 * x1[j] / (x1[j] + b1) - d1 - (a2[j] * x3[j] / (x2[j] + b2) - (a2[j] * x2[j] * x3[j]) / (x2[j] + b2)^2),
                           -(a2[j] * x2[j] / (x2[j] + b2)),
                           0,
                           a2[j] * e2 * x3[j] / (x2[j] + b2) - (a2[j] * e2 * x2[j] * x3[j]) / (x2[j] + b2)^2,
                           a2[j] * e2 * x2[j] / (x2[j] + b2) - d2), nrow = 3, byrow = TRUE)
      }
    }
    names(J) <- sprintf("%.5f", ts$time)
    # loop over values of tau
    for (j in 1:length(tau_list)) {
      print(j)
      # percentage of recurrence time to evolve perturbations
      perc_rec_time <- tau_list[j]
      tau <- time_step * n_points * perc_rec_time
      
      # compute perturbation growth rate along trajectory ------------------------------
      mean_ana_results_df <- data.frame()
      Sigma <- diag(rep(1, n_sp))
      min_growth_rate <- c() 
      max_growth_rate <- c() 
      expect_avg_growth_rate <- c()
      for (k in 1:nrow(ts_sub)) {
        curr_time <- ts_sub$time[k]
        curr_state <- as.numeric(ts_sub[k, sp_names])
        names(curr_state) <- sp_names
        J_curr <- J[sprintf("%.5f", seq(curr_time, (curr_time + tau), by = time_step))]
        Phi <- Reduce("%*%", rev(lapply(J_curr, function(A) expm(time_step * A))))
        min_singular <- min(Re(eigen(Phi %*% t(Phi))$values))
        max_singular <- max(Re(eigen(Phi %*% t(Phi))$values))
        min_growth_rate[k] <- log(min_singular) / (2 * tau)
        max_growth_rate[k] <- log(max_singular) / (2 * tau)
        expect_avg_growth_rate[k] <- (log(sum(diag(Phi %*% Sigma %*% t(Phi)))) - log(sum(diag(Sigma)))) / (2 * tau)
      }
      # save current results
      plot_df <- ts_sub
      plot_df$min_growth_rate <- min_growth_rate
      plot_df$max_growth_rate <- max_growth_rate
      plot_df$expect_avg_growth_rate <- expect_avg_growth_rate
      plot_df$tau <- perc_rec_time
      results_df <- rbind(results_df, plot_df)
    }
    if (n_sp == 2) {
      results_df$x3 <- NA
    }
    results_df$model <- func_name
    full_results_df <- rbind(full_results_df, results_df)
  }
  # save results
  write.csv(x = full_results_df, file = "tables/perturbation_growth_rate_multiple_tau_all_models.csv", row.names = FALSE)
} else {
  # load results from previous analysis
  full_results_df <- read.csv(file = "tables/perturbation_growth_rate_multiple_tau_all_models.csv")
}

# plots of how expected and maximum growth rate change with tau  ------------------------------
# compute maximum perturbation growth rate under for each combination of model and tau
summ_results_df <- ddply(full_results_df, c("tau", "model"), summarise,
                         max_max_growth_rate = max(max_growth_rate),
                         max_expect_avg_growth_rate = max(expect_avg_growth_rate))
# rescale maximum expected growth rate by the maximum within each model
summ_summ_results_df <- ddply(summ_results_df, c("model"), summarise,
                              max_expect_avg_growth_rate = max_expect_avg_growth_rate / max(max_expect_avg_growth_rate))
summ_summ_results_df$tau <- tau_list
summ_summ_results_df$model <- factor(summ_summ_results_df$model, 
                                     levels = c("rosenzweig_macarthur_2sp_forced_fixed_point",
                                                "rosenzweig_macarthur_2sp_limit_cycle",
                                                "hastings_powell_3sp_forced_cycle",
                                                "hastings_powell_3sp_chaos"))
# plot for expected growth rate
fig <- ggplot() +
  geom_point(data = summ_summ_results_df, aes(x = tau, y = max_expect_avg_growth_rate, color = model), 
             size = 4) +
  geom_line(data = summ_summ_results_df, aes(x = tau, y = max_expect_avg_growth_rate, color = model), 
             size = 0.7) +
  xlab(latex2exp::TeX(r"(Fraction of recurrence time)")) +
  ylab(latex2exp::TeX(r"(Relative maximum $\bar{r}_{\tau}$ along trajectory)")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.8, "cm"))
# rescale maximum maximum growth rate by the maximum within each model
summ_summ_results_df <- ddply(summ_results_df, c("model"), summarise,
                              max_max_growth_rate = max_max_growth_rate / max(max_max_growth_rate))
summ_summ_results_df$tau <- tau_list
summ_summ_results_df$model <- factor(summ_summ_results_df$model, 
                                     levels = c("rosenzweig_macarthur_2sp_forced_fixed_point",
                                                "rosenzweig_macarthur_2sp_limit_cycle",
                                                "hastings_powell_3sp_forced_cycle",
                                                "hastings_powell_3sp_chaos"))
# plot for maximum growth rate
fig <- ggplot() +
  geom_point(data = summ_summ_results_df, aes(x = tau, y = max_max_growth_rate, color = model), 
             size = 4) +
  geom_line(data = summ_summ_results_df, aes(x = tau, y = max_max_growth_rate, color = model), 
            size = 0.7) +
  xlab(latex2exp::TeX(r"(Fraction of recurrence time)")) +
  ylab(latex2exp::TeX(r"(Relative maximum $\max (r_{\tau})$ along trajectory)")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.8, "cm"))
