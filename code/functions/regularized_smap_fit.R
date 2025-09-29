# Function that fits a locally weighted linear regularized regression (Regularized S-map) 
# to infer the Jacobian matrix of a dynamical system at different points in time using 
# a multivariate time series

# Arguments:
# ts: time series with time in first column and variables in other columns
# points: which time series points to use as target points
# theta: state kernel exponential decay strength (theta=0 corresponds to an AR(1) model)
# alpha: proportion of each regularization (alpha=1 is lasso only and alpha=0 is ridge only)
# lambda: strength of regularization penalty
# response: type of response variable (abundance, continuous_growth_rate, or discrete_growth_rate)
# delta_t: time step between observations, if response is in continuous time
# method: which method to use to fit regression (analytical or glmnet)
# scale: whether to scale (TRUE) or not (FALSE) variables to compute weights (variables are always scaled to fit the regression)
# structure_transform: how to transform coefficients according to structure of dynamics (none, growth_rate, or model_prediction)

# Output list elements:
# J: list with inferred Jacobian matrices
# intercepts: matrix with inferred intercept values

regularized_smap_fit <- function(ts, points, theta, alpha, lambda, response, delta_t, 
                                 method, scale, structure_transform) {
  # number of points
  n <- length(points)
  # number of variables
  n_var <- ncol(ts) - 1
  # to save Jacobian matrices and intercepts
  J <- lapply(rep(NA, n), matrix, nrow = n_var, ncol = n_var)
  intercepts <- matrix(NA, nrow = n, ncol = n_var)
  # time series at time t 
  ts_past <- as.matrix(ts[ , -1])
  # time series at time t+1
  ts_future <- as.matrix(ts[-1, -1])
  ts_future <- rbind(ts_future, rep(NA, n_var))
  # matrices with predictor and response variables
  if (response == "abundance") {
    Y <- ts_future
    X <- cbind(1, ts_past)
  }
  if (response == "continuous_growth_rate") {
    ts_growth_rate <- log(ts_future / ts_past) / delta_t
    Y <- ts_growth_rate
    X <- cbind(1, ts_past)
  }
  if (response == "discrete_growth_rate") {
    ts_growth_rate <- log(ts_future / ts_past)
    Y <- ts_growth_rate
    X <- cbind(1, ts_past)
  }
  # scale variables to mean 0 and standard deviation 1 to compute weights
  Y_mean <- apply(Y, 2, mean, na.rm = TRUE)
  Y_sd <- apply(Y, 2, sd, na.rm = TRUE)
  X_mean <- apply(X[ , -1], 2, mean, na.rm = TRUE)
  X_sd <- apply(X[ , -1], 2, sd, na.rm = TRUE)
  if (scale) {
    Y_std <- t((t(Y) - Y_mean) / Y_sd)
    X_std <- X
    X_std[ , -1] <- t((t(X[ , -1]) - X_mean) / X_sd)
  } else {
    Y_std <- Y
    X_std <- X
  }
  # add zero to rows with NA for glmnet to work (weights for these rows will also be zero)
  which_na_X <- apply(X, 2, function(x) which(!is.finite(x)))
  which_na_Y <- apply(Y, 2, function(x) which(!is.finite(x)))
  which_na <- sort(unique(c(unlist(which_na_X), unlist(which_na_Y))))
  Y[which_na, ] <- 0
  X[which_na, ] <- 0
  # fit regression for each point in the time series
  for (t in 1:n) {
    # current state (note that X_std still has NAs)
    curr_state <- X_std[points[t], -1]
    # proceed with regression if none of variables is NA
    if (all(!is.na(curr_state))) {
      # weights for each point based on distance in state space
      state_dist <- as.numeric(apply(X_std[ , -1], 1, function(x) sqrt(sum((curr_state - x)^2))))
      avg_state_dist <- mean(state_dist, na.rm = TRUE)
      weights <- exp(-theta * (state_dist / avg_state_dist))
      # set weights for rows with NA to zero
      weights[which_na] <- 0
      # analytical ridge regression solution
      if ((method == "analytical") & (alpha == 0)) {
        Y_std <- t((t(Y) - Y_mean) / Y_sd)
        X_std[ , -1] <- t((t(X[ , -1]) - X_mean) / X_sd)
        Y_std <- Y_std * weights
        X_std <- X_std * weights
        # fit regression
        regression_coefs <- t(solve(t(X_std) %*% X_std + lambda*diag(1, n_var+1)) %*% t(X_std) %*% Y_std)
        # transform intercepts and regression coefficients back to original scale
        for (i in 1:n_var) {
          intercepts[t, i] <- as.numeric(Y_sd[i] * regression_coefs[i, 1] + Y_mean[i] - 
                                           sum((Y_sd[i] / X_sd) * regression_coefs[i, -1] * X_mean))
          J[[t]][i, ] <- (Y_sd[i] / X_sd) * regression_coefs[i, -1]
        }
      }
      # elastic net solution
      if (method == "glmnet") {
        # fit regression
        elastic_net_fit <- glmnet(x = X[ , -1], y = Y, family = "mgaussian", weights = weights, 
                                  alpha = alpha, lambda = lambda, standardize = TRUE, intercept = TRUE)
        # extract coefficients
        regression_coefs <- matrix(NA, nrow = n_var, ncol = n_var + 1)
        regression_coefs[ , 1] <- elastic_net_fit$a0
        for (i in 1:length(elastic_net_fit$beta)) {
          regression_coefs[i, -1] <- as.numeric(unlist(elastic_net_fit$beta[[i]]))
        }
        intercepts[t, ] <- as.numeric(regression_coefs[ , 1])
        J[[t]] <- regression_coefs[ , -1]
      }
      # transform regression coefficients according to structure of dynamics
      if (structure_transform == "growth_rate") {
        if (response == "continuous_growth_rate") {
          J[[t]] <- J[[t]] * as.numeric(ts[points[t], -1])
          diag(J[[t]]) <- as.numeric(diag(J[[t]]) + log(ts[points[t]+1, -1] / ts[points[t], -1]) / delta_t)
        }
        if (response == "discrete_growth_rate") {
          J[[t]] <- J[[t]] * as.numeric(ts[points[t]+1, -1])
          diag(J[[t]]) <- as.numeric(diag(J[[t]]) + (ts[points[t]+1, -1] / ts[points[t], -1]))
        }
      }
      if (structure_transform == "model_prediction") {
        if (response == "continuous_growth_rate") {
          J[[t]] <- J[[t]] * as.numeric(ts[points[t], -1])
          diag(J[[t]]) <- as.numeric(diag(J[[t]]) + regularized_smap_forecast(state = as.numeric(ts[points[t], -1]), 
                                                                              coefs = regression_coefs, response = "abundance")) 
        }
        if (response == "discrete_growth_rate") {
          J[[t]] <- J[[t]] * as.numeric(ts[points[t]+1, -1])
          diag(J[[t]]) <- as.numeric(diag(J[[t]]) + exp(regularized_smap_forecast(state = as.numeric(ts[points[t], -1]), 
                                                                                  coefs = regression_coefs, response = "abundance")))
        }
      }
    }
  }
  # return results
  return(list(J, intercepts))
}
