# Function that performs leave-future-out cross-validation with inferred coefficients
# from a locally weighted linear regularized regression (Regularized S-map) and
# computes the leave-future-out R2 value for each variable

# Arguments:
# training_ts: training set with time in first column and variables in other columns
# test_ts: test set with time in first column and variables in other columns
# theta: state kernel exponential decay strength (theta=0 corresponds to an AR(1) model)
# alpha: proportion of each regularization (alpha=1 is lasso only and alpha=0 is ridge only)
# lambda: strength of regularization penalty
# response: type of response variable (abundance, continuous_growth_rate, or discrete_growth_rate)
# delta_t: time step between observations, if response is in continuous time
# method: which method to use to fit regression (analytical or glmnet)
# scale: whether to scale (TRUE) or not (FALSE) variables to compute weights (variables are always scaled to fit the regression)
# output: whether to output forecasts, R2 for all variables, or mean R2

# Output:
# forecasts: forecasts for each variable over time
# or
# R2: R2 value for each variable
# or
# mean_R2: mean R2 value across all variables

regularized_smap_cv <- function(training_ts, test_ts, theta, alpha, lambda, response, delta_t,
                                method, scale, output) {
  # for each training set, fit regression, forecast future point and compute error of forecast
  forecasts <- matrix(NA, nrow = nrow(test_ts), ncol = ncol(test_ts))
  error_matrix <- matrix(NA, nrow = nrow(test_ts), ncol = ncol(test_ts)-1)
  for (t in 1:nrow(test_ts)) {
    # fit regression (note that we use structure_transform = FALSE for forecasting)
    smap_results <- regularized_smap_fit(ts = training_ts, points = nrow(training_ts), theta = theta, 
                                         alpha = alpha, lambda = lambda, response = response, delta_t = delta_t, 
                                         method = method, scale = scale, structure_transform = FALSE)
    # extract last coefficient matrix
    coefs <- cbind(smap_results[[2]][1, ], smap_results[[1]][[1]])
    # perform forecast
    forecast <- regularized_smap_forecast(state = as.numeric(training_ts[nrow(training_ts), -1]), 
                                          coefs = coefs, response = response) 
    forecasts[t, ] <- c(t, forecast)
    # compute forecast squared error for each variable
    error_matrix[t, ] <- as.numeric((forecast - test_ts[t, -1])^2)
    # increase training set by adding one new point
    training_ts <- rbind(training_ts, test_ts[t, ]) 
  }
  # sum forecast squared error for each variable across time
  error_sum_forecast <- apply(error_matrix, 2, sum, na.rm = TRUE)
  # compute sum of squared deviations from the mean for the test set
  test_ts_mean <- apply(test_ts[ , -1], 2, mean, na.rm = TRUE)
  error_sum_test <- apply(t(t(test_ts[ , -1]) - test_ts_mean)^2, 2, sum, na.rm = TRUE)
  # compute R2 for each variable
  R2 <- 1 - (error_sum_forecast / error_sum_test)
  # compute mean RMSE 
  mean_R2 <- mean(R2, na.rm = TRUE)
  # return results
  if (output == "forecasts") {
    return(forecasts)
  }
  if (output == "R2") {
    return(R2)
  }
  if (output == "mean_R2") {
    return(mean_R2)
  }
}
