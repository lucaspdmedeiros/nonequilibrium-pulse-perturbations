# Function that performs a one-step-ahead forecast of different variables using a linear 
# model, where the coefficients of the model are supplied as inputs to the function

# Arguments:
# state: n-dimensional vector at time t, where n is the number of variables in the system
# coefs: n x n+1 matrix of coefficients, where row i has the coefficients to forecast variable i
# response: type of response variable (abundance, continuous_growth_rate, or discrete_growth_rate)

# Output:
# forecast: n-dimensional vector at time t+1, giving the forecast of each of the n variables

regularized_smap_forecast <- function(state, coefs, response) {
  # perform forecast
  if (response == "abundance") {
    forecast <- as.numeric(coefs %*% c(1, state))
  }
  if (response == "continuous_growth_rate" | response == "discrete_growth_rate") {
    forecast <- state * exp(as.numeric(coefs %*% c(1, state)))
  }
  # return results
  return(forecast)
}
