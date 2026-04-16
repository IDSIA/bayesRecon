# Analytical reconciliation of Gaussian base forecasts

Closed form computation of the reconciled forecasts in case of Gaussian
base forecasts.

## Usage

``` r
reconc_gaussian(
  A,
  base_fc_mean,
  base_fc_cov = NULL,
  residuals = NULL,
  return_upper = FALSE
)
```

## Arguments

- A:

  aggregation matrix (n_upper x n_bottom).

- base_fc_mean:

  a vector containing the means of the base forecasts.

- base_fc_cov:

  a matrix containing the covariance matrix of the base forecasts.

- residuals:

  a matrix with the residuals of the base forecasts, with n_upper +
  n_bottom columns. The covariance matrix of the base forecasts is
  computed from the residuals using the Schäfer Strimmer shrinkage
  estimator. If base_fc_cov is provided, residuals are ignored.

- return_upper:

  logical, whether to return the reconciled parameters for the upper
  variables (default is FALSE).

## Value

A list containing the reconciled forecasts. The list has the following
named elements:

- `bottom_rec_mean`: reconciled mean for the bottom forecasts;

- `bottom_rec_cov`: reconciled covariance for the bottom forecasts;

- `upper_rec_mean`: (only if `return_upper = TRUE`) reconciled mean for
  the upper forecasts;

- `upper_rec_cov`: (only if `return_upper = TRUE`) reconciled covariance
  for the upper forecasts.

## Details

In the vector of the means of the base forecasts the order must be:
first the upper, then the bottom; the order within the uppers is given
by the rows of A, the order within the bottoms by the columns of A. The
order of the rows of the covariance matrix of the base forecasts is the
same.

Unless `return_upper = TRUE`, the function returns only the reconciled
parameters of the bottom variables.

## References

Corani, G., Azzimonti, D., Augusto, J.P.S.C., Zaffalon, M. (2021).
*Probabilistic Reconciliation of Hierarchical Forecast via Bayes' Rule*.
ECML PKDD 2020. Lecture Notes in Computer Science, vol 12459.
[doi:10.1007/978-3-030-67664-3_13](https://doi.org/10.1007/978-3-030-67664-3_13)
.

Zambon, L., Agosto, A., Giudici, P., Corani, G. (2024). *Properties of
the reconciled distributions for Gaussian and count forecasts*.
International Journal of Forecasting (in press).
[doi:10.1016/j.ijforecast.2023.12.004](https://doi.org/10.1016/j.ijforecast.2023.12.004)
.

## See also

[`reconc_t()`](https://idsia.github.io/bayesRecon/reference/reconc_t.md),
[`reconc_BUIS()`](https://idsia.github.io/bayesRecon/reference/reconc_BUIS.md)

## Examples

``` r
library(bayesRecon)

#' # ---- Example 1: base forecasts are given ----

# Create a minimal hierarchy with 2 bottom and 1 upper variable
A <- get_reconc_matrices(agg_levels = c(1, 2), h = 2)$A

# Set the parameters of the Gaussian base forecast distributions
mu1 <- 2
mu2 <- 4
muY <- 9
base_fc_mean <- c(muY, mu1, mu2)  # vector of means

sigma1 <- 2
sigma2 <- 2
sigmaY <- 3
sigmas <- c(sigmaY, sigma1, sigma2)
base_fc_cov <- diag(sigmas^2)  # covariance matrix

analytic_rec <- reconc_gaussian(A,
  base_fc_mean = base_fc_mean,
  base_fc_cov = base_fc_cov
)

bottom_mean_rec <- analytic_rec$bottom_rec_mean
bottom_cov_rec <- analytic_rec$bottom_rec_cov

# To obtain reconciled samples for the entire hierarchy, sample from the reconciled 
# bottom distribution and then aggregate using A. 

# Sample from the reconciled bottom-level Gaussian distribution
# First, compute the Cholesky decomposition of the reconciled covariance matrix:
chol_decomp <- chol(bottom_cov_rec)
# Then, sample from the standard normal distribution and apply the transformation:
Z <- matrix(stats::rnorm(n = 2000), nrow = 2) 
B <- t(chol_decomp) %*% Z + matrix(rep(bottom_mean_rec, 1000), nrow = 2) 

# Aggregate bottom samples to get upper samples, then stack
U <- A %*% B
Y_reconc <- rbind(U, B)

cat("Dimensions of reconciled samples (upper + bottom):", dim(Y_reconc), "\n")
#> Dimensions of reconciled samples (upper + bottom): 3 1000 


# ---- Example 2: using residuals from fitted ETS models ----
# \donttest{
if (requireNamespace("forecast", quietly = TRUE)) {

  # Simulate 2 bottom series from AR(1) processes
  set.seed(1234)
  n_obs <- 200
  y1 <- arima.sim(model = list(ar = 0.8), n = n_obs)
  y2 <- arima.sim(model = list(ar = 0.5), n = n_obs)

  # Upper series is the sum of the two bottom series
  y_upper <- y1 + y2

  # Aggregation matrix A:
  A <- matrix(c(1, 1), nrow = 1)

  # Fit additive ETS models
  fit1 <- forecast::ets(y1, additive.only = TRUE)
  fit2 <- forecast::ets(y2, additive.only = TRUE)
  fit_upper <- forecast::ets(y_upper, additive.only = TRUE)

  # Point forecasts (h = 1):
  fc1 <- forecast::forecast(fit1, h = 1)$mean
  fc2 <- forecast::forecast(fit2, h = 1)$mean
  fc_upper <- forecast::forecast(fit_upper, h = 1)$mean
  base_fc_mean <- c(fc_upper, fc1, fc2)

  # Residuals matrix (T x n, columns in same order as base_fc_mean)
  res <- cbind(residuals(fit_upper),
               residuals(fit1),
               residuals(fit2))

  # Reconcile (covariance estimated internally via Schafer-Strimmer)
  result <- reconc_gaussian(A, base_fc_mean = base_fc_mean, residuals = res, return_upper = TRUE)

  bottom_mean <- result$bottom_rec_mean
  bottom_cov <- result$bottom_rec_cov
  upper_mean <- result$upper_rec_mean
  upper_cov <- result$upper_rec_cov

  # Print reconciled means
  cat("Reconciled bottom means:", round(bottom_mean, 3), "\n")
  cat("Reconciled upper mean:", round(upper_mean, 3), "\n")

  # Print 95% predictions intervals
  cat("Reconciled bottom 95% prediction intervals:\n")
  for (i in 1:length(bottom_mean)) {
    lower <- bottom_mean[i] - 1.96 * sqrt(bottom_cov[i, i])
    upper <- bottom_mean[i] + 1.96 * sqrt(bottom_cov[i, i])
    cat(paste0("Bottom ", i, ": [", round(lower, 3), ", ", round(upper, 3), "]\n"))
  }
  cat("Reconciled upper 95% prediction interval:\n")
  lower <- upper_mean - 1.96 * sqrt(upper_cov[1, 1])
  upper <- upper_mean + 1.96 * sqrt(upper_cov[1, 1])
  cat(paste0("Upper: [", round(lower, 3), ", ", round(upper, 3), "]\n"))

}
#> Reconciled bottom means: 4.354 1.632 
#> Reconciled upper mean: 5.986 
#> Reconciled bottom 95% prediction intervals:
#> Bottom 1: [2.288, 6.419]
#> Bottom 2: [-0.477, 3.742]
#> Reconciled upper 95% prediction interval:
#> Upper: [3, 8.972]
# }
```
