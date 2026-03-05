# t-Rec: Reconciliation via Conditioning with uncertain covariance via Multivariate Student-t

Reconciles base forecasts in a hierarchy by conditioning on the
hierarchical constraints, specified by the aggregation matrix A. The
base forecasts are assumed to be jointly Gaussian, conditionally on the
covariance matrix, a Bayesian approach is adopted using an
Inverse-Wishart prior, leading to a Multivariate Student-t distribution
for the base forecasts. The reconciliation is in closed-form, yielding a
multivariate Student-t reconciled distribution.

## Usage

``` r
reconc_t(
  A,
  base_fc_mean,
  y_train = NULL,
  residuals = NULL,
  ...,
  return_upper = FALSE,
  return_parameters = FALSE
)
```

## Arguments

- A:

  Matrix (n_upp x n_bott) defining the hierarchy (u = Ab).

- base_fc_mean:

  Vector of base forecasts (length n = n_upp + n_bott).

- y_train:

  mts (or matrix) of historical training data (T x n) used for setting
  prior parameters.

- residuals:

  Matrix (T x n) of base forecast residuals.

- ...:

  Additional arguments for advanced usage: see details.

- return_upper:

  Logical; if TRUE, also returns parameters for the upper level
  reconciled distribution.

- return_parameters:

  Logical; if TRUE, returns internal parameters like C and posterior nu.

## Value

A list containing:

- `bottom_rec_mean`: Reconciled bottom-level mean forecasts.

- `bottom_rec_scale_matrix`: Reconciled bottom-level scale matrix.

- `bottom_rec_df`: Reconciled degrees of freedom.

If `return_upper` is TRUE, also returns:

- `upper_rec_mean`: Reconciled upper-level mean forecasts.

- `upper_rec_scale_matrix`: Reconciled upper-level scale matrix.

- `upper_rec_df`: Reconciled upper-level degrees of freedom.

If `return_parameters` is TRUE, also returns:

- `prior_nu`: Prior degrees of freedom.

- `posterior_nu`: Posterior degrees of freedom.

- `posterior_Psi`: Posterior scale matrix.

- `C`: Scaling factor for the scale matrix.

## Details

**Standard Usage and Parameter Estimation:** The standard workflow for
this function is to provide the in-sample `residuals` and the historical
training data `y_train`.

- **Prior Scale (\\\Psi_0\\):** Set as the covariance of the residuals
  of naive (or seasonal naive, a criterion is used to choose between
  the 2) forecasts computed on `y_train`.

- **Prior Degrees of Freedom (\\\nu_0\\):** Estimated via Bayesian
  Leave-One-Out Cross-Validation (LOOCV) to maximize out-of-sample
  performance.

**Advanced Options:** Users can bypass the automated estimation by:

1.  Directly passing the `prior` parameters as a list with entries 'nu'
    and 'Psi'. This skips the LOOCV step for \\\nu_0\\ and the
    covariance estimation from `y_train`. It requires `residuals` to
    compute the posterior.

2.  Directly passing the `posterior` parameters as a list with entries
    'nu' and 'Psi'. This skips all internal estimation and updating
    logic.

Moreover, users can specify:

- `freq`: positive integer, used as frequency of data for the seasonal
  naive forecast in the specification of \\\Psi_0\\. By default, if
  `y_train` is a multivariate time series, the frequency of the data is
  used; otherwise, it is set to 1 (no seasonality).

- `criterion`: either 'RSS' (default) or 'seas-test', specifying which
  criterior is used to choose between the naive and seasonal naive
  forecasts for the specification of \\\Psi_0\\. 'RSS' computes the
  residual sum of squares for both methods and chooses the one with
  lower RSS, while 'seas-test' uses a statistical test for seasonality
  (currently implemented using the number of seasonal differences
  suggested by the `forecast` package, which must be installed).

**The Reconciled Bottom Distribution:** The reconciliation yields a
distribution: \$\$\tilde{\mathbf{b}} \sim t(\hat{\mathbf{b}}\_{tilde},
\tilde{\Sigma}\_B, \tilde{\nu})\$\$ where the reconciled mean is:
\$\$\hat{\mathbf{b}}\_{tilde} = \hat{\mathbf{b}} + ((\Psi'\_{UB})^\top -
\Psi'\_B A^\top) Q^{-1} (A\hat{\mathbf{b}} - \hat{\mathbf{u}})\$\$ and
the scale matrix is: \$\$\tilde{\Sigma}\_B = C \[\Psi'\_B -
((\Psi'\_{UB})^\top - \Psi'\_B A^\top) Q^{-1} ((\Psi'\_{UB})^\top -
\Psi'\_B A^\top)^\top\]\$\$ with scalar \$\$C = \frac{1 +
(A\hat{\mathbf{b}} - \hat{\mathbf{u}})^\top Q^{-1} (A\hat{\mathbf{b}} -
\hat{\mathbf{u}})}{\tilde{\nu}}.\$\$

## References

Carrara, C., Corani, G., Azzimonti, D., & Zambon, L. (2025). Modeling
the uncertainty on the covariance matrix for probabilistic forecast
reconciliation. arXiv preprint arXiv:2506.19554.
<https://arxiv.org/abs/2506.19554>

## See also

[`reconc_gaussian()`](https://idsia.github.io/bayesRecon/reference/reconc_gaussian.md)

## Examples

``` r
# \donttest{
library(bayesRecon)

if (requireNamespace("forecast", quietly = TRUE)) {

  set.seed(1234)
  n_obs <- 100

  # Simulate 2 bottom series from AR(1) processes
  y1 <- arima.sim(model = list(ar = 0.8), n = n_obs)
  y2 <- arima.sim(model = list(ar = 0.5), n = n_obs)

  y_upper <- y1 + y2   # upper series is the sum of the two bottoms
  A <- matrix(c(1, 1), nrow = 1)  # Aggregation matrix

  # Fit additive ETS models
  fit1 <- forecast::ets(y1, additive.only = TRUE)
  fit2 <- forecast::ets(y2, additive.only = TRUE)
  fit_upper <- forecast::ets(y_upper, additive.only = TRUE)

  # Point forecasts (h = 1)
  fc_upper <- as.numeric(forecast::forecast(fit_upper, h = 1)$mean)
  fc1 <- as.numeric(forecast::forecast(fit1, h = 1)$mean)
  fc2 <- as.numeric(forecast::forecast(fit2, h = 1)$mean)
  base_fc_mean <- c(fc_upper, fc1, fc2)

  # Residuals and training data (n_obs x n matrices, columns in same order as base_fc_mean)
  res <- cbind(residuals(fit_upper), residuals(fit1), residuals(fit2))
  y_train <- cbind(y_upper, y1, y2)

  # --- 1) Generate joint reconciled samples ---
  result <- reconc_t(A, base_fc_mean, y_train = y_train, residuals = res)

  # Sample from the reconciled bottom-level Student-t distribution
  n_samples <- 2000
  L_chol <- t(chol(result$bottom_rec_scale_matrix))
  z <- matrix(rt(ncol(A) * n_samples, df = result$bottom_rec_df), nrow = ncol(A))
  bottom_samples <- result$bottom_rec_mean + L_chol %*% z  # 2 x n_samples

  # Aggregate bottom samples to get upper samples
  upper_samples <- A %*% bottom_samples              
  joint_samples <- rbind(upper_samples, bottom_samples)
  rownames(joint_samples) <- c("upper", "bottom_1", "bottom_2")

  cat("Reconciled means (from samples):\n")
  print(round(rowMeans(joint_samples), 3))

  cat("Reconciled standard deviations (from samples):\n")
  print(round(apply(joint_samples, 1, sd), 3))

  # --- 2) 95% prediction intervals via t-distribution quantiles ---
  result2 <- reconc_t(A, base_fc_mean, y_train = y_train,
                      residuals = res, return_upper = TRUE)

  alpha <- 0.05
  # Bottom series intervals
  for (i in seq_len(ncol(A))) {
    s_i <- sqrt(result2$bottom_rec_scale_matrix[i, i])
    lo  <- result2$bottom_rec_mean[i] + s_i * qt(alpha / 2,     df = result2$bottom_rec_df)
    hi  <- result2$bottom_rec_mean[i] + s_i * qt(1 - alpha / 2, df = result2$bottom_rec_df)
    cat(sprintf("Bottom %d: 95%% PI = [%.3f, %.3f]\n", i, lo, hi))
  }
  # Upper series interval
  s_u <- sqrt(result2$upper_rec_scale_matrix[1, 1])
  lo  <- result2$upper_rec_mean[1] + s_u * qt(alpha / 2,     df = result2$upper_rec_df)
  hi  <- result2$upper_rec_mean[1] + s_u * qt(1 - alpha / 2, df = result2$upper_rec_df)
  cat(sprintf("Upper:    95%% PI = [%.3f, %.3f]\n", lo, hi))
}
#> Reconciled means (from samples):
#>    upper bottom_1 bottom_2 
#>   -2.370   -0.056   -2.314 
#> Reconciled standard deviations (from samples):
#>    upper bottom_1 bottom_2 
#>    1.574    0.947    1.212 
#> Bottom 1: 95% PI = [-1.917, 1.777]
#> Bottom 2: 95% PI = [-4.640, 0.062]
#> Upper:    95% PI = [-5.377, 0.659]
# }
```
