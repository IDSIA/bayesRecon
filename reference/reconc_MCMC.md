# MCMC for Probabilistic Reconciliation of forecasts via conditioning

Uses Markov Chain Monte Carlo algorithm to draw samples from the
reconciled forecast distribution, which is obtained via conditioning.

This is a bare-bones implementation of the Metropolis-Hastings
algorithm, we suggest the usage of tools to check the convergence. The
function only works with Poisson or Negative Binomial base forecasts.

The function
[`reconc_BUIS()`](https://idsia.github.io/bayesRecon/reference/reconc_BUIS.md)
is generally faster on most hierarchies.

## Usage

``` r
reconc_MCMC(
  A,
  base_fc,
  distr,
  num_samples = 10000,
  tuning_int = 100,
  init_scale = 1,
  burn_in = 1000,
  return_upper = TRUE,
  seed = NULL
)
```

## Arguments

- A:

  aggregation matrix (n_upper x n_bottom).

- base_fc:

  list of the parameters of the base forecast distributions, see
  details.

- distr:

  a string describing the type of predictive distribution.

- num_samples:

  number of samples to draw using MCMC.

- tuning_int:

  number of iterations between scale updates of the proposal.

- init_scale:

  initial scale of the proposal.

- burn_in:

  number of initial samples to be discarded.

- return_upper:

  Logical. If `TRUE` (default), also returns the reconciled samples for
  the upper time series. Default is `TRUE`.

- seed:

  seed for reproducibility.

## Value

A list containing the reconciled forecasts. The list has the following
named elements:

- `bottom_rec_samples`: a matrix (n_bottom x `num_samples`) containing
  reconciled samples for the bottom time series;

- `upper_rec_samples`: a matrix (n_upper x `num_samples`) containing
  reconciled samples for the upper time series (only if
  `return_upper = TRUE`).

## Details

The parameter `base_fc` is a list containing n = n_upper + n_bottom
elements. Each element is a list containing the estimated:

- mean and sd for the Gaussian base forecast, see
  [Normal](https://rdrr.io/r/stats/Normal.html), if `distr`='gaussian';

- lambda for the Poisson base forecast, see
  [Poisson](https://rdrr.io/r/stats/Poisson.html), if `distr`='poisson';

- size and prob (or mu) for the negative binomial base forecast, see
  [NegBinomial](https://rdrr.io/r/stats/NegBinomial.html), if
  `distr`='nbinom'.

The first n_upper elements of the list are the upper base forecasts, in
the order given by the rows of A. The elements from n_upper+1 until the
end of the list are the bottom base forecasts, in the order given by the
columns of A.

## References

Corani, G., Azzimonti, D., Rubattu, N. (2024). *Probabilistic
reconciliation of count time series*. International Journal of
Forecasting 40 (2), 457-469.
[doi:10.1016/j.ijforecast.2023.04.003](https://doi.org/10.1016/j.ijforecast.2023.04.003)
.

## See also

[`reconc_BUIS()`](https://idsia.github.io/bayesRecon/reference/reconc_BUIS.md)

## Examples

``` r
library(bayesRecon)

# Create a minimal hierarchy with 2 bottom and 1 upper variable
rec_mat <- get_reconc_matrices(agg_levels = c(1, 2), h = 2)
A <- rec_mat$A

# Set the parameters of the Poisson base forecast distributions
lambda1 <- 2
lambda2 <- 4
lambdaY <- 9
lambdas <- c(lambdaY, lambda1, lambda2)

base_fc <- list()
for (i in 1:length(lambdas)) {
  base_fc[[i]] <- list(lambda = lambdas[i])
}

# Sample from the reconciled forecast distribution using MCMC
mcmc <- reconc_MCMC(A, base_fc,
  distr = "poisson",
  num_samples = 30000, seed = 42
)
samples_mcmc <- rbind(mcmc$upper_rec_samples, mcmc$bottom_rec_samples)

# Compare the reconciled means with those obtained via BUIS
buis <- reconc_BUIS(A, base_fc,
  in_type = "params",
  distr = "poisson", num_samples = 100000, seed = 42
)
samples_buis <- rbind(buis$upper_rec_samples, buis$bottom_rec_samples)

print(rowMeans(samples_mcmc))
#> [1] 7.091733 2.363633 4.728100
print(rowMeans(samples_buis))
#> [1] 7.09171 2.36542 4.72629
```
