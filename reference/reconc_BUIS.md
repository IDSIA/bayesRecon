# BUIS for Probabilistic Reconciliation of forecasts via conditioning

Uses the Bottom-Up Importance Sampling algorithm to draw samples from
the reconciled forecast distribution, obtained via conditioning.

## Usage

``` r
reconc_BUIS(
  A,
  base_fc,
  in_type,
  distr,
  num_samples = 20000,
  suppress_warnings = FALSE,
  return_upper = TRUE,
  seed = NULL
)
```

## Arguments

- A:

  aggregation matrix (n_upper x n_bottom).

- base_fc:

  A list containing the base_forecasts, see details.

- in_type:

  A string or a list of length n_upper + n_bottom. If it is a list the
  i-th element is a string with two possible values:

  - 'samples' if the i-th base forecasts are in the form of samples;

  - 'params' if the i-th base forecasts are in the form of estimated
    parameters.

  If it `in_type` is a string it is assumed that all base forecasts are
  of the same type.

- distr:

  A string or a list of length n_upper + n_bottom describing the type of
  base forecasts. If it is a list the i-th element is a string with two
  possible values:

  - 'continuous' or 'discrete' if `in_type[[i]]`='samples';

  - 'gaussian', 'poisson' or 'nbinom' if `in_type[[i]]`='params'.

  If `distr` is a string it is assumed that all distributions are of the
  same type.

- num_samples:

  Number of samples drawn from the reconciled distribution. This is
  ignored if `bottom_in_type='samples'`; in this case, the number of
  reconciled samples is equal to the number of samples of the base
  forecasts.

- suppress_warnings:

  Logical. If `TRUE`, no warnings about effective sample size are
  triggered. If `FALSE`, warnings are generated. Default is `FALSE`. See
  Details.

- return_upper:

  Logical, whether to return the reconciled parameters for the upper
  variables (default is TRUE).

- seed:

  Seed for reproducibility.

## Value

A list containing the reconciled forecasts. The list has the following
named elements:

- `bottom_rec_samples`: a matrix (n_bottom x `num_samples`) containing
  the reconciled samples for the bottom time series;

- `upper_rec_samples`: (only if `return_upper = TRUE`) a matrix (n_upper
  x `num_samples`) containing the reconciled samples for the upper time
  series.

## Details

The parameter `base_fc` is a list containing n = n_upper + n_bottom
elements. The first n_upper elements of the list are the upper base
forecasts, in the order given by the rows of A. The elements from
n_upper+1 until the end of the list are the bottom base forecasts, in
the order given by the columns of A.

The i-th element depends on the values of `in_type[[i]]` and
`distr[[i]]`.

If `in_type[[i]]`='samples', then `base_fc[[i]]` is a vector containing
samples from the base forecast distribution.

If `in_type[[i]]`='params', then `base_fc[[i]]` is a list containing the
estimated:

- mean and sd for the Gaussian base forecast if `distr[[i]]`='gaussian',
  see [Normal](https://rdrr.io/r/stats/Normal.html);

- lambda for the Poisson base forecast if `distr[[i]]`='poisson', see
  [Poisson](https://rdrr.io/r/stats/Poisson.html);

- size and prob (or mu) for the negative binomial base forecast if
  `distr[[i]]`='nbinom', see
  [NegBinomial](https://rdrr.io/r/stats/NegBinomial.html).

See the description of the parameters `in_type` and `distr` for more
details.

Warnings are triggered from the Importance Sampling step if:

- weights are all zeros, then the upper is ignored during
  reconciliation;

- the effective sample size is \< 200;

- the effective sample size is \< 1% of the sample size (`num_samples`
  if `in_type` is 'params' or the size of the base forecast if if
  `in_type` is 'samples').

Note that warnings are an indication that the base forecasts might have
issues. Please check the base forecasts in case of warnings.

## References

Zambon, L., Azzimonti, D. & Corani, G. (2024). *Efficient probabilistic
reconciliation of forecasts for real-valued and count time series*.
Statistics and Computing 34 (1), 21.
[doi:10.1007/s11222-023-10343-y](https://doi.org/10.1007/s11222-023-10343-y)
.

## See also

[`reconc_gaussian()`](https://idsia.github.io/bayesRecon/reference/reconc_gaussian.md)

## Examples

``` r
library(bayesRecon)

# Create a minimal hierarchy with 2 bottom and 1 upper variable
rec_mat <- get_reconc_matrices(agg_levels = c(1, 2), h = 2)
A <- rec_mat$A
S <- rec_mat$S


# 1) Gaussian base forecasts

# Set the parameters of the Gaussian base forecast distributions
mu1 <- 2
mu2 <- 4
muY <- 9
mus <- c(muY, mu1, mu2)

sigma1 <- 2
sigma2 <- 2
sigmaY <- 3
sigmas <- c(sigmaY, sigma1, sigma2)

base_fc <- list()
for (i in 1:length(mus)) {
  base_fc[[i]] <- list(mean = mus[[i]], sd = sigmas[[i]])
}


# Sample from the reconciled forecast distribution using the BUIS algorithm
buis <- reconc_BUIS(A, base_fc,
  in_type = "params",
  distr = "gaussian", num_samples = 100000, seed = 42
)

samples_buis <- rbind(buis$upper_rec_samples, buis$bottom_rec_samples)

# In the Gaussian case, the reconciled distribution is still Gaussian and can be
# computed in closed form
Sigma <- diag(sigmas^2) # transform into covariance matrix
analytic_rec <- reconc_gaussian(A,
  base_fc_mean = mus,
  base_fc_cov = Sigma
)

# Compare the reconciled means obtained analytically and via BUIS
print(c(S %*% analytic_rec$bottom_rec_mean))
#> [1] 7.411765 2.705882 4.705882
print(rowMeans(samples_buis))
#> [1] 7.413147 2.707427 4.705720


# 2) Poisson base forecasts

# Set the parameters of the Poisson base forecast distributions
lambda1 <- 2
lambda2 <- 4
lambdaY <- 9
lambdas <- c(lambdaY, lambda1, lambda2)

base_fc <- list()
for (i in 1:length(lambdas)) {
  base_fc[[i]] <- list(lambda = lambdas[i])
}

# Sample from the reconciled forecast distribution using the BUIS algorithm
buis <- reconc_BUIS(A, base_fc,
  in_type = "params",
  distr = "poisson", num_samples = 100000, seed = 42
)
samples_buis <- rbind(buis$upper_rec_samples, buis$bottom_rec_samples)

# Print the reconciled means
print(rowMeans(samples_buis))
#> [1] 7.09171 2.36542 4.72629
```
