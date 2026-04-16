# Probabilistic forecast reconciliation of mixed hierarchies

`reconc_MixCond()` uses importance sampling to draw samples from the
reconciled forecast distribution, obtained via conditioning, in the case
of a mixed hierarchy.

`reconc_TDcond()` uses a top-down conditioning algorithm: first, upper
base forecasts are reconciled via conditioning using only the
hierarchical constraints between the upper; then, the bottom
distributions are updated via a probabilistic top-down procedure.

## Usage

``` r
reconc_MixCond(
  A,
  base_fc_bottom,
  base_fc_upper,
  bottom_in_type = "pmf",
  distr = NULL,
  num_samples = 20000,
  return_type = "pmf",
  return_upper = TRUE,
  suppress_warnings = FALSE,
  seed = NULL
)

reconc_TDcond(
  A,
  base_fc_bottom,
  base_fc_upper,
  bottom_in_type = "pmf",
  distr = NULL,
  num_samples = 20000,
  return_type = "pmf",
  return_upper = TRUE,
  suppress_warnings = TRUE,
  seed = NULL
)
```

## Arguments

- A:

  Aggregation matrix (n_upper x n_bottom).

- base_fc_bottom:

  A list containing the bottom base forecasts, see details.

- base_fc_upper:

  A list containing the upper base forecasts, see details.

- bottom_in_type:

  A string with three possible values:

  - 'pmf' if the bottom base forecasts are in the form of pmf, see
    details;

  - 'samples' if the bottom base forecasts are in the form of samples;

  - 'params' if the bottom base forecasts are in the form of estimated
    parameters.

- distr:

  A string describing the type of bottom base forecasts ('poisson' or
  'nbinom').

  This is only used if `bottom_in_type='params'`.

- num_samples:

  Number of samples drawn from the reconciled distribution. This is
  ignored if `bottom_in_type='samples'`; in this case, the number of
  reconciled samples is equal to the number of samples of the base
  forecasts.

- return_type:

  The return type of the reconciled distributions. A string with three
  possible values:

  - 'pmf' returns a list containing the reconciled marginal pmf objects;

  - 'samples' returns a list containing the reconciled multivariate
    samples;

  - 'all' returns a list with both pmf objects and samples.

- return_upper:

  Logical, whether to return the reconciled parameters for the upper
  variables (default is TRUE).

- suppress_warnings:

  Logical. If `TRUE`, no warnings about samples are triggered; if
  `FALSE`, warnings are generated. Default is `FALSE` for
  `reconc_MixCond` and `TRUE` for `reconc_TDcond`. See the respective
  sections above.

- seed:

  Seed for reproducibility.

## Value

A list containing the reconciled forecasts. The list has the following
named elements:

- `bottom_rec_pmf`: a list of PMF objects for each bottom series (only
  if `return_type` is `'pmf'` or `'all'`);

- `bottom_rec_samples`: a matrix (n_bottom x `num_samples`) of
  reconciled bottom samples (only if `return_type` is `'samples'` or
  `'all'`);

- `upper_rec_pmf`: a list of PMF objects for each upper series (only if
  `return_type` is `'pmf'` or `'all'`, and `return_upper = TRUE`);

- `upper_rec_samples`: a matrix (n_upper x `num_samples`) of reconciled
  upper samples (only if `return_type` is `'samples'` or `'all'`, and
  `return_upper = TRUE`).

## Details

The base bottom forecasts `base_fc_bottom` must be a list of length
n_bottom, where each element is either

- a PMF object (see details below), if `bottom_in_type='pmf'`;

- a vector of samples, if `bottom_in_type='samples'`;

- a list of parameters, if `bottom_in_type='params'`:

  - lambda for the Poisson base forecast if `distr`='poisson', see
    [Poisson](https://rdrr.io/r/stats/Poisson.html);

  - size and prob (or mu) for the negative binomial base forecast if
    `distr`='nbinom', see
    [NegBinomial](https://rdrr.io/r/stats/NegBinomial.html).

The base upper forecasts `base_fc_upper` must be a list containing the
parameters of the multivariate Gaussian distribution of the upper
forecasts. The list must contain only the named elements `mean` (vector
of length n_upper) and `cov` (n_upper x n_upper matrix).

The order of the upper and bottom base forecasts must match the order of
(respectively) the rows and the columns of A.

A PMF object is a numerical vector containing the probability mass
function of a discrete distribution. Each element corresponds to the
probability of the integers from 0 to the last value of the support. See
also [PMF](https://idsia.github.io/bayesRecon/reference/PMF.md) for
functions that handle PMF objects.

**Warnings and errors.**

In `reconc_MixCond`, warnings are triggered from the importance sampling
step if:

- weights are all zeros, then the upper forecast is ignored during
  reconciliation;

- the effective sample size is \< 200;

- the effective sample size is \< 1% of the sample size.

These warnings are an indication that the base forecasts might have
issues. Please check the base forecasts in case of warnings.

In `reconc_TDcond`, if some of the reconciled upper samples lie outside
the support of the bottom-up distribution, those samples are discarded;
the remaining ones are resampled with replacement, so that the number of
output samples is equal to `num_samples`. In this case, a warning is
issued if `suppress_warnings=FALSE` (default is `TRUE`). If the fraction
of discarded samples is above 50%, the function returns an error.

## References

Zambon, L., Azzimonti, D., Rubattu, N., Corani, G. (2024).
*Probabilistic reconciliation of mixed-type hierarchical time series*.
Proceedings of the Fortieth Conference on Uncertainty in Artificial
Intelligence, PMLR 244:4078-4095.
<https://proceedings.mlr.press/v244/zambon24a.html>.

## See also

[`reconc_BUIS()`](https://idsia.github.io/bayesRecon/reference/reconc_BUIS.md),
[`reconc_gaussian()`](https://idsia.github.io/bayesRecon/reference/reconc_gaussian.md),
[PMF](https://idsia.github.io/bayesRecon/reference/PMF.md)

## Examples

``` r
library(bayesRecon)

# Consider a simple hierarchy with two bottom and one upper
A <- matrix(c(1, 1), nrow = 1)
# The bottom forecasts are Poisson with lambda=15
lambda <- 15
n_tot <- 60
base_fc_bottom <- list()
base_fc_bottom[[1]] <- apply(matrix(seq(0, n_tot)), MARGIN = 1,
                             FUN = \(x) dpois(x, lambda = lambda))
base_fc_bottom[[2]] <- apply(matrix(seq(0, n_tot)), MARGIN = 1,
                             FUN = \(x) dpois(x, lambda = lambda))

# The upper forecast is a Normal with mean 40 and std 5
base_fc_upper <- list(mean = 40, cov = matrix(5^2))

# Reconcile with reconc_MixCond
res.mixCond <- reconc_MixCond(A, base_fc_bottom, base_fc_upper)

# Note that the bottom distributions are slightly shifted to the right
PMF_summary(res.mixCond$bottom_rec_pmf[[1]])
#>   Min. 1st Qu. Median     Mean 3rd Qu. Max
#> 1    0      15     18 17.76425      20  31
PMF_summary(base_fc_bottom[[1]])
#>   Min. 1st Qu. Median Mean 3rd Qu. Max
#> 1    0      12     15   15      18  43

PMF_summary(res.mixCond$bottom_rec_pmf[[2]])
#>   Min. 1st Qu. Median     Mean 3rd Qu. Max
#> 1    0      15     18 17.74365      20  31
PMF_summary(base_fc_bottom[[2]])
#>   Min. 1st Qu. Median Mean 3rd Qu. Max
#> 1    0      12     15   15      18  43

# The upper distribution is slightly shifted to the left
PMF_summary(res.mixCond$upper_rec_pmf[[1]])
#>   Min. 1st Qu. Median    Mean 3rd Qu. Max
#> 1    0      33     36 35.5079      38  54
PMF_get_var(res.mixCond$upper_rec_pmf[[1]])
#>          [,1]
#> [1,] 14.58635


library(bayesRecon)

# Consider a simple hierarchy with two bottom and one upper
A <- matrix(c(1, 1), nrow = 1)
# The bottom forecasts are Poisson with lambda=15
lambda <- 15
n_tot <- 60
base_fc_bottom <- list()
base_fc_bottom[[1]] <- apply(matrix(seq(0, n_tot)), MARGIN = 1,
                             FUN = \(x) dpois(x, lambda = lambda))
base_fc_bottom[[2]] <- apply(matrix(seq(0, n_tot)), MARGIN = 1,
                             FUN = \(x) dpois(x, lambda = lambda))

# The upper forecast is a Normal with mean 40 and std 5
base_fc_upper <- list(mean = 40, cov = matrix(c(5^2)))

# Reconcile with reconc_TDcond
res.TDcond <- reconc_TDcond(A, base_fc_bottom, base_fc_upper)

# Note that the bottom distributions are shifted to the right
PMF_summary(res.TDcond$bottom_rec_pmf[[1]])
#>   Min. 1st Qu. Median     Mean 3rd Qu. Max
#> 1    0      17     20 19.97955      23  37
PMF_summary(base_fc_bottom[[1]])
#>   Min. 1st Qu. Median Mean 3rd Qu. Max
#> 1    0      12     15   15      18  43

PMF_summary(res.TDcond$bottom_rec_pmf[[2]])
#>   Min. 1st Qu. Median     Mean 3rd Qu. Max
#> 1    0      17     20 19.98025      23  37
PMF_summary(base_fc_bottom[[2]])
#>   Min. 1st Qu. Median Mean 3rd Qu. Max
#> 1    0      12     15   15      18  43

# The upper distribution remains similar
PMF_summary(res.TDcond$upper_rec_pmf[[1]])
#>   Min. 1st Qu. Median    Mean 3rd Qu. Max
#> 1    0      37     40 39.9598      43  60
PMF_get_var(res.TDcond$upper_rec_pmf[[1]])
#>          [,1]
#> [1,] 25.58181

## Example 2: reconciliation with unbalanced hierarchy
# We consider the example in Fig. 9 of Zambon et al. (2024).

# The hierarchy has 5 bottoms and 3 uppers
A <- matrix(c(
  1, 1, 1, 1, 1,
  1, 1, 0, 0, 0,
  0, 0, 1, 1, 0
), nrow = 3, byrow = TRUE)
# Note that the 5th bottom only appears in the highest level, this is an unbalanced hierarchy.
n_upper <- nrow(A)
n_bottom <- ncol(A)

# The bottom forecasts are Poisson with lambda=15
lambda <- 15
n_tot <- 60
base_fc_bottom <- list()
for (i in seq(n_bottom)) {
  base_fc_bottom[[i]] <- apply(matrix(seq(0, n_tot)), MARGIN = 1,
                               FUN = \(x) dpois(x, lambda = lambda))
}

# The upper forecasts are a multivariate Gaussian
mean <- c(75, 30, 30)
cov <- matrix(c(
  5^2, 5, 5,
  5, 10, 0,
  5, 0, 10
), nrow = 3, byrow = TRUE)

base_fc_upper <- list(mean = mean, cov = cov)
if (FALSE) { # \dontrun{
# If we reconcile with reconc_TDcond it won't work (unbalanced hierarchy)
res.TDcond <- reconc_TDcond(A, base_fc_bottom, base_fc_upper)
} # }

# We can balance the hierarchy by duplicating the node b5:
# i) consider the time series observations for b5 as the upper u4,
# ii) fit the multivariate ts model for u1, u2, u3, u4.

# In this example we simply assume that the forecast for u1-u4 is
# Gaussian with the mean and variance of u4 given by the parameters in b5.
mean_b5 <- lambda
var_b5 <- lambda
mean <- c(75, 30, 30, mean_b5)
cov <- matrix(c(
  5^2, 5, 5, 5,
  5, 10, 0, 0,
  5, 0, 10, 0,
  5, 0, 0, var_b5
), nrow = 4, byrow = TRUE)
base_fc_upper <- list(mean = mean, cov = cov)

# We also need to update the aggregation matrix
A <- matrix(c(
  1, 1, 1, 1, 1,
  1, 1, 0, 0, 0,
  0, 0, 1, 1, 0,
  0, 0, 0, 0, 1
), nrow = 4, byrow = TRUE)

# We can now reconcile with TDcond
res.TDcond <- reconc_TDcond(A, base_fc_bottom, base_fc_upper)

# Note that the reconciled distribution of b5 and u4 are identical,
# keep this in mind when using the results of your reconciliation!
max(abs(res.TDcond$bottom_rec_pmf[[5]] - res.TDcond$upper_rec_pmf[[4]]))
#> [1] 0
```
