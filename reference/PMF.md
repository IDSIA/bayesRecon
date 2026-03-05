# PMF operations and utilities

A set of functions for working with Probability Mass Functions (PMFs) in
bayesRecon. A PMF is represented as a normalized numeric vector where
element `v[j+1]` represents the probability of value `j` (support starts
from 0).

These functions provide utilities for:

- Drawing samples from PMFs

- Computing summary statistics (mean, variance, quantiles)

- Summarizing PMF distributions

## Usage

``` r
PMF_sample(pmf, N_samples)
PMF_get_mean(pmf)
PMF_get_var(pmf)
PMF_get_quantile(pmf, p)
PMF_summary(pmf, Ltoll = .TOLL, Rtoll = .RTOLL)

PMF_get_mean(pmf)

PMF_get_var(pmf)

PMF_get_quantile(pmf, p)

PMF_summary(pmf, Ltoll = .TOLL, Rtoll = .RTOLL)
```

## Arguments

- pmf:

  A PMF object (numeric vector where element j+1 is the probability of
  value j).

- N_samples:

  Number of samples to draw from the PMF.

- p:

  Probability level for quantile computation (between 0 and 1).

- Ltoll:

  Lower tolerance for computing the minimum of the PMF (default: 1e-15).

- Rtoll:

  Upper tolerance for computing the maximum of the PMF (default: 1e-9).

## Functions

- `PMF_sample(pmf, N_samples)`:

  Draws random samples from the probability distribution specified by
  the PMF. Uses sampling with replacement from the discrete support
  values, weighted by their probabilities. This is useful for generating
  synthetic data or for Monte Carlo simulations.

- `PMF_get_mean(pmf)`:

  Computes the expected value (mean) of the distribution represented by
  the PMF. The mean is calculated as the sum of each value in the
  support weighted by its probability: \\\sum\_{x} x \cdot P(X=x)\\.

- `PMF_get_var(pmf)`:

  Computes the variance of the distribution represented by the PMF. The
  variance measures the spread of the distribution and is calculated as
  \\E\[X^2\] - (E\[X\])^2\\, where \\E\[X\]\\ is the mean.

- `PMF_get_quantile(pmf, p)`:

  Computes the quantile of the distribution at probability level `p`.
  Returns the smallest value `x` such that the cumulative probability up
  to `x` is greater than or equal to `p`. For example, `p=0.5` gives the
  median.

- `PMF_summary(pmf, Ltoll, Rtoll)`:

  Provides a comprehensive summary of the distribution including
  minimum, maximum, quartiles, median, and mean. The minimum and maximum
  are determined based on probability thresholds to handle the
  potentially infinite tails of discrete distributions.

## Examples

``` r
library(bayesRecon)

# Let's build the pmf of a Binomial distribution with parameters n and p
n <- 10
p <- 0.6
pmf_binomial <- apply(matrix(seq(0, n)), MARGIN = 1, FUN = \(x) dbinom(x, size = n, prob = p))

# Draw samples from the PMF object
set.seed(1)
samples <- PMF_sample(pmf = pmf_binomial, N_samples = 1e4)

# Compute statistics
PMF_get_mean(pmf_binomial) # Mean: should be close to n*p = 6
#>      [,1]
#> [1,]    6
PMF_get_var(pmf_binomial) # Variance: should be close to n*p*(1-p) = 2.4
#>      [,1]
#> [1,]  2.4
PMF_get_quantile(pmf_binomial, 0.5) # Median
#> [1] 6
PMF_summary(pmf_binomial) # Full summary
#>   Min. 1st Qu. Median Mean 3rd Qu. Max
#> 1    0       5      6    6       7  10
```
