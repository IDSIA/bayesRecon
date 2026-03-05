# Core Reconciliation via Importance Sampling for Mixed Hierarchies

Internal function that performs the core reconciliation logic using
importance sampling (IS) to reconcile mixed-type hierarchies. The base
bottom forecasts (provided as samples) are reweighted according to their
fit to the upper multivariate Gaussian forecasts.

## Usage

``` r
.core_reconc_MixCond(
  A,
  B,
  mean_upper,
  cov_upper,
  num_samples,
  return_type,
  return_ESS = TRUE,
  return_upper = TRUE,
  suppress_warnings = FALSE
)
```

## Arguments

- A:

  Matrix (n_upper x n_bottom) defining the hierarchy where upper = A
  %\*% bottom.

- B:

  Matrix (n_samples x n_bottom) of bottom base forecast samples to be
  reconciled.

- mean_upper:

  Vector of upper level means.

- cov_upper:

  Covariance matrix of upper level.

- num_samples:

  Number of samples to draw/resample from.

- return_type:

  Character string specifying return format: 'pmf', 'samples', or 'all'.

- return_ESS:

  Logical, whether to return the Effective Sample Size (ESS) from
  importance sampling weights (default TRUE).

- return_upper:

  Logical, whether to return the reconciled parameters for the upper
  variables (default TRUE).

- suppress_warnings:

  Logical. If TRUE, suppresses warnings about sample quality (default
  FALSE).

## Value

A list containing:

- `bottom_rec`: List with reconciled bottom forecasts (pmf and/or
  samples).

- `upper_rec`: (only if `return_upper = TRUE`) List with reconciled
  upper forecasts (pmf and/or samples).

- `ESS`: Effective Sample Size resulting from importance sampling
  reweighting (only if `return_ESS = TRUE`).
