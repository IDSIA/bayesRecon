# Core Reconciliation via Top-Down Conditioning for Mixed Hierarchies

Internal function that performs the core reconciliation logic using
top-down conditioning (TD) for mixed hierarchies. First, upper forecasts
are reconciled analytically via conditioning (if necessary), then bottom
distributions are updated through probabilistic top-down procedure by
conditioning on the reconciled upper values.

## Usage

``` r
.core_reconc_TDcond(
  A,
  mean_upper,
  cov_upper,
  L_pmf,
  num_samples,
  return_type,
  return_upper = TRUE,
  suppress_warnings = TRUE,
  min_fraction_samples_ok = .MIN_FRACTION_SAMPLES_OK
)
```

## Arguments

- A:

  Matrix (n_upper x n_bottom) defining the hierarchy where upper = A
  %\*% bottom.

- mean_upper:

  Vector of mean forecasts for upper level.

- cov_upper:

  Covariance matrix of upper level forecasts.

- L_pmf:

  List of PMF objects representing the bottom level base forecasts.

- num_samples:

  Number of samples to draw from the reconciled distribution.

- return_type:

  Character string specifying return format: 'pmf', 'samples', or 'all'.

- return_upper:

  Logical, whether to return reconciled upper forecasts (default TRUE).

- suppress_warnings:

  Logical, whether to suppress warnings about samples outside support
  (default TRUE).

- min_fraction_samples_ok:

  Numeric between 0 and 1, minimum fraction of reconciled upper samples
  that must lie in the support of the bottom-up distribution (default
  0.5). If the fraction is below this threshold, the function returns an
  error.

## Value

A list containing:

- `bottom_rec`: List with reconciled bottom forecasts (pmf and/or
  samples).

- `upper_rec`: (only if `return_upper = TRUE`) List with reconciled
  upper forecasts (pmf and/or samples).

## Details

The function internally:

1.  Identifies the "lowest upper" nodes in the hierarchy.

2.  If all uppers are lowest-uppers, samples directly from the upper
    MVN. Otherwise, analytically reconciles the upper hierarchy and
    samples from the lowest level.

3.  Reconciles bottom distributions by conditioning on the
    sampled/reconciled upper values using the probabilistic top-down
    algorithm.

4.  Discards samples that fall outside the support of the bottom-up
    distribution.
