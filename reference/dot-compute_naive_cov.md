# Estimate covariance from naive or seasonal-naive residuals

Estimates via shrinkage the covariance matrix of the residuals of the
naive or seasonal naive forecasts. If the frequency of the time series
is \> 1, the function chooses between the two methods series-by-series
according to a selection criterion. If the frequency is 1 or is not
provided, only the naive residuals are used.

## Usage

``` r
.compute_naive_cov(y_train, freq = NULL, criterion = "RSS")
```

## Arguments

- y_train:

  Multivariate time series object or numeric matrix of historical
  observations, with dimensions `T x n` (rows are time points, columns
  are series).

- freq:

  Positive integer seasonal frequency (optional). If not provided and
  `y_train` is a multivariate time series, the frequency of the data is
  used.

- criterion:

  Character string used when `freq > 1` to choose residuals. Supported
  values are `"RSS"` (default) and `"seas-test"`. `"RSS"` chooses the
  method with the lower residual sum of squares (RSS), while
  `"seas-test"` uses a statistical test for seasonality (requires
  `forecast` package).

## Value

A numeric `n x n` shrinkage covariance matrix estimated with
[`schaferStrimmer_cov()`](https://idsia.github.io/bayesRecon/reference/schaferStrimmer_cov.md).

## See also

[`schaferStrimmer_cov()`](https://idsia.github.io/bayesRecon/reference/schaferStrimmer_cov.md),
[`reconc_t()`](https://idsia.github.io/bayesRecon/reference/reconc_t.md)
