# Build hierarchy matrices

Creates the aggregation and summing matrices for a temporal hierarchy of
time series from a user-selected list of aggregation levels.

## Usage

``` r
get_reconc_matrices(agg_levels, h)
```

## Arguments

- agg_levels:

  user-selected list of aggregation levels.

- h:

  number of steps ahead for the bottom level forecasts.

## Value

A list containing the named elements:

- `A` the aggregation matrix;

- `S` the summing matrix.

## See also

[`temporal_aggregation()`](https://idsia.github.io/bayesRecon/reference/temporal_aggregation.md)

## Examples

``` r
library(bayesRecon)

# Create monthly hierarchy
agg_levels <- c(1, 2, 3, 4, 6, 12)
h <- 12
rec_mat <- get_reconc_matrices(agg_levels, h)
S <- rec_mat$S
A <- rec_mat$A
```
