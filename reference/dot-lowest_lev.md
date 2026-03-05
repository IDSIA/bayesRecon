# Find rows corresponding to the lowest level

Identifies the rows of aggregation matrix A that correspond to the
lowest hierarchical level (i.e., the most disaggregated upper variables
that collectively sum to all bottom variables).

## Usage

``` r
.lowest_lev(A)
```

## Arguments

- A:

  Aggregation matrix (must be hierarchical)

## Value

Integer vector of row indices corresponding to the lowest level
