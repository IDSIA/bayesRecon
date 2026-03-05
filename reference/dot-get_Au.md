# Get aggregation matrix of upper sub-hierarchy

Constructs the aggregation matrix Au for the sub-hierarchy composed only
of the upper variables. This matrix relates the "upper uppers" to the
"lowest uppers" (the most disaggregated level of the upper hierarchy).

## Usage

``` r
.get_Au(A, lowest_rows = NULL)
```

## Arguments

- A:

  Aggregation matrix (must be hierarchical)

- lowest_rows:

  Optional. Integer vector of row indices corresponding to the lowest
  upper level. If NULL, computed using `.lowest_lev(A)`.

## Value

Aggregation matrix for the upper sub-hierarchy, or NULL if all uppers
are at the lowest level.
