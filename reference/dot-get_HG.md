# Extract and sort hierarchy rows

Function that extracts the "best hierarchy rows" from A (see
.get_hier_rows), and sorts them in the correct order (i.e. bottom-up).
Also sorts accordingly the vectors v, d, it (e.g. of parameters).

## Usage

``` r
.get_HG(A, v, d, it)
```

## Arguments

- A:

  Aggregation matrix

- v:

  Vector of values

- d:

  Distribution parameters

- it:

  Input type parameters

## Value

List containing sorted hierarchy rows and remaining rows
