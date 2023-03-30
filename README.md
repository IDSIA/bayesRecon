
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesReco

<!-- badges: start -->
<!-- badges: end -->

The goal of BayesReco is to …

## Installation

You can install the development version of BayesReco like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
# This will be install though our github/gitlab
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(BayesReco)
## basic example code
```

- Basic example of BUIS as in the python code, should run quick and
  ideally have nice plots.

## NOTES ON README.Rmd vs README.md (TO REMOVE FROM FINAL VERSIONS)

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
