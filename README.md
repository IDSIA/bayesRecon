
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayesReco

<!-- badges: start -->
<!-- badges: end -->

The goal of bayesReco is to …

## Installation

You can install the development version of bayesReco like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
# This will be install though our github/gitlab
```

## Example

Let us consider the minimal temporal hierarchy in the figure, where the
bottom variables are the two 6-monthly forecasts and the upper variable
is the yearly forecast.

**Insert figure**

The hierarchy is described by the *aggregating matrix* S, which can be
obtained using the function ADD LINK TO get_reconc_matrices FUNCTION

``` r
library(bayesReco)

rec_mat <- get_reconc_matrices(aggf=c(1,2), bottom.f=2, bottom.h=2)
S <- rec_mat$S
print(S)
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    1    0
#> [3,]    0    1
```

We assume that the base forecasts are Poisson distributed, with
parameters given by $\lambda_{Y} = 9$, $\lambda_{S_1} = 2$, and
$\lambda_{S_2} = 4$.

``` r
lambda1 <- 2
lambda2 <- 4
lambdaY <- 9
lambdas <- c(lambdaY,lambda1,lambda2)

base_forecasts = list()
for (i in 1:nrow(S)) {
  base_forecasts[[i]] = lambdas[i]
}
```

We perform probabilistic reconciliation based on conditioning, as
explained in ADD REFERENCES

We recommend using the BUIS algorithm to sample from the reconciled
distribution. (REF?)

``` r
buis <- reconc_BUIS(S, base_forecasts, in_type="params",
                   distr="poisson", num_samples=100000, seed=42)

samples_buis <- buis$reconciled_samples
```

Nevertheless, we also provide a function for sampling using
Metropolis-Hastings.

``` r
mcmc = reconc_MCMC(S,base_forecasts,distr="poisson",
                   num_samples=30000, seed=42)

samples_mcmc <- mcmc$reconciled_samples
```

-   Plots?
