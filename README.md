
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayesRecon

<!-- badges: start -->
<!-- badges: end -->

bayesRecon is a package for the probabilistic reconciliation of
hierarchical time series forecasts.

## Installation

You can install the development version of bayesRecon like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
# This will be install though our github/gitlab
```

## Example

Let us consider the minimal temporal hierarchy in the figure, where the
bottom variables are the two 6-monthly forecasts and the upper variable
is the yearly forecast.

<img src="./man/figures/minimal_hierarchy.png" width="50%" style="display: block; margin: auto;" />

The hierarchy is described by the *aggregating matrix* S, which can be
obtained using the function `get_reconc_matrices`.

``` r
library(bayesRecon)

rec_mat <- get_reconc_matrices(aggf=c(1,2), h=2)
S <- rec_mat$S
print(S)
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    1    0
#> [3,]    0    1
```

### Poisson base forecasts

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

Nevertheless, we also provide a function for sampling using Markov Chain
Monte Carlo.

``` r
mcmc = reconc_MCMC(S,base_forecasts,distr="poisson",
                   num_samples=30000, seed=42)

samples_mcmc <- mcmc$reconciled_samples
```

- Plots?

## References

Corani, G., Azzimonti, D., Augusto, J.P.S.C., Zaffalon, M. (2021).
*Probabilistic Reconciliation of Hierarchical Forecast via Bayes’ Rule*.
In: Hutter, F., Kersting, K., Lijffijt, J., Valera, I. (eds) Machine
Learning and Knowledge Discovery in Databases. ECML PKDD 2020. Lecture
Notes in Computer Science(), vol 12459. Springer, Cham.
[DOI:10.1007/978-3-030-67664-3_13](https://doi.org/10.1007/978-3-030-67664-3_13).

Corani, G., Rubattu, N., Azzimonti, D., Antonucci, A. (2022).
*Probabilistic reconciliation of count time series*.
[arXiv.2207.09322](https://doi.org/10.48550/arXiv.2207.09322).

Zambon, L., Azzimonti, D. & Corani, G. (2022). *Efficient probabilistic
reconciliation of forecasts for real-valued and count time series*.
[arXiv.2210.02286](https://doi.org/10.48550/arXiv.2210.02286).

Zambon, L., Agosto, A., Giudici, P., Corani, G. (2023). *Properties of
the reconciled distributions for Gaussian and count forecasts*.
[arXiv.2303.15135](https://doi.org/10.48550/arXiv.2303.15135).
