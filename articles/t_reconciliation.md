# t-Rec: Reconciliation of Gaussian forecasts by modelling the uncertainty on the covariance matrix

## Introduction

This vignette showcases the method *t-Rec*, a Bayesian approach to
forecast reconciliation that explicitly accounts for uncertainty in the
covariance matrix of residuals.

By exploiting the flexibility of Bayesian models to incorporate
parameter uncertainty, t-Rec assumes an Inverse-Wishart prior on the
covariance matrix. This choice allows the reconciliation to be derived
in closed form, resulting in a reconciled predictive distribution that
follows a multivariate t. See Carrara et al. (2025) for a detailed
explanation.

``` r
# load the packages
library(bayesRecon)
library(forecast) # base forecasts
library(ggplot2)  # plots
```

## Loading the data

We consider the monthly Swiss overnight tourist stays data divided by
Canton and aggregated over the whole country. This is a cross-sectional
hierarchical structure, consisting of two levels: the national total and
the division by canton.

The data is available in this package and can be loaded as
[`bayesRecon::swiss_tourism`](https://idsia.github.io/bayesRecon/reference/swiss_tourism.md).
This is a list that contains the `mts` object, the aggregation matrix
and the number of upper and bottom time series. The raw data is also
publicly available on the official website of the Swiss Confederation at
this [link](https://www.bfs.admin.ch/asset/en/px-x-1003020000_102).

The data set spans the period from January 2005 to January 2025,
comprising 241 monthly observations.

The time series exhibit strong seasonality. See, in **Figure 1**, the
top-level (aggregate) time series.

``` r
# save all time series 
Y = swiss_tourism$ts

# plot the first (top) time series
autoplot(Y[,1], ylab = "Overnight stays in Switzerland",linewidth=0.9)+
  scale_y_continuous(labels = function(x) paste0(formatC(x / 1e6, format = "g"), "M"))
```

![\*\*Figure 1\*\*: Swiss tourism - monthly Swiss overnight
stays.](t_reconciliation_files/figure-html/swiss-tourism-plot-1.png)

**Figure 1**: Swiss tourism - monthly Swiss overnight stays.

We load the aggregation matrix and select the length of the training
set. The length of the training set can be selected within the range
\[26, 240\] observations.

``` r
# Save aggregation matrix
A = swiss_tourism$agg_mat

# Number of bottom and upper time series
n_b = ncol(A)
n_u = nrow(A)
n = n_b + n_u

# Frequency is monthly:
print(frequency(Y))
#> [1] 12

# Select the length of the training set and the forecast horizon
L = 60    
h = 1    

# Select the training set and the actuals for the forecast horizon
train = window(Y, end = time(Y)[L])
actuals = window(Y, start= time(Y)[L + 1], end = time(Y)[L + h])
```

## Forecasts

The base forecasts are generated using the
[`ets()`](https://pkg.robjhyndman.com/forecast/reference/ets.html)
function from the `forecast` package, which fits individual exponential
smoothing models to each time series.

For each time series forecast, we save the base forecast mean
(`base_fc`) and the residuals from the fitted model (`res`), which are
used to estimate the covariance matrix of the forecast errors.

``` r
# Compute base forecasts and residuals for each time series
base_fc = rep(NA, n)
res = matrix(NA, ncol = n, nrow = L)
for (i in 1:n){
    models = forecast::ets(train[,i], model = "AZZ")
    f = forecast(models, h = h)
    base_fc[i] = f$mean
    res[, i] = models$residuals
}
```

## Reconciliation

We can now reconcile the forecasts with the function
[`reconc_t()`](https://idsia.github.io/bayesRecon/reference/reconc_t.md),
which implements the t-Rec method.

This function takes as input the aggregation matrix `A`, the base
forecast means `base_fc_mean`, the training data `y_train` and the model
residuals `residuals` and returns the parameters of the reconciled
forecasts, which are distributed as a multivariate-t. The flag
`return_upper = TRUE` forces the function to return also the parameters
of the upper-level reconciled forecasts. The flag
`return_parameters = TRUE` forces the function to return also the
parameters of the posterior distribution of the covariance matrix, which
we will use to compare the covariance estimates obtained with t-Rec and
MinT.

``` r
t_rec_results = reconc_t(A, base_fc_mean = base_fc, 
                         y_train = train, 
                         residuals = res,
                         return_parameters = TRUE,
                         return_upper = TRUE)
```

For the selected data window, we save the base forecasts (Base).

``` r
# Base forecasts
Base_mean = base_fc
Base_cov_mat = crossprod(res)/nrow(res) # covariance of the residuals
```

We further compute the reconciled forecasts with the standard Gaussian
reconciliation (MinT). Note the flag `return_upper = TRUE`, which forces
the function to return the parameters of the upper-level reconciled
forecasts.

``` r
# Gaussian/MinT: compute reconciliation with bayesRecon
gauss_results = reconc_gaussian(A, base_fc, residuals = res, return_upper = TRUE)
# Reconciled mean for the whole hierarchy:
MinT_reconciled_mean = c(gauss_results$upper_rec_mean, 
                         gauss_results$bottom_rec_mean)
```

## Comparison of results

We now compare the predictive densities for the 1-step ahead forecasts
of the upper-level time series (Switzerland) obtained by the three
methods above: the base forecasts (Base), Minimum Trace reconciliation
(MinT), and t-Rec.

![\*\*Figure 2\*\*: Predictive densities of the upper time series
obtained with MinT (purple), t-Rec (green) and Base (blue). The black
triangle indicates the actual
value.](t_reconciliation_files/figure-html/comparison-plot-1.png)

**Figure 2**: Predictive densities of the upper time series obtained
with MinT (purple), t-Rec (green) and Base (blue). The black triangle
indicates the actual value.

The base forecast and the MinT reconciled forecast are normally
distributed, while t-Rec has a Student’s t distribution. Both
reconciliation methods improve the forecast: the mean of the densities
is closer to the actual value. However, the t-Rec method returns a wider
density, which is more likely to contain the actual value. This is
because t-Rec accounts for the uncertainty in the covariance matrix of
the forecast errors, which leads to a more realistic representation of
the forecast uncertainty.

## Comparison of the covariance matrix estimates

In this section we showcase the main strength of t-Rec: it provides an
estimate of the covariance which accounts for the uncertainty in the
estimation. To illustrate this point, we compare the covariance matrix
estimates obtained with t-Rec and with the standard Gaussian
reconciliation method (MinT).

We focus on the covariance matrix between the upper-level series,
denoted as CH, and the bottom-level time series with the largest average
values “Graubünden”, denoted as GR. Analogous considerations apply also
for other series. The code in the vignette is parametric so that a user
could manually visualize different comparisons.

While MinT only returns a point estimate for the covariance matrix,
t-Rec provides a distribution for the covariance matrix. This is
obtained in a Bayesian way: starting from a prior on the covariance, we
include the data and obtain a posterior distribution which is an Inverse
Wishart with known parameters $\nu$ and $\Psi$. Those parameters are
stored in `t_rec_results$posterior_nu` and
`t_rec_results$posterior_Psi`.

Since we know analytically the distribution of the covariance, we can
also compute the marginal distribution of the variances in closed-form.
This is a scaled inverse Gamma distribution with the following
parameters:
$$\Sigma_{ii} \sim \text{Inv-Gamma}\left( \frac{\nu - n + 1}{2},\frac{\Psi_{ii}}{2} \right).$$

Instead of visualizing the variance, we show the standard deviation
which produces a more intuitive interpretation. The standard deviation
is obtained by applying the square root transformation to the variance,
which results in a non-standard distribution. The density of the
standard deviation can be computed using the change of variable formula,
which leads to the following expression:

``` r
# Select which series to plot
i = 1     # CH
j = 19    # GR

# density of the inverse gamma
dinvgamma <- function(x, shape, rate) {
  dgamma(1/x, shape = shape, rate = rate) / x^2
}

# density of the standard deviation (square root transform of variance)
d_std_dev <- function(x, shape,rate){
  dinvgamma(x^2, shape = shape, rate = rate)*2*x
}

# compute the density of the variance parameters for CH
mean_ch <- t_rec_results$posterior_Psi[i, i]/(t_rec_results$posterior_nu - n - 1)
x_ch <- sqrt(seq(mean_ch*0.5, mean_ch*1.7, length.out = 1000))
shape_ch <- (t_rec_results$posterior_nu - n + 1) /2
rate_ch  <- t_rec_results$posterior_Psi[i, i] / 2
dens_ch <- d_std_dev(x_ch, shape = shape_ch, rate = rate_ch)

# compute the density of the variance parameters for GR
mean_gr <- t_rec_results$posterior_Psi[j, j]/(t_rec_results$posterior_nu - n - 1)
x_gr <- sqrt(seq(mean_gr*0.5, mean_gr*1.7, length.out = 1000))
shape_gr <- (t_rec_results$posterior_nu - n + 1) /2
rate_gr  <- t_rec_results$posterior_Psi[j, j] / 2
dens_gr <- d_std_dev(x_gr, shape = shape_gr, rate = rate_gr)
```

**Figure 3** shows the density of the posterior standard deviation of
the forecasts for the upper time series (CH) and the bottom time series
(GR).

![\*\*Figure 3\*\*: Density of the posterior standard deviation of the
forecasts.](t_reconciliation_files/figure-html/density%20plot-1.png)

**Figure 3**: Density of the posterior standard deviation of the
forecasts.

The posterior distribution for the covariance and for the correlation
values is not available in closed form but it can be obtained via
sampling. Since the posterior distribution is an Inverse-Wishart
distribution, we can sample from this distribution using the custom
function `rinvwishart()`, defined below.

``` r
# generate k samples from an IW(Psi, nu) distribution
rinvwishart <- function(k, nu, Psi, seed=42) {
  p <- nrow(Psi)
  Sigma <- solve(Psi)
  
  set.seed(seed)
  all_W <- rWishart(k, df = nu, Sigma = Sigma)
  
  W <- array(NA, dim = c(p, p, k))
  for (i in 1:k) {
    W[,,i] <- solve(all_W[,,i])
  }
  return(W)
}

IW_post_samples <- rinvwishart(k=1000, nu = t_rec_results$posterior_nu,
                            Psi = t_rec_results$posterior_Psi)
```

**Figure 4** shows the posterior density of the correlation between CH
and GR. The value estimated by MinT, plotted as a vertical dashed line,
differs from the posterior mode estimated by t-Rec, illustrating how
t-Rec captures a different view of the dependence structure.

![\*\*Figure 4\*\*: Density of the posterior correlation between CH and
GR obtained with
t-Rec.](t_reconciliation_files/figure-html/plot%20densities-1.png)

**Figure 4**: Density of the posterior correlation between CH and GR
obtained with t-Rec.

Finally, we show in **Table 1** the standard deviation and correlation
estimates for the upper-level series (CH) and the bottom-level series
(GR) obtained with Base, MinT and t-Rec methods. While Base and MinT
provide only point estimates for the standard deviation and correlation,
t-Rec provides a distribution for these parameters. In the table, we
report the mean of the posterior distribution for the standard deviation
and the correlation.

| Method | ${\widehat{\sigma}}_{\text{CH}}$ | ${\widehat{\sigma}}_{\text{GR}}$ | ${\widehat{\rho}}_{\text{CH,GR}}$ |
|:-------|---------------------------------:|---------------------------------:|----------------------------------:|
| Base   |                           96,618 |                           34,747 |                              0.79 |
| MinT   |                           84,580 |                           34,433 |                              0.79 |
| t-Rec  |                          121,587 |                           39,589 |                              0.73 |

**Table 1**: Standard deviation and correlation estimates for CH and GR.

Note that the standard deviation estimated by t-Rec is higher than both
the base and the MinT estimates. This is because the posterior
Inverse-Wishart distribution depends on the prior, which is initialized
using the residuals of the naive forecasts. Since these residuals
capture the full uncertainty of the non-reconciled forecasts, the prior
inflates the posterior variance relative to the point estimate provided
by MinT. Moreover, the correlation between CH and GR is also estimated
differently by t-Rec compared to MinT.

## References

Carrara, Chiara, Dario Azzimonti, Giorgio Corani, and Lorenzo Zambon.
2025. “Modeling the Uncertainty on the Covariance Matrix for
Probabilistic Forecast Reconciliation.”
<https://arxiv.org/abs/2506.19554>.

Zambon, Lorenzo, Arianna Agosto, Paolo Giudici, and Giorgio Corani.
2024. “Properties of the Reconciled Distributions for Gaussian and Count
Forecasts.” *International Journal of Forecasting* 40 (4): 1438–48.
https://doi.org/<https://doi.org/10.1016/j.ijforecast.2023.12.004>.
