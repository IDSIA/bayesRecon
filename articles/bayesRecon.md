# Probabilistic Reconciliation via Conditioning with \`bayesRecon\`

## Introduction

This vignette shows how to perform *probabilistic reconciliation* with
the `bayesRecon` package. We provide three examples:

1.  *Temporal hierarchy for a count time series*: we build a temporal
    hierarchy over a count time series, produce the base forecasts using
    `glarma` and reconcile them via Bottom-Up Importance Sampling
    (BUIS).

2.  *Temporal hierarchy for a smooth time series*: we build a temporal
    hierarchy over a smooth time series, compute the base forecasts
    using `ets` and we reconcile them in closed form using Gaussian
    reconciliation. The covariance matrix is diagonal.

3.  *Hierarchical of smooth time series*: this is an example of a
    cross-sectional hierarchy. We generate the base forecasts using
    `ets` and we reconcile them via Gaussian reconciliation. The
    covariance matrix is full and estimated via shrinkage.

## Installation

The package, available at this [CRAN
page](https://cran.r-project.org/package=bayesRecon), can be installed
and loaded with the usual commands:

``` r
install.packages("bayesRecon", dependencies = TRUE)
```

Load the package:

``` r
library(bayesRecon)
```

## Temporal hierarchy over a count time series

We select a monthly time series of counts from the *carparts* dataset,
available from the expsmooth package (R. J. Hyndman 2015). The data set
contains time series of sales of cars part from Jan. 1998 to Mar. 2002.
For this example we select time series \#2655, which we make available
as
[`bayesRecon::carparts_example`](https://idsia.github.io/bayesRecon/reference/carparts_example.md).

This time series has a skewed distribution of values.

``` r
layout(mat = matrix(c(1, 2), nrow = 1, ncol = 2), widths = c(2, 1))
plot(carparts_example, xlab = "Time", ylab = "Car part sales", main = NULL)
hist(carparts_example, xlab = "Car part sales", main = NULL)
```

![\*\*Figure 1\*\*: Carpart - monthly car part
sales.](bayesRecon_files/figure-html/carpart-plot-1.png)

**Figure 1**: Carpart - monthly car part sales.

  
  
We divide the time series into train and test; the test set contains the
last 12 months.

``` r
train <- window(carparts_example, end = c(2001, 3))
test <- window(carparts_example, start = c(2001, 4))
```

We build the temporal hierarchy using the `temporal aggregation`
function. We specify the aggregation levels using the `agg_levels`
argument; in this case they are *2-Monthly*, *Quarterly*, *4-Monthly*,
*Biannual*, and *Annual*.

``` r
train.agg <- bayesRecon::temporal_aggregation(train, agg_levels = c(2, 3, 4, 6, 12))
levels <- c("Annual", "Biannual", "4-Monthly", "Quarterly", "2-Monthly", "Monthly")
names(train.agg) <- levels
```

The function returns a list of aggregated time series, ordered from the
most aggregated (top of the hierarchy) to the most disagreggated (bottom
of the hierarchy). We plot them below.

``` r
par(mfrow = c(2, 3), mai = c(0.6, 0.6, 0.5, 0.5))
for (l in levels) {
  plot(train.agg[[l]], xlab = "Time", ylab = "Car part sales", main = l)
}
```

![\*\*Figure 2\*\*: The aggregated time series of the temporal
hierarchy.](bayesRecon_files/figure-html/temp-agg-plot-1.png)

**Figure 2**: The aggregated time series of the temporal hierarchy.

  
  
We compute the *base forecasts* using the package
[`glarma`](https://cran.r-project.org/package=glarma), a package
specific for forecasting count time series. We forecast 12 steps ahead
at the monthly level, 4 steps ahead at the quarterly level, etc. by
iterating over the levels of the hierarchy, At each level, we fit a
`glarma` model with Poisson predictive distribution and we forecast;
each forecast is constituted by 20000 integer samples.

Eventually we collect the samples of the 28 predictive distributions
(one at the *Annual* level, two at the *Biannual* level, etc.) in a
list. The code below takes about 30 seconds on a standard computer.

``` r
# install.packages("glarma", dependencies = TRUE)
# library(glarma)

fc.samples <- list()
D <- 20000
fc.count <- 1

# iterating over the temporal aggregation levels
for (l in seq_along(train.agg)) {
  f.level <- frequency(train.agg[[l]])
  print(paste("Forecasting at ", levels[l], "...", sep = ""))
  # fit an independent model for each aggregation level
  model <- glarma::glarma(
    train.agg[[l]],
    phiLags = if (f.level == 1) 1 else 1:(min(6, f.level - 1)),
    thetaLags = if (f.level == 1) NULL else f.level,
    X = cbind(intercept = rep(1, length(train.agg[[l]]))),
    offset = cbind(intercept = rep(0, length(train.agg[[l]]))),
    type = "Poi"
  )
  # forecast 1 year ahead
  h <- f.level
  tmp <- matrix(data = NA, nrow = h, ncol = D)
  for (s in (1:D)) {
    # each call to 'forecast.glarma' returns a simulation path
    tmp[, s] <- glarma::forecast(
      model,
      n.ahead = h,
      newdata = cbind(intercept = rep(1, h)),
      newoffset = rep(0, h)
    )$Y
  }
  # collect the forecasted samples
  for (i in 1:h) {
    fc.samples[[fc.count]] <- tmp[i, ]
    fc.count <- fc.count + 1
  }
}
#> [1] "Forecasting at Annual..."
#> [1] "Forecasting at Biannual..."
#> [1] "Forecasting at 4-Monthly..."
#> [1] "Forecasting at Quarterly..."
#> [1] "Forecasting at 2-Monthly..."
#> [1] "Forecasting at Monthly..."
```

Reconciliation requires the aggregation matrix $\mathbf{A}$, which we
obtain using the function `get_reconc_matrices`. It requires:

- the aggregation factors of the hierarchy, which in this example are
  $\{ 2,3,4,6,12\}$;
- the length of the forecasting horizon at the bottom level, which is 12
  in this example.

``` r
recon.matrices <- bayesRecon::get_reconc_matrices(agg_levels = c(2, 3, 4, 6, 12), h = 12)
# Aggregation matrix
A <- recon.matrices$A
```

To reconcile using Bottom-Up Important Sampling (BUIS) we we use the
function `reconc_BUIS`, passing to it the $\mathbf{A}$ matrix, the *base
forecasts*, the type of the base forecasts (`in_type`=“samples”) and
whether the samples are continuous or discrete (`distr` = “discrete”).

``` r
recon.res <- bayesRecon::reconc_BUIS(
  A,
  base_fc = fc.samples,
  in_type = "samples",
  distr = "discrete",
  seed = 42
)
```

Here we obtain samples from the reconciled forecast distribution.

``` r
reconciled_samples <- recon.res$bottom_rec_samples
dim(reconciled_samples)
#> [1]    12 20000
```

We now compute the Mean Absolute Error (MAE) and the Continuous Ranked
Probability Score (CRPS) for the bottom (i.e., *monthly*) time series.
For computing CRPS, we use the package
[`scoringRules`](https://cran.r-project.org/package=scoringRules).

``` r
# install.packages("scoringRules", dependencies = TRUE)
library(scoringRules)

ae.fc <- list()
ae.reconc <- list()
crps.fc <- list()
crps.reconc <- list()
for (h in 1:length(test)) {
  y.hat_ <- median(fc.samples[[nrow(A) + h]])
  y.reconc_ <- median(recon.res$bottom_rec_samples[h, ])
  # Compute Absolute Errors
  ae.fc[[h]] <- abs(test[h] - y.hat_)
  ae.reconc[[h]] <- abs(test[h] - y.reconc_)
  # Compute Continuous Ranked Probability Score (CRPS)
  crps.fc[[h]] <-
    scoringRules::crps_sample(y = test[h], dat = fc.samples[[nrow(A) + h]])
  crps.reconc[[h]] <-
    scoringRules::crps_sample(y = test[h], dat = recon.res$bottom_rec_samples[h, ])
}

mae.fc <- mean(unlist(ae.fc))
mae.reconc <- mean(unlist(ae.reconc))
crps.fc <- mean(unlist(crps.fc))
crps.reconc <- mean(unlist(crps.reconc))
metrics <- data.frame(
  row.names = c("MAE", "CRPS"),
  base.forecasts = c(mae.fc, crps.fc),
  reconciled.forecasts = c(mae.reconc, crps.reconc)
)
metrics
#>      base.forecasts reconciled.forecasts
#> MAE       1.2500000            1.0000000
#> CRPS      0.7093284            0.6509691
```

## Temporal hierarchy over a smooth time series

In this second example, we select a smooth monthly time series (N1485)
from the M3 forecasting competition (Makridakis and Hibon 2000). The
data set is available in the Mcomp package (R. Hyndman 2018) and we make
available this time series as
[`bayesRecon::M3_example`](https://idsia.github.io/bayesRecon/reference/M3_example.md).

``` r
plot(M3_example$train, xlab = "Time", ylab = "y", main = "N1485")
```

![\*\*Figure 3\*\*: M3 - N1485 time
series.](bayesRecon_files/figure-html/m3-plot-1.png)

**Figure 3**: M3 - N1485 time series.

  
We build the temporal hierarchy using the `temporal_aggregation`
function.

Without specifying `agg_levels`, the function generates by default all
the feasible aggregation: 2-Monthly, Quarterly, 4-Monthly, Biannual, and
Annual.

``` r
train.agg <- bayesRecon::temporal_aggregation(M3_example$train)
levels <- c("Annual", "Biannual", "4-Monthly", "Quarterly", "2-Monthly", "Monthly")
names(train.agg) <- levels
```

We generate the base forecasts using `ets` from the
[forecast](https://cran.r-project.org/package=forecast) package (R. J.
Hyndman and Khandakar 2008). At every level we predict 18 months ahead.
All the predictive distributions are Gaussian.

``` r
# install.packages("forecast", dependencies = TRUE)
library(forecast)

H <- length(M3_example$test)
H
#> [1] 18

fc <- list()
level.idx <- 1
fc.idx <- 1
for (level in train.agg) {
  level.name <- names(train.agg)[level.idx]
  # fit an ETS model for each temporal level
  model <- ets(level)
  # generate forecasts for each level within 18 months
  h <- floor(H / (12 / frequency(level)))
  print(paste("Forecasting at ", level.name, ", h=", h, "...", sep = ""))
  level.fc <- forecast(model, h = h)
  # save mean and sd of the gaussian predictive distribution
  for (i in 1:h) {
    fc[[fc.idx]] <- list(
      mean = level.fc$mean[[i]],
      sd = (level.fc$upper[, "95%"][[i]] - level.fc$mean[[i]]) / qnorm(0.975)
    )
    fc.idx <- fc.idx + 1
  }
  level.idx <- level.idx + 1
}
#> [1] "Forecasting at Annual, h=1..."
#> [1] "Forecasting at Biannual, h=3..."
#> [1] "Forecasting at 4-Monthly, h=4..."
#> [1] "Forecasting at Quarterly, h=6..."
#> [1] "Forecasting at 2-Monthly, h=9..."
#> [1] "Forecasting at Monthly, h=18..."
```

Using the function `get_reconc_matrices`, we get matrix $\mathbf{A}$.

``` r
rmat <- get_reconc_matrices(agg_levels = c(2, 3, 4, 6, 12), h = 18)

par(mai = c(1, 1, 0.5, 0.5))
image(1:ncol(rmat$A), 1:nrow(rmat$A),
  t(apply(t(rmat$A), 1, rev)),
  xaxt = "n", yaxt = "n", ylab = "", xlab = levels[6]
)
axis(1, at = 1:ncol(rmat$A), label = 1:ncol(rmat$A), las = 1)
axis(2, at = c(23, 22, 19, 15, 9), label = levels[1:5], las = 2)
```

![\*\*Figure 4\*\*: M3 - The aggregation matrix A (red=1,
yellow=0).](bayesRecon_files/figure-html/m3-rmat-1.png)

**Figure 4**: M3 - The aggregation matrix A (red=1, yellow=0).

  
The function `reconc_gaussian` implements the exact Gaussian
reconciliation. We also run `reconc_BUIS`, to check the consistency
between the two approaches.

``` r
recon.gauss <- bayesRecon::reconc_gaussian(
  A = rmat$A,
  base_fc_mean = sapply(fc, "[[", 1),
  base_fc_cov = diag(sapply(fc, "[[", 2))^2
)

reconc.buis <- bayesRecon::reconc_BUIS(
  A = rmat$A,
  base_fc = fc,
  in_type = "params",
  distr = "gaussian",
  num_samples = 20000,
  seed = 42
)

# check that the algorithms return consistent results
round(rbind(
  c(rmat$S %*% recon.gauss$bottom_rec_mean),
  rowMeans(rbind(reconc.buis$upper_rec_samples, reconc.buis$bottom_rec_samples))
))
#>       [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
#> [1,] 74977 35913 39063 41491 23520 25136 26321 27174 17464 18450 19251 19812
#> [2,] 74946 35897 39049 41470 23532 25091 26323 27134 17462 18435 19231 19818
#>      [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
#> [1,] 20412 21079 11527 11993 12393 12743 13047 13274 13531 13643 14317  5694
#> [2,] 20425 21046 11547 11985 12365 12726 13041 13282 13482 13652 14336  5694
#>      [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35] [,36]
#> [1,]  5833  5936  6056  6138  6255  6324  6419  6508  6538  6619  6655  6759
#> [2,]  5853  5914  6070  6127  6237  6308  6418  6505  6536  6612  6670  6729
#>      [,37] [,38] [,39] [,40] [,41]
#> [1,]  6772  6881  6762  7160  7157
#> [2,]  6753  6943  6709  7152  7184
```

We now compare *base forecasts* and *reconciled forecasts*:

``` r
yhat.mu <- tail(sapply(fc, "[[", 1), 18)
yhat.sigma <- tail(sapply(fc, "[[", 2), 18)
yhat.hi95 <- qnorm(0.975, mean = yhat.mu, sd = yhat.sigma)
yhat.lo95 <- qnorm(0.025, mean = yhat.mu, sd = yhat.sigma)
yreconc.mu <- rowMeans(reconc.buis$bottom_rec_samples)
yreconc.hi95 <- apply(
  reconc.buis$bottom_rec_samples, 1,
  function(x) quantile(x, 0.975)
)
yreconc.lo95 <- apply(
  reconc.buis$bottom_rec_samples, 1,
  function(x) quantile(x, 0.025)
)

ylim_ <- c(
  min(M3_example$train, M3_example$test, yhat.lo95, yreconc.lo95) - 1,
  max(M3_example$train, M3_example$test, yhat.hi95, yreconc.hi95) + 1
)

plot(M3_example$train,
  xlim = c(1990, 1995.6), ylim = ylim_,
  ylab = "y", main = "N1485 Forecasts"
)
lines(M3_example$test, lty = "dashed")
lines(ts(yhat.mu, start = start(M3_example$test), frequency = 12),
  col = "coral", lwd = 2
)
lines(ts(yreconc.mu, start = start(M3_example$test), frequency = 12),
  col = "blue2", lwd = 2
)
xtest <- time(M3_example$test)
polygon(c(xtest, rev(xtest)), c(yhat.mu, rev(yhat.hi95)),
  col = "#FF7F5066", border = "#FF7F5066"
)
polygon(c(xtest, rev(xtest)), c(yhat.mu, rev(yhat.lo95)),
  col = "#FF7F5066", border = "#FF7F5066"
)
polygon(c(xtest, rev(xtest)), c(yreconc.mu, rev(yreconc.hi95)),
  col = "#0000EE4D", border = "#0000EE4D"
)
polygon(c(xtest, rev(xtest)), c(yreconc.mu, rev(yreconc.lo95)),
  col = "#0000EE4D", border = "#0000EE4D"
)
```

![\*\*Figure 5\*\*: M3 - Base and reconciled forecasts. The black line
shows the actual data (dashed in the test). The orange line is the mean
of the base forecasts, the blu line is the reconciled mean. We also show
the 95% prediction
intervals.](bayesRecon_files/figure-html/m3-plotfore-1.png)

**Figure 5**: M3 - Base and reconciled forecasts. The black line shows
the actual data (dashed in the test). The orange line is the mean of the
base forecasts, the blu line is the reconciled mean. We also show the
95% prediction intervals.

## Gaussian reconciliation of a cross-sectional hierarchy

In this example, we consider the hierarchical time series *infantgts*,
which is available from the `hts` package (R. Hyndman et al. 2021) We
make it available also in our package as
[`bayesRecon::infantMortality`](https://idsia.github.io/bayesRecon/reference/infantMortality.md).

It contains counts of infant mortality (deaths) in Australia,
disaggregated by state and sex (male and female).

We forecast one year ahead using `auto.arima` from the
[`forecast`](https://cran.r-project.org/package=forecast) package. We
collect the residuals, which we will later use to compute the covariance
matrix.

``` r
# install.packages("forecast", dependencies = TRUE)
library(forecast)

fc <- list()
residuals <- matrix(NA,
  nrow = length(infantMortality$total),
  ncol = length(infantMortality)
)
fc.idx <- 1
for (s in infantMortality) {
  s.name <- names(infantMortality)[fc.idx]
  print(paste("Forecasting at ", s.name, "...", sep = ""))
  # fit an auto.arima model and forecast with h=1
  model <- auto.arima(s)
  s.fc <- forecast(model, h = 1)
  # save mean and sd of the gaussian predictive distribution
  fc[[s.name]] <- c(
    s.fc$mean,
    (s.fc$upper[, "95%"][[1]] - s.fc$mean) / qnorm(0.975)
  )
  residuals[, fc.idx] <- s.fc$residuals
  fc.idx <- fc.idx + 1
}
#> [1] "Forecasting at total..."
#> [1] "Forecasting at NSW..."
#> [1] "Forecasting at VIC..."
#> [1] "Forecasting at QLD..."
#> [1] "Forecasting at SA..."
#> [1] "Forecasting at WA..."
#> [1] "Forecasting at NT..."
#> [1] "Forecasting at ACT..."
#> [1] "Forecasting at TAS..."
#> [1] "Forecasting at male..."
#> [1] "Forecasting at female..."
#> [1] "Forecasting at NSW male..."
#> [1] "Forecasting at NSW female..."
#> [1] "Forecasting at VIC male..."
#> [1] "Forecasting at VIC female..."
#> [1] "Forecasting at QLD male..."
#> [1] "Forecasting at QLD female..."
#> [1] "Forecasting at SA male..."
#> [1] "Forecasting at SA female..."
#> [1] "Forecasting at WA male..."
#> [1] "Forecasting at WA female..."
#> [1] "Forecasting at NT male..."
#> [1] "Forecasting at NT female..."
#> [1] "Forecasting at ACT male..."
#> [1] "Forecasting at ACT female..."
#> [1] "Forecasting at TAS male..."
#> [1] "Forecasting at TAS female..."
```

Now we build the $\mathbf{A}$ matrix.

``` r
# we have 16 bottom time series, and 11 upper time series
A <- matrix(data = c(
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
  1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
  0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1
), byrow = TRUE, ncol = 16)

# plot of A
par(mai = c(1.5, 1, 0.5, 0.5))
image(1:ncol(A), 1:nrow(A),
  t(apply(t(A), 1, rev)),
  xaxt = "n", yaxt = "n", ann = FALSE
)
axis(1, at = 1:ncol(A), label = names(infantMortality)[12:27], las = 2)
axis(2, at = c(1:11), label = rev(names(infantMortality)[1:11]), las = 2)
```

![\*\*Figure 6\*\*: Infants mortality - The aggregation matrix A (red=1,
yellow=0).](bayesRecon_files/figure-html/infants-s-1.png)

**Figure 6**: Infants mortality - The aggregation matrix A (red=1,
yellow=0).

  
We use
[`bayesRecon::schaferStrimmer_cov`](https://idsia.github.io/bayesRecon/reference/schaferStrimmer_cov.md)
to estimate the covariance matrix of the residuals with shrinkage
(Schäfer and Strimmer 2005).

``` r
# means
mu <- sapply(fc, "[[", 1)
# Shrinkage covariance
shrink.res <- bayesRecon::schaferStrimmer_cov(residuals)
print(paste("The estimated shrinkage intensity is", round(shrink.res$lambda_star, 3)))
#> [1] "The estimated shrinkage intensity is 0.153"
Sigma <- shrink.res$shrink_cov
```

We now perform Gaussian reconciliation:

``` r
recon.gauss <- bayesRecon::reconc_gaussian(A,
  base_fc_mean = mu,
  base_fc_cov = Sigma
)

bottom_mu_reconc <- recon.gauss$bottom_rec_mean
bottom_Sigma_reconc <- recon.gauss$bottom_rec_covariance

# Obtain reconciled mu and Sigma for the upper variable
upper_mu_reconc <- A %*% bottom_mu_reconc
upper_Sigma_reconc <- A %*% bottom_Sigma_reconc %*% t(A)

upper_mu_reconc
#>            [,1]
#>  [1,] 688.82408
#>  [2,] 190.05371
#>  [3,] 165.26902
#>  [4,] 173.39776
#>  [5,]  40.95861
#>  [6,]  63.82008
#>  [7,]  20.37299
#>  [8,]  17.24244
#>  [9,]  17.70948
#> [10,] 424.30287
#> [11,] 264.52121
```

## References

Corani, Giorgio, Dario Azzimonti, João P. S. C. Augusto, and Marco
Zaffalon. 2021. “Probabilistic Reconciliation of Hierarchical Forecast
via Bayes’ Rule.” In *Machine Learning and Knowledge Discovery in
Databases*, edited by Frank Hutter, Kristian Kersting, Jefrey Lijffijt,
and Isabel Valera, 211–26. Springer International Publishing.
<https://doi.org/10.1007/978-3-030-67664-3_13>.

Corani, Giorgio, Dario Azzimonti, and Nicolò Rubattu. 2024.
“Probabilistic Reconciliation of Count Time Series.” *International
Journal of Forecasting* 40 (2): 457–69.
https://doi.org/<https://doi.org/10.1016/j.ijforecast.2023.04.003>.

Hyndman, Rob. 2018. *Mcomp: Data from the m-Competitions*.
<https://CRAN.R-project.org/package=Mcomp>.

Hyndman, Rob J. 2015. *Expsmooth: Data Sets from "Forecasting with
Exponential Smoothing"*. <https://CRAN.R-project.org/package=expsmooth>.

Hyndman, Rob J, and Yeasmin Khandakar. 2008. “Automatic Time Series
Forecasting: The Forecast Package for R.” *Journal of Statistical
Software* 26 (3): 1–22. <https://doi.org/10.18637/jss.v027.i03>.

Hyndman, Rob, Alan Lee, Earo Wang, and Shanika Wickramasuriya. 2021.
*Hts: Hierarchical and Grouped Time Series*.
<https://CRAN.R-project.org/package=hts>.

Makridakis, Spyros, and Michele Hibon. 2000. “The M3-Competition:
Results, Conclusions and Implications.” *International Journal of
Forecasting* 16 (4): 451–76.

Schäfer, Juliane, and Korbinian Strimmer. 2005. “A Shrinkage Approach to
Large-Scale Covariance Matrix Estimation and Implications for Functional
Genomics.” *Statistical Applications in Genetics and Molecular Biology*
4 (1).

Zambon, Lorenzo, Arianna Agosto, Paolo Giudici, and Giorgio Corani.
2024. “Properties of the Reconciled Distributions for Gaussian and Count
Forecasts.” *International Journal of Forecasting* 40 (4): 1438–48.
https://doi.org/<https://doi.org/10.1016/j.ijforecast.2023.12.004>.

Zambon, Lorenzo, Dario Azzimonti, and Giorgio Corani. 2024. “Efficient
Probabilistic Reconciliation of Forecasts for Real-Valued and Count Time
Series.” *Statistics and Computing* 34 (1): 21.
<https://doi.org/10.1007/s11222-023-10343-y>.
