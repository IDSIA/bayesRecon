# Base forecasts for the extreme market events dataset

Base forecasts for the `extr_mkt_events` dataset, computed using the
model by Agosto, A. (2022). *Multivariate Score-Driven Models for Count
Time Series to Assess Financial Contagion*.
[doi:10.2139/ssrn.4119895](https://doi.org/10.2139/ssrn.4119895) .

## Usage

``` r
extr_mkt_events_basefc
```

## Format

A list `extr_mkt_events_basefc` containing

- `extr_mkt_events_basefc$mu`:

  data frame of the base forecast means, for each day

- `extr_mkt_events_basefc$size`:

  data frame of the static base forecast size parameters

## Source

Agosto, A. (2022). *Multivariate Score-Driven Models for Count Time
Series to Assess Financial Contagion*.
[doi:10.2139/ssrn.4119895](https://doi.org/10.2139/ssrn.4119895)

## Details

The predictive distribution for the bottom time series is a multivariate
negative binomial with a static vector of dispersion parameters and a
time-varying vector of location parameters following a score-driven
dynamics. The base forecasts for the upper time series are computed
using a univariate version of this model. They are in-sample forecasts:
for each training instant, they are computed for time t+1 by
conditioning on the counts observed up to time t.

## References

Agosto, A. (2022). *Multivariate Score-Driven Models for Count Time
Series to Assess Financial Contagion*.
[doi:10.2139/ssrn.4119895](https://doi.org/10.2139/ssrn.4119895)

Zambon, L., Agosto, A., Giudici, P., Corani, G. (2024). *Properties of
the reconciled distributions for Gaussian and count forecasts*.
International Journal of Forecasting (in press).
[doi:10.1016/j.ijforecast.2023.12.004](https://doi.org/10.1016/j.ijforecast.2023.12.004)
.
