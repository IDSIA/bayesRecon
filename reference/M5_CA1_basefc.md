# Example of hierarchical forecasts for a store from the M5 competition

This dataset contains forecasts for the hierarchy of time series related
to the store `CA_1` from the M5 competition.

## Usage

``` r
M5_CA1_basefc
```

## Format

A list containing:

- `upper`: a list of 11 elements each representing an aggregation level.
  Each element contains: `mu`, `sigma` the mean and standard deviation
  of the Gaussian forecast, `actual` the actual value, `residuals` the
  residuals of the model used to estimate forecasts covariance.

- `lower`: a list of 3049 elements each representing a forecast for each
  item. Each element contains `pmf` the probability mass function of the
  item level forecast, `actual` the actual value.

- `A`: the aggregation matrix for A.

- `S`: the S matrix for the hierarchy.

- `Q_u`: scaling factors for computing MASE on the upper forecasts.

- `Q_b`: scaling factors for computing MASE on the bottom forecasts.

## Source

Makridakis, Spyros & Spiliotis, Evangelos & Assimakopoulos, Vassilis.
(2020). *The M5 Accuracy competition: Results, findings and
conclusions.* International Journal of Forecasting 38(4) 1346-1364.
[doi:10.1016/j.ijforecast.2021.10.009](https://doi.org/10.1016/j.ijforecast.2021.10.009)

## Details

The store `CA_1` contains 3049 item level time series and 11 aggregate
time series:

- Store level aggregation (`CA_1`)

- Category level aggregations (`HOBBIES`, `HOUSEHOLD`, `FOODS`)

- Department level aggregations (`HOBBIES_1`, `HOBBIES_2`,
  `HOUSEHOLD_1`, `HOUSEHOLD_2`, `FOODS_1`, `FOODS_2`, `FOODS_3`)

Forecasts are generated with the function `forecast` and the model
`adam` from the package `smooth`.

- The models for the bottom time series are selected with multiplicative
  Gamma error term (`MNN`);

- The models for the upper time series (`AXZ`) is selected with Gaussian
  additive error term, seasonality selected based on information
  criterion.

The raw data was downloaded with the package `m5`.

## References

Joachimiak K (2022). *m5: 'M5 Forecasting' Challenges Data*. R package
version 0.1.1, <https://CRAN.R-project.org/package=m5>.

Makridakis, Spyros & Spiliotis, Evangelos & Assimakopoulos, Vassilis.
(2020). *The M5 Accuracy competition: Results, findings and
conclusions.* International Journal of Forecasting 38(4) 1346-1364.
[doi:10.1016/j.ijforecast.2021.10.009](https://doi.org/10.1016/j.ijforecast.2021.10.009)

Svetunkov I (2023). *smooth: Forecasting Using State Space Models*. R
package version 4.0.0, <https://CRAN.R-project.org/package=smooth>.
