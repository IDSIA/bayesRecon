#' Example of a time series from carparts 
#'
#' A monthly time series from the `carparts` dataset, 51 observations, Jan 1998 - Mar 2002.
#'
#' @format
#' Univariate time series of class \link[stats]{ts}.
#' 
#' @references
#' Hyndman, R.J., Koehler, A.B., Ord, J.K., and Snyder, R.D., (2008).
#' Forecasting with exponential smoothing: the state space approach.
#' Springer Science & Business Media.
#' 
#' Godahewa, R., Bergmeir, C., Webb, G., Hyndman, R., & Montero-Manso, P. (2020). 
#' Car Parts Dataset (without Missing Values) (Version 2) 
#' \doi{10.5281/zenodo.4656021}
#' 
#' @source 
#' Godahewa, R., Bergmeir, C., Webb, G., Hyndman, R.J., & Montero-Manso, P. (2020). 
#' Car Parts Dataset (without Missing Values) (Version 2) 
#' \doi{10.5281/zenodo.4656021}
"carparts_example"


#' Infant Mortality grouped time series dataset
#'
#' A yearly grouped time series dataset, from 1901 to 2003, of infant mortality counts (deaths) in Australia; 
#' disaggregated by state (see below), and sex (male and female).
#' 
#' States: New South Wales (NSW), Victoria (VIC), Queensland (QLD), South Australia (SA), Western Australia 
#' (WA), Northern Territory (NT), Australian Capital Territory (ACT), and Tasmania (TAS).
#'
#' @format
#' List of time series of class \link[stats]{ts}.
#' 
#' @source 
#' hts package [CRAN](https://cran.r-project.org/package=hts)
#' 
#' @references 
#' Hyndman, R.J., Ahmed, R.A., Athanasopoulos, G., Shang, H.L. (2011).
#' Optimal combination forecasts for hierarchical time series. 
#' Computational Statistics and Data Analysis, 55(9), 2579-2589.
"infantMortality"


#' Example of a time series from the M3 forecasting competition
#'
#' A monthly time series, from the M3 forecasting competition ("N1485").
#'
#' @format
#' List of time series of class \link[stats]{ts}.
#' 
#' 
#' @source [https://forecasters.org/resources/time-series-data/m3-competition/](https://forecasters.org/resources/time-series-data/m3-competition/)
"M3_example"


#' Extreme market events dataset
#'
#' Count time series of extreme market events in five economic sectors.
#' The data refer to the trading days between 2004/12/31 and 2018/12/19 (3508 trading days in total).
#' 
#' The counts are computed by considering 29 companies included in the Euro Stoxx 
#' 50 index and observing if the value of the CDS spread on a given day exceeds 
#' the 90-th percentile of its distribution in the last trading year.
#' The companies are divided in the following  sectors: Financial (FIN), Information 
#' and Communication Technology (ICT), Manufacturing (MFG), Energy (ENG), and Trade (TRD). 
#' 
#' There are 6 time series: 
#' - 5 bottom time series, corresponding to the daily counts for each sector
#' - 1 upper time series, which is the sum of all the bottom (ALL)
#'
#' @format
#' A multivariate time series of class \link[stats]{ts}.
#' 
#' @source 
#' Zambon, L., Agosto, A., Giudici, P., Corani, G. (2024). 
#' *Properties of the reconciled distributions for Gaussian and count forecasts*. 
#' International Journal of Forecasting (in press).
#' \doi{10.1016/j.ijforecast.2023.12.004}.
#' 
#' @references 
#' Zambon, L., Agosto, A., Giudici, P., Corani, G. (2024). 
#' *Properties of the reconciled distributions for Gaussian and count forecasts*. 
#' International Journal of Forecasting (in press).
#' \doi{10.1016/j.ijforecast.2023.12.004}.
#'
#' Agosto, A. (2022). 
#' *Multivariate Score-Driven Models for Count Time Series to Assess Financial Contagion*. 
#' \doi{10.2139/ssrn.4119895}
"extr_mkt_events"


#' Base forecasts for the extreme market events dataset 
#'
#' Base forecasts for the `extr_mkt_events` dataset, computed using the model by 
#' Agosto, A. (2022). 
#' *Multivariate Score-Driven Models for Count Time Series to Assess Financial Contagion*. 
#' \doi{10.2139/ssrn.4119895}. 
#' 
#' The predictive distribution for the bottom time series is a multivariate negative 
#' binomial with a static vector of dispersion parameters and a time-varying vector 
#' of location parameters following a score-driven dynamics. 
#' The base forecasts for the upper time series are computed using a univariate version of this model.
#' They are in-sample forecasts: for each training instant, they are computed for 
#' time t+1 by conditioning on the counts observed up to time t. 
#'
#' @format
#' A list `extr_mkt_events_basefc` containing
#' \describe{
#'    \item{`extr_mkt_events_basefc$mu`}{data frame of the base forecast means, for each day}
#'    \item{`extr_mkt_events_basefc$size`}{data frame of the static base forecast size parameters}
#' }
#' 
#' @source 
#' Agosto, A. (2022). 
#' *Multivariate Score-Driven Models for Count Time Series to Assess Financial Contagion*. 
#' \doi{10.2139/ssrn.4119895}
#' 
#' @references 
#' Agosto, A. (2022). 
#' *Multivariate Score-Driven Models for Count Time Series to Assess Financial Contagion*. 
#' \doi{10.2139/ssrn.4119895}
#' 
#' Zambon, L., Agosto, A., Giudici, P., Corani, G. (2024). 
#' *Properties of the reconciled distributions for Gaussian and count forecasts*. 
#' International Journal of Forecasting (in press).
#' \doi{10.1016/j.ijforecast.2023.12.004}. 
"extr_mkt_events_basefc"


#' Example of hierarchical forecasts for a store from the M5 competition
#'
#' This dataset contains forecasts for the hierarchy of time series related to the store `CA_1` from the M5 competition. 
#' 
#' The store `CA_1` contains 3049 item level time series and 11 aggregate time series:
#' * Store level aggregation (`CA_1`)
#' * Category level aggregations  (`HOBBIES`, `HOUSEHOLD`, `FOODS`)
#' * Department level aggregations (`HOBBIES_1`, `HOBBIES_2`, `HOUSEHOLD_1`, `HOUSEHOLD_2`, `FOODS_1`, `FOODS_2`, `FOODS_3`)
#' 
#' Forecasts are generated with the function `forecast` and the model `adam` from the package `smooth`.
#' * The models for the bottom time series are selected with multiplicative Gamma error term (`MNN`);
#' * The models for the upper time series (`AXZ`) is selected with Gaussian additive error term, seasonality selected based on information criterion.
#' 
#' The raw data was downloaded with the package `m5`. 
#'
#' @format
#' A list containing:
#' * `upper`: a list of 11 elements each representing an aggregation level. Each element contains: `mu`, `sigma` the mean and standard deviation of the Gaussian forecast, `actual` the actual value, `residuals` the residuals of the model used to estimate forecasts covariance. 
#' * `lower`: a list of 3049 elements each representing a forecast for each item. Each element contains `pmf` the probability mass function of the item level forecast, `actual` the actual value. 
#' * `A`: the aggregation matrix for A.
#' * `S`: the S matrix for the hierarchy. 
#' * `Q_u`: scaling factors for computing MASE on the upper forecasts.
#' * `Q_b`: scaling factors for computing MASE on the bottom forecasts.
#' 
#' @references
#' Joachimiak K (2022). *m5: 'M5 Forecasting' Challenges Data*. R package version 0.1.1, <https://CRAN.R-project.org/package=m5>.
#' 
#' Makridakis, Spyros & Spiliotis, Evangelos & Assimakopoulos, Vassilis. (2020). *The M5 Accuracy competition: Results, findings and conclusions.* International Journal of Forecasting 38(4) 1346-1364. \doi{10.1016/j.ijforecast.2021.10.009}
#' 
#' Svetunkov I (2023). *smooth: Forecasting Using State Space Models*. R package version 4.0.0, <https://CRAN.R-project.org/package=smooth>.
#' 
#' @source 
#' Makridakis, Spyros & Spiliotis, Evangelos & Assimakopoulos, Vassilis. (2020). *The M5 Accuracy competition: Results, findings and conclusions.* International Journal of Forecasting 38(4) 1346-1364. \doi{10.1016/j.ijforecast.2021.10.009}
"M5_CA1_basefc"