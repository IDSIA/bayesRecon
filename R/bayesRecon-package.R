#' @section Learn more: 
#' 
#' To learn more about `bayesRecon`, start with the vignettes: `browseVignettes(package = "bayesRecon")`
#' 
#' @section Main functions:
#'
#' The package implements reconciliation via conditioning for probabilistic forecasts 
#' of hierarchical time series. The main functions are:
#'
#' * [reconc_gaussian()]: reconciliation via conditioning of multivariate Gaussian 
#'    base forecasts; this is done analytically;
#' * [reconc_BUIS()]: reconciliation via conditioning of any probabilistic forecast 
#'    via importance sampling; this is the recommended option for non-Gaussian base forecasts;
#' * [reconc_MCMC()]: reconciliation via conditioning of discrete probabilistic 
#'    forecasts via Markov Chain Monte Carlo;
#' * [reconc_MixCond()]: reconciliation via conditioning of mixed hierarchies, where 
#'    the upper forecasts are multivariate Gaussian and the bottom forecasts are discrete distributions;
#' * [reconc_TDcond()]: reconciliation via top-down conditioning of mixed hierarchies, where 
#'    the upper forecasts are multivariate Gaussian and the bottom forecasts are discrete distributions.
#' 
#' @section Utility functions:
#'
#' * [temporal_aggregation()]: temporal aggregation of a given time series object of class \link[stats]{ts};
#' * [get_reconc_matrices()]: aggregation and summing matrices for a temporal hierarchy 
#'    of time series from user-selected list of aggregation levels;
#' * [schaferStrimmer_cov()]: computes the Schäfer-Strimmer shrinkage estimator for the covariance matrix;
#' * [PMF.get_mean()], [PMF.get_var()], [PMF.get_quantile()], [PMF.summary()], [PMF.sample()]: 
#'    functions for handling PMF objects.
#'
#' @references
#' Corani, G., Azzimonti, D., Augusto, J.P.S.C., Zaffalon, M. (2021). 
#' *Probabilistic Reconciliation of Hierarchical Forecast via Bayes' Rule*. 
#' ECML PKDD 2020. Lecture Notes in Computer Science, vol 12459.
#' \doi{10.1007/978-3-030-67664-3_13}.
#' 
#' Corani, G., Azzimonti, D., Rubattu, N. (2024). 
#' *Probabilistic reconciliation of count time series*. 
#' International Journal of Forecasting 40 (2), 457-469.
#' \doi{10.1016/j.ijforecast.2023.04.003}.
#' 
#' Zambon, L., Azzimonti, D. & Corani, G. (2024). 
#' *Efficient probabilistic reconciliation of forecasts for real-valued and count time series*. 
#' Statistics and Computing 34 (1), 21.
#' \doi{10.1007/s11222-023-10343-y}.
#'
#' Zambon, L., Agosto, A., Giudici, P., Corani, G. (2024). 
#' *Properties of the reconciled distributions for Gaussian and count forecasts*. 
#' International Journal of Forecasting (in press).
#' \doi{10.1016/j.ijforecast.2023.12.004}.
#'
#' Zambon, L., Azzimonti, D., Rubattu, N., Corani, G. (2024). 
#' *Probabilistic reconciliation of mixed-type hierarchical time series*. 
#' The 40th Conference on Uncertainty in Artificial Intelligence, accepted.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
