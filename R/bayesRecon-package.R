#' @section Learn more:
#'
#' To learn more about `bayesRecon`, start with the vignettes: `browseVignettes(package = "bayesRecon")`
#'
#' @section Reconciliation functions:
#'
#' The package implements reconciliation via conditioning for probabilistic forecasts
#' of hierarchical time series. The reconciliation functions are:
#'
#' * [reconc_gaussian()]: reconciliation via conditioning assuming multivariate Gaussian
#'    base forecasts; this is done analytically;
#' * [reconc_t()]: reconciliation via conditioning assuming multivariate Gaussian
#'    base forecasts with uncertain covariance matrix; 
#'    the reconciled forecasts are multivariate Student-t; this is done analytically;
#' * [reconc_BUIS()]: reconciliation via conditioning of any probabilistic forecast
#'    via bottom-up importance sampling; an alternative method for discrete forecasts 
#'    is implemented in [reconc_MCMC()], but we recommend using `reconc_BUIS`;
#' * [reconc_MixCond()] and [reconc_TDcond()]: reconciliation of mixed hierarchies, where
#'    the upper forecasts are multivariate Gaussian and the bottom forecasts are discrete distributions;
#'    `reconc_MixCond` implements conditioning via importance sampling, 
#'    while `reconc_TDcond` implements top-down conditioning.
#'
#' @section Utility functions:
#'
#' * [temporal_aggregation()]: temporal aggregation of a given time series object of class \link[stats]{ts};
#' * [get_reconc_matrices()]: aggregation and summing matrices for a temporal hierarchy
#'    of time series from user-selected list of aggregation levels;
#' * [schaferStrimmer_cov()]: computes the Schäfer-Strimmer shrinkage estimator for the covariance matrix;
#' * [multi_log_score_optimization()]: estimates the optimal degrees of freedom for [reconc_t()] 
#'    by maximizing the leave-one-out (LOO) multivariate log density;
#' * [PMF]: functions for handling PMF objects (sampling, computing statistics like mean, variance, quantiles, and summaries).
#'
#' @references
#' 
#' Carrara, C., Corani, G., Azzimonti, D., & Zambon, L. (2025).
#' *Modeling the uncertainty on the covariance matrix for probabilistic forecast reconciliation*.
#' arXiv preprint arXiv:2506.19554. <https://arxiv.org/abs/2506.19554>.
#' 
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
#' Proceedings of the Fortieth Conference on Uncertainty in Artificial Intelligence,
#' PMLR 244:4078-4095. <https://proceedings.mlr.press/v244/zambon24a.html>.
#'
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
