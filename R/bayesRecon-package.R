#' @section Main functions:
#'
#' The package implements reconciliation via conditioning for probabilistic forecasts of hierarchical time series. The main functions are
#'
#' * [reconc_gaussian()]: analytical reconciliation of Gaussian base forecasts;
#' * [reconc_BUIS()]: reconciliation of any probabilistic base forecast via importance sampling;
#'    this is the recommended option for non-Gaussian base forecasts;
#' * [reconc_MCMC()]: reconciliation of probabilistic discrete base forecasts via Markov Chain Monte Carlo.
#'
#' @section Utility functions:
#'
#' * [temporal_aggregation()]: temporal aggregation of a given time series object of class \link[stats]{ts};
#' * [get_reconc_matrices()]: aggregation and summing matrices for a temporal hierarchy of time series from user-selected list of aggregation levels.
#'
#' @references
#' Corani, G., Azzimonti, D., Augusto, J.P.S.C., Zaffalon, M. (2021). *Probabilistic Reconciliation of Hierarchical Forecast via Bayes’ Rule*. In: Hutter, F., Kersting, K., Lijffijt, J., Valera, I. (eds) Machine Learning and Knowledge Discovery in Databases. ECML PKDD 2020. Lecture Notes in Computer Science(), vol 12459. Springer, Cham. \doi{10.1007/978-3-030-67664-3_13}.
#'
#' Corani, G., Rubattu, N., Azzimonti, D., Antonucci, A. (2022). *Probabilistic reconciliation of count time series*. \doi{10.48550/arXiv.2207.09322}.
#'
#' Zambon, L., Azzimonti, D. & Corani, G. (2022). *Efficient probabilistic reconciliation of forecasts for real-valued and count time series*. \doi{10.48550/arXiv.2210.02286}.
#'
#' Zambon, L., Agosto, A., Giudici, P., Corani, G. (2023). *Properties of the reconciled distributions for Gaussian and count forecasts*. \doi{10.48550/arXiv.2303.15135}.
#'
#'
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
