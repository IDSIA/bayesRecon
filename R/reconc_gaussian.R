#'
#' @title Analytical reconciliation of Gaussian base forecasts
#'
#' @description
#' Closed form computation of the reconciled forecasts in case of Gaussian base forecasts.
#'
#' @param A aggregation matrix (n_upper x n_bottom).
#' @param base_forecasts.mu a vector containing the means of the base forecasts.
#' @param base_forecasts.Sigma a matrix containing the covariance matrix of the base forecasts.
#' @param residuals a matrix with the residuals of the base forecasts, with n_upper + n_bottom columns.
#' The covariance matrix of the base forecasts is computed from the residuals using the Schäfer Strimmer shrinkage estimator.
#' If base_forecasts.Sigma is provided, residuals are ignored.
#' @param return_uppers logical, whether to return the reconciled parameters for the upper variables (default is FALSE).
#'
#' @details
#' In the vector of the means of the base forecasts the order must be: first the upper,
#' then the bottom; the order within the uppers is given by the rows of A,
#' the order within the bottoms by the columns of A.
#' The order of the rows of the covariance matrix of the base forecasts is the same.
#'
#' Unless `return_uppers = TRUE`, the function returns only the reconciled parameters of the bottom variables.
#'
#' @return A list containing the reconciled forecasts. The list has the following named elements:
#'
#' * `bottom_reconciled_mean`: reconciled mean for the bottom forecasts;
#' * `bottom_reconciled_covariance`: reconciled covariance for the bottom forecasts;
#' * `upper_reconciled_mean`: reconciled mean for the upper forecasts (if `return_uppers = TRUE`);
#' * `upper_reconciled_covariance`: reconciled covariance for the upper forecasts (if `return_uppers = TRUE`).
#'
#'
#' @examples
#'
#' library(bayesRecon)
#'
#' #' # ---- Example 1: base forecasts are given ----
#' 
#' # Create a minimal hierarchy with 2 bottom and 1 upper variable
#' A <- get_reconc_matrices(agg_levels = c(1, 2), h = 2)$A
#'
#' # Set the parameters of the Gaussian base forecast distributions
#' mu1 <- 2
#' mu2 <- 4
#' muY <- 9
#' mus <- c(muY, mu1, mu2)  # vector of means
#'
#' sigma1 <- 2
#' sigma2 <- 2
#' sigmaY <- 3
#' sigmas <- c(sigmaY, sigma1, sigma2)
#' Sigma <- diag(sigmas^2)  # covariance matrix
#' 
#' analytic_rec <- reconc_gaussian(A,
#'   base_forecasts.mu = mus,
#'   base_forecasts.Sigma = Sigma
#' )
#'
#' bottom_mu_reconc <- analytic_rec$bottom_reconciled_mean
#' bottom_Sigma_reconc <- analytic_rec$bottom_reconciled_covariance
#'
#' # To obtain reconciled samples for the entire hierarchy, sample from the reconciled 
#' # bottom distribution and then aggregate using A. 
#' 
#' # Sample from the reconciled bottom-level Gaussian distribution
#' # First, compute the Cholesky decomposition of the reconciled covariance matrix:
#' chol_decomp <- chol(bottom_Sigma_reconc)
#' # Then, sample from the standard normal distribution and apply the transformation:
#' Z <- matrix(stats::rnorm(n = 2000), nrow = 2) 
#' B <- t(chol_decomp) %*% Z + matrix(rep(bottom_mu_reconc, 1000), nrow = 2) 
#'
#' # Aggregate bottom samples to get upper samples, then stack
#' U <- A %*% B
#' Y_reconc <- rbind(U, B)
#' 
#' cat("Dimensions of reconciled samples (upper + bottom):", dim(Y_reconc), "\n")
#'
#' 
#' # ---- Example 2: using residuals from fitted ETS models ----
#' \donttest{
#' if (requireNamespace("forecast", quietly = TRUE)) {
#'
#'   # Simulate 2 bottom series from AR(1) processes
#'   set.seed(1234)
#'   n_obs <- 200
#'   y1 <- arima.sim(model = list(ar = 0.8), n = n_obs)
#'   y2 <- arima.sim(model = list(ar = 0.5), n = n_obs)
#'
#'   # Upper series is the sum of the two bottom series
#'   y_upper <- y1 + y2
#'
#'   # Aggregation matrix A:
#'   A <- matrix(c(1, 1), nrow = 1)
#'
#'   # Fit additive ETS models
#'   fit1 <- forecast::ets(y1, additive.only = TRUE)
#'   fit2 <- forecast::ets(y2, additive.only = TRUE)
#'   fit_upper <- forecast::ets(y_upper, additive.only = TRUE)
#'
#'   # Point forecasts (h = 1):
#'   fc1 <- forecast::forecast(fit1, h = 1)$mean
#'   fc2 <- forecast::forecast(fit2, h = 1)$mean
#'   fc_upper <- forecast::forecast(fit_upper, h = 1)$mean
#'   point_fc <- c(fc_upper, fc1, fc2)
#'
#'   # Residuals matrix (T x n, columns in same order as point_fc)
#'   res <- cbind(residuals(fit_upper),
#'                residuals(fit1),
#'                residuals(fit2))
#'
#'   # Reconcile (covariance estimated internally via Schafer-Strimmer)
#'   result <- reconc_gaussian(A, base_forecasts.mu = point_fc, residuals = res, return_uppers = TRUE)
#'
#'   bottom_mu <- result$bottom_reconciled_mean
#'   bottom_Sigma <- result$bottom_reconciled_covariance
#'   upper_mu <- result$upper_reconciled_mean
#'   upper_Sigma <- result$upper_reconciled_covariance
#'
#'   # Print reconciled means
#'   cat("Reconciled bottom means:", round(bottom_mu, 3), "\n")
#'   cat("Reconciled upper mean:", round(A %*% bottom_mu, 3), "\n")
#' 
#'   # Print 95% predictions intervals
#'   cat("Reconciled bottom 95% prediction intervals:\n")
#'   for (i in 1:length(bottom_mu)) {
#'     lower <- bottom_mu[i] - 1.96 * sqrt(bottom_Sigma[i, i])
#'     upper <- bottom_mu[i] + 1.96 * sqrt(bottom_Sigma[i, i])
#'     cat(paste0("Bottom ", i, ": [", round(lower, 3), ", ", round(upper, 3), "]\n"))
#'   }
#'   cat("Reconciled upper 95% prediction interval:\n")
#'   lower <- upper_mu - 1.96 * sqrt(upper_Sigma[1, 1])
#'   upper <- upper_mu + 1.96 * sqrt(upper_Sigma[1, 1])
#'   cat(paste0("Upper: [", round(lower, 3), ", ", round(upper, 3), "]\n"))
#' 
#' }
#' }
#'
#' @references
#' Corani, G., Azzimonti, D., Augusto, J.P.S.C., Zaffalon, M. (2021).
#' *Probabilistic Reconciliation of Hierarchical Forecast via Bayes' Rule*.
#' ECML PKDD 2020. Lecture Notes in Computer Science, vol 12459.
#' \doi{10.1007/978-3-030-67664-3_13}.
#'
#' Zambon, L., Agosto, A., Giudici, P., Corani, G. (2024).
#' *Properties of the reconciled distributions for Gaussian and count forecasts*.
#' International Journal of Forecasting (in press).
#' \doi{10.1016/j.ijforecast.2023.12.004}.
#'
#' @seealso [reconc_BUIS()]
#'
#' @export
reconc_gaussian <- function(A, base_forecasts.mu,
                            base_forecasts.Sigma = NULL,
                            residuals = NULL,
                            return_uppers = FALSE) {
  # Check matrix A
  .check_A(A)
  k <- nrow(A) # number of upper TS
  m <- ncol(A) # number of bottom TS
  n <- length(base_forecasts.mu) # total number of TS
  if (!(k + m == n)) {
    stop("Input error: the shape of A is not correct")
  }

  # If residuals are not provided, base_forecasts.Sigma must be provided
  if (is.null(residuals)) {
    if (is.null(base_forecasts.Sigma)) {
      stop("Input error: either residuals or base_forecasts.Sigma must be provided")
    }
    if (!(nrow(base_forecasts.Sigma) == n)) {
      stop("Input error: nrow(base_forecasts.Sigma) != length(base_forecasts.mu)")
    }
    .check_cov(base_forecasts.Sigma, "Sigma", pd_check = FALSE, symm_check = TRUE)
  } else {
    if (!is.null(base_forecasts.Sigma)) {
      warning("Input warning: both residuals and base_forecasts.Sigma are provided, ignoring residuals")
      .check_cov(base_forecasts.Sigma, "Sigma", pd_check = FALSE, symm_check = TRUE)
    } else if (ncol(residuals) != n) {
      stop("Input error: ncol(residuals) != length(base_forecasts.mu)")
    } else {
      # Compute the covariance matrix of the base forecasts from the residuals
      base_forecasts.Sigma <- schaferStrimmer_cov(residuals)$shrink_cov
    }
  }

  Sigma_u <- base_forecasts.Sigma[1:k, 1:k]
  Sigma_b <- base_forecasts.Sigma[(k + 1):n, (k + 1):n]
  Sigma_ub <- base_forecasts.Sigma[1:k, (k + 1):n, drop = FALSE]
  mu_u <- base_forecasts.mu[1:k]
  mu_b <- base_forecasts.mu[(k + 1):n]

  # Formulation from:
  # Zambon, L., et al. "Properties of the reconciled distributions for
  # Gaussian and count forecasts." (2023)
  Sigma_ub_At <- tcrossprod(Sigma_ub, A)
  Sigma_b_At <- tcrossprod(Sigma_b, A)
  Q <- Sigma_u - Sigma_ub_At - t(Sigma_ub_At) + A %*% Sigma_b_At
  # we only need to check if Q is p.d.
  .check_cov(Q, "Q", pd_check = TRUE, symm_check = FALSE)
  temp_diff <- t(Sigma_ub) - Sigma_b_At
  K <- t(solve(Q, t(temp_diff))) # equal to temp_diff %*% Q^-1

  mu_b_tilde <- mu_b + K %*% (A %*% mu_b - mu_u)
  Sigma_b_tilde <- Sigma_b - K %*% t(temp_diff)

  out <- list(
    bottom_reconciled_mean = as.vector(mu_b_tilde),
    bottom_reconciled_covariance = Sigma_b_tilde
  )

  if (return_uppers) {
    out$upper_reconciled_mean <- as.vector(A %*% mu_b_tilde)
    out$upper_reconciled_covariance <- A %*% Sigma_b_tilde %*% t(A)
  }

  return(out)
}
