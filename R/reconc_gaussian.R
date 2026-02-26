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
#' # Create a minimal hierarchy with 2 bottom and 1 upper variable
#' A <- get_reconc_matrices(agg_levels = c(1, 2), h = 2)$A
#'
#' # Set the parameters of the Gaussian base forecast distributions
#' mu1 <- 2
#' mu2 <- 4
#' muY <- 9
#' mus <- c(muY, mu1, mu2)
#'
#' sigma1 <- 2
#' sigma2 <- 2
#' sigmaY <- 3
#' sigmas <- c(sigmaY, sigma1, sigma2)
#'
#' Sigma <- diag(sigmas^2) # need to transform into covariance matrix
#' analytic_rec <- reconc_gaussian(A,
#'   base_forecasts.mu = mus,
#'   base_forecasts.Sigma = Sigma
#' )
#'
#' bottom_mu_reconc <- analytic_rec$bottom_reconciled_mean
#' bottom_Sigma_reconc <- analytic_rec$bottom_reconciled_covariance
#'
#' # Obtain reconciled mu and Sigma for the upper variable
#' upper_mu_reconc <- A %*% bottom_mu_reconc
#' upper_Sigma_reconc <- A %*% bottom_Sigma_reconc %*% t(A)
#'
#' # Obtain reconciled mu and Sigma for the entire hierarchy
#' S <- rbind(A, diag(2)) # first, get summing matrix S
#' Y_mu_reconc <- S %*% bottom_mu_reconc
#' Y_Sigma_reconc <- S %*% bottom_Sigma_reconc %*% t(S) # note that this is a singular matrix
#'
#' # Obtain reconciled samples for the entire hierarchy:
#' # i.e., sample from the reconciled bottoms and multiply by S
#' chol_decomp <- chol(bottom_Sigma_reconc) # Compute the Cholesky Decomposition
#' Z <- matrix(stats::rnorm(n = 2000), nrow = 2) # Sample from standard normal
#' B <- t(chol_decomp) %*% Z + matrix(rep(bottom_mu_reconc, 1000), nrow = 2) # Apply the transformation
#'
#' U <- S %*% B
#' Y_reconc <- rbind(U, B)
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
