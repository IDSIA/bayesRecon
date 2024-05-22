#'
#' @title Analytical reconciliation of Gaussian base forecasts
#'
#' @description
#' Closed form computation of the reconciled forecasts in case of Gaussian base forecasts.
#'
#' @param S summing matrix (n x n_bottom).
#' @param base_forecasts.mu a vector containing the means of the base forecasts.
#' @param base_forecasts.Sigma a matrix containing the covariance matrix of the base forecasts.
#'
#' @details
#' The order of the base forecast means and covariance is given by the order of the time series in the summing matrix.
#' 
#' The function returns only the reconciled parameters of the bottom variables.
#' The reconciled upper parameters and the reconciled samples for the entire hierarchy can be obtained from the reconciled bottom parameters. 
#' See the example section.
#'
#'
#' @return A list containing the bottom reconciled forecasts. The list has the following named elements:
#'
#' * `bottom_reconciled_mean`: reconciled mean for the bottom forecasts;
#' * `bottom_reconciled_covariance`: reconciled covariance for the bottom forecasts.
#' 
#'
#' @examples
#'
#'library(bayesRecon)
#'
#'# Create a minimal hierarchy with 2 bottom and 1 upper variable
#'rec_mat <- get_reconc_matrices(agg_levels=c(1,2), h=2)
#'S <- rec_mat$S
#'A <- rec_mat$A
#'
#'#Set the parameters of the Gaussian base forecast distributions
#'mu1 <- 2
#'mu2 <- 4
#'muY <- 9
#'mus <- c(muY,mu1,mu2)
#'
#'sigma1 <- 2
#'sigma2 <- 2
#'sigmaY <- 3
#'sigmas <- c(sigmaY,sigma1,sigma2)
#'
#'Sigma <- diag(sigmas^2)  #need to transform into covariance matrix
#'analytic_rec <- reconc_gaussian(S, base_forecasts.mu = mus,
#'                                base_forecasts.Sigma = Sigma)
#'
#'bottom_mu_reconc <- analytic_rec$bottom_reconciled_mean
#'bottom_Sigma_reconc <- analytic_rec$bottom_reconciled_covariance
#'
#'# Obtain reconciled mu and Sigma for the upper variable
#'upper_mu_reconc <- A %*% bottom_mu_reconc
#'upper_Sigma_reconc <- A %*% bottom_Sigma_reconc %*% t(A)
#'
#'# Obtain reconciled mu and Sigma for the entire hierarchy
#'Y_mu_reconc <- S %*% bottom_mu_reconc
#'Y_Sigma_reconc <- S %*% bottom_Sigma_reconc %*% t(S)  # note: singular matrix
#'
#'# Obtain reconciled samples for the entire hierarchy:
#'# i.e., sample from the reconciled bottoms and multiply by S
#'chol_decomp = chol(bottom_Sigma_reconc) # Compute the Cholesky Decomposition
#'Z = matrix(stats::rnorm(n = 2000), nrow = 2) # Sample from standard normal
#'B = t(chol_decomp) %*% Z + matrix(rep(bottom_mu_reconc, 1000), nrow=2) # Apply the transformation
#'
#'U = S %*% B
#'Y_reconc = rbind(U, B)
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
reconc_gaussian <- function(S, base_forecasts.mu,
                            base_forecasts.Sigma) {
  # Check if S contains only 0s and 1s. 
  .check_S(S)
  hier = .get_A_from_S(S)
  A = hier$A
  k = nrow(A)    #number of upper TS
  m = ncol(A)    #number of bottom TS
  n = length(base_forecasts.mu) #total number of TS

  # Ensure that data inputs are valid
  if (!(nrow(base_forecasts.Sigma) == n)) {
    stop("Input error: nrow(base_forecasts.Sigma) != length(base_forecasts.mu)")
  }
  if (!(k + m == n)) {
    stop("Input error: the shape of S is not correct")
  }
  .check_cov(base_forecasts.Sigma, "Sigma", pd_check=FALSE, symm_check=TRUE)
  Sigma_u = base_forecasts.Sigma[hier$upper_idxs, hier$upper_idxs]
  Sigma_b = base_forecasts.Sigma[hier$bottom_idxs, hier$bottom_idxs]
  Sigma_ub = matrix(base_forecasts.Sigma[hier$upper_idxs, hier$bottom_idxs],
                    nrow = length(hier$upper_idxs))
  mu_u = base_forecasts.mu[hier$upper_idxs]
  mu_b = base_forecasts.mu[hier$bottom_idxs]

  # Formulation from:
  # Zambon, L., et al. "Properties of the reconciled distributions for
  # Gaussian and count forecasts." (2023)
  Q = Sigma_u - Sigma_ub %*% t(A) - A %*% t(Sigma_ub) + A %*% Sigma_b %*% t(A)
  # we only need to check if Q is p.d.
  .check_cov(Q, "Q", pd_check=TRUE, symm_check=FALSE)
  invQ = solve(Q)
  mu_b_tilde = mu_b + (t(Sigma_ub) - Sigma_b %*% t(A)) %*% invQ %*% (A %*% mu_b - mu_u)
  Sigma_b_tilde = Sigma_b - (t(Sigma_ub) - Sigma_b %*% t(A)) %*% invQ %*% t(t(Sigma_ub) - Sigma_b %*% t(A))

  out = list(
    bottom_reconciled_mean = mu_b_tilde,
    bottom_reconciled_covariance = Sigma_b_tilde
  )
  return(out)
}
