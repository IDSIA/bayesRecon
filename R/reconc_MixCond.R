###############################################################################
# Reconciliation with mixed-conditioning (Mix-Cond)
###############################################################################


#' @title Probabilistic forecast reconciliation of mixed hierarchies 
#'
#' @description
#'
#' `reconc_MixCond()` uses importance sampling to draw samples from the reconciled
#' forecast distribution, obtained via conditioning, in the case of a mixed hierarchy.
#'
#' `reconc_TDcond()` uses a top-down conditioning algorithm: first, upper base forecasts are
#' reconciled via conditioning using only the hierarchical constraints between the
#' upper; then, the bottom distributions are updated via a probabilistic top-down procedure.
#'
#' @details
#'
#' The base bottom forecasts `base_fc_bottom` must be a list of length n_bottom, where each element is either
#' * a PMF object (see details below), if `bottom_in_type='pmf'`;
#' * a vector of samples, if `bottom_in_type='samples'`;
#' * a list of parameters, if `bottom_in_type='params'`:
#'    * lambda for the Poisson base forecast if `distr`='poisson', see \link[stats]{Poisson};
#'    * size and prob (or mu) for the negative binomial base forecast if `distr`='nbinom',
#'      see \link[stats]{NegBinomial}.
#'
#' The base upper forecasts `base_fc_upper` must be a list containing the parameters of
#' the multivariate Gaussian distribution of the upper forecasts.
#' The list must contain only the named elements `mean` (vector of length n_upper)
#' and `cov` (n_upper x n_upper matrix).
#'
#' The order of the upper and bottom base forecasts must match the order of (respectively) the rows and the columns of A.
#'
#' A PMF object is a numerical vector containing the probability mass function of a discrete distribution.
#' Each element corresponds to the probability of the integers from 0 to the last value of the support.
#' See also [PMF] for functions that handle PMF objects.
#'
#' @section Warnings and errors:
#'
#' In `reconc_MixCond`, warnings are triggered from the importance sampling step if:
#' * weights are all zeros, then the upper forecast is ignored during reconciliation;
#' * the effective sample size is < 200;
#' * the effective sample size is < 1% of the sample size.
#' 
#' These warnings are an indication that the base forecasts might have issues.
#' Please check the base forecasts in case of warnings.
#' 
#' In `reconc_TDcond`, if some of the reconciled upper samples lie outside the support of the bottom-up
#' distribution, those samples are discarded; the remaining ones are resampled with
#' replacement, so that the number of output samples is equal to `num_samples`.
#' In this case, a warning is issued if `suppress_warnings=FALSE` (default is `TRUE`).
#' If the fraction of discarded samples is above 50%, the function returns an error.
#'
#' @param A Aggregation matrix (n_upper x n_bottom).
#' @param base_fc_bottom A list containing the bottom base forecasts, see details.
#' @param base_fc_upper A list containing the upper base forecasts, see details.
#' @param bottom_in_type A string with three possible values:
#'
#' * 'pmf' if the bottom base forecasts are in the form of pmf, see details;
#' * 'samples' if the bottom base forecasts are in the form of samples;
#' * 'params'  if the bottom base forecasts are in the form of estimated parameters.
#'
#' @param distr A string describing the type of bottom base forecasts ('poisson' or 'nbinom').
#'
#' This is only used if `bottom_in_type='params'`.
#'
#' @param num_samples Number of samples drawn from the reconciled distribution.
#'        This is ignored if `bottom_in_type='samples'`; in this case, the number of
#'        reconciled samples is equal to the number of samples of the base forecasts.
#'
#' @param return_type The return type of the reconciled distributions.
#'        A string with three possible values:
#'
#' * 'pmf' returns a list containing the reconciled marginal pmf objects;
#' * 'samples' returns a list containing the reconciled multivariate samples;
#' * 'all' returns a list with both pmf objects and samples.
#'
#' @param return_upper Logical, whether to return the reconciled parameters for the upper variables (default is TRUE).
#'
#' @param suppress_warnings Logical. If \code{TRUE}, no warnings about samples are triggered;
#'        if \code{FALSE}, warnings are generated. Default is \code{FALSE} for `reconc_MixCond`
#'        and \code{TRUE} for `reconc_TDcond`. See the respective sections above.
#' 
#' @param seed Seed for reproducibility.
#'
#' @return A list containing the reconciled forecasts. The list has the following named elements:
#'
#' * `bottom_rec`: a list containing the pmf, the samples (matrix n_bottom x `num_samples`) or both,
#'    depending on the value of `return_type`;
#' * `upper_rec`: a list containing the pmf, the samples (matrix n_upper x `num_samples`) or both,
#'    depending on the value of `return_type` (only if `return_upper = TRUE`).
#'
#' @references
#' Zambon, L., Azzimonti, D., Rubattu, N., Corani, G. (2024).
#' *Probabilistic reconciliation of mixed-type hierarchical time series*.
#' Proceedings of the Fortieth Conference on Uncertainty in Artificial Intelligence,
#' PMLR 244:4078-4095. <https://proceedings.mlr.press/v244/zambon24a.html>.
#'
#' @seealso [reconc_BUIS()], [reconc_gaussian()], [PMF]
#'
#' @name Mixed_reconciliation
NULL

#' @rdname Mixed_reconciliation
#' @examples
#'
#' library(bayesRecon)
#'
#' # Consider a simple hierarchy with two bottom and one upper
#' A <- matrix(c(1, 1), nrow = 1)
#' # The bottom forecasts are Poisson with lambda=15
#' lambda <- 15
#' n_tot <- 60
#' base_fc_bottom <- list()
#' base_fc_bottom[[1]] <- apply(matrix(seq(0, n_tot)), MARGIN = 1,
#'                              FUN = \(x) dpois(x, lambda = lambda))
#' base_fc_bottom[[2]] <- apply(matrix(seq(0, n_tot)), MARGIN = 1,
#'                              FUN = \(x) dpois(x, lambda = lambda))
#'
#' # The upper forecast is a Normal with mean 40 and std 5
#' base_fc_upper <- list(mean = 40, cov = matrix(5^2))
#'
#' # Reconcile with reconc_MixCond
#' res.mixCond <- reconc_MixCond(A, base_fc_bottom, base_fc_upper)
#'
#' # Note that the bottom distributions are slightly shifted to the right
#' PMF_summary(res.mixCond$bottom_rec$pmf[[1]])
#' PMF_summary(base_fc_bottom[[1]])
#'
#' PMF_summary(res.mixCond$bottom_rec$pmf[[2]])
#' PMF_summary(base_fc_bottom[[2]])
#'
#' # The upper distribution is slightly shifted to the left
#' PMF_summary(res.mixCond$upper_rec$pmf[[1]])
#' PMF_get_var(res.mixCond$upper_rec$pmf[[1]])
#'
#' @export
reconc_MixCond <- function(A, base_fc_bottom, base_fc_upper,
                           bottom_in_type = "pmf", distr = NULL,
                           num_samples = 2e4, return_type = "pmf",
                           return_upper = TRUE, suppress_warnings = FALSE, 
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Check inputs
  .check_input_TD(
    A, base_fc_bottom, base_fc_upper,
    bottom_in_type, distr,
    return_type
  )

  # Prepare samples from the base bottom distribution
  if (bottom_in_type == "pmf") {
    B <- lapply(base_fc_bottom, PMF_sample, N_samples = num_samples)
    B <- do.call("cbind", B) # matrix of bottom samples (N_samples x n_bottom)
  } else if (bottom_in_type == "samples") {
    B <- do.call("cbind", base_fc_bottom)
    num_samples <- nrow(B)
  } else if (bottom_in_type == "params") {
    L_pmf <- lapply(base_fc_bottom, PMF_from_params, distr = distr)
    B <- lapply(L_pmf, PMF_sample, N_samples = num_samples)
    B <- do.call("cbind", B) # matrix of bottom samples (N_samples x n_bottom)
  }

  # Get mean and covariance matrix of the MVN upper base forecasts
  mean_upper <- base_fc_upper$mean
  cov_upper <- as.matrix(base_fc_upper$cov)

  out <- .core_reconc_MixCond(
    A, B, mean_upper, cov_upper, num_samples,
    return_type = return_type, 
    return_ESS = FALSE,
    return_upper = return_upper,
    suppress_warnings = suppress_warnings
  )

  return(out)
}


#' Core Reconciliation via Importance Sampling for Mixed Hierarchies
#'
#' Internal function that performs the core reconciliation logic using importance sampling
#' (IS) to reconcile mixed-type hierarchies. The base bottom forecasts (provided as samples)
#' are reweighted according to their fit to the upper multivariate Gaussian forecasts.
#'
#' @param A Matrix (n_upper x n_bottom) defining the hierarchy where upper = A %*% bottom.
#' @param B Matrix (n_samples x n_bottom) of bottom base forecast samples to be reconciled.
#' @param mean_upper Vector of upper level means.
#' @param cov_upper Covariance matrix of upper level.
#' @param num_samples Number of samples to draw/resample from.
#' @param return_type Character string specifying return format: 'pmf', 'samples', or 'all'.
#' @param return_ESS Logical, whether to return the Effective Sample Size (ESS) from importance sampling weights (default TRUE).
#' @param return_upper Logical, whether to return the reconciled parameters for the upper variables (default TRUE).
#' @param suppress_warnings Logical. If TRUE, suppresses warnings about sample quality (default FALSE).

#' @return A list containing:
#'   \itemize{
#'     \item `bottom_rec`: List with reconciled bottom forecasts (pmf and/or samples).
#'     \item `upper_rec`: (only if `return_upper = TRUE`) List with reconciled upper forecasts (pmf and/or samples).
#'     \item `ESS`: Effective Sample Size resulting from importance sampling reweighting (only if `return_ESS = TRUE`).
#'   }
#'
#' @keywords internal
#' @export
.core_reconc_MixCond <- function(A, B, mean_upper, cov_upper, num_samples, return_type, 
                                 return_ESS = TRUE, return_upper = TRUE,
                                 suppress_warnings = FALSE) {
  # Get dimensions
  n_u <- nrow(A)
  n_b <- ncol(A)

  # IS using MVN
  U <- B %*% t(A)
  weights <- .MVN_density(x = U, mu = mean_upper, Sigma = cov_upper)


  check_weights_res <- .check_weights(weights)
  if (check_weights_res$warning & !suppress_warnings) {
    warning_msg <- check_weights_res$warning_msg
    warning(warning_msg)
  }
  if (!(check_weights_res$warning & (1 %in% check_weights_res$warning_code))) {
    B <- .resample(B, weights, num_samples)
  }

  B <- t(B)
  U <- A %*% B

  # Prepare output: include the marginal pmfs and/or the samples (depending on "return" inputs)
  out <- list(bottom_rec = list(), upper_rec = list())
  if (return_type %in% c("pmf", "all")) {
    bottom_pmf <- lapply(1:n_b, function(i) PMF_from_samples(B[i, ]))
    out$bottom_rec$pmf <- bottom_pmf
    if (return_upper) {
      upper_pmf <- lapply(1:n_u, function(i) PMF_from_samples(U[i, ]))
      out$upper_rec$pmf <- upper_pmf
    }
  }
  if (return_type %in% c("samples", "all")) {
    out$bottom_rec$samples <- B
    if (return_upper) {
      out$upper_rec$samples <- U
    }
  }
  
  if (return_ESS) {
    out$ESS <- sum(weights)**2 / sum(weights**2)
  }

  return(out)
}
