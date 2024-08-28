###############################################################################
# Reconciliation with mixed-conditioning (Mix-Cond)
###############################################################################



#' @title Probabilistic forecast reconciliation of mixed hierarchies via conditioning
#'
#' @description
#'
#' Uses importance sampling to draw samples from the reconciled
#' forecast distribution, obtained via conditioning, in the case of a mixed hierarchy. 
#'
#' @details
#' 
#' The base bottom forecasts `fc_bottom` must be a list of length n_bottom, where each element is either
#' * a PMF object (see details below), if `bottom_in_type='pmf'`;
#' * a vector of samples, if `bottom_in_type='samples'`;
#' * a list of parameters, if `bottom_in_type='params'`:
#'    * lambda for the Poisson base forecast if `distr`='poisson', see \link[stats]{Poisson};
#'    * size and prob (or mu) for the negative binomial base forecast if `distr`='nbinom', 
#'      see \link[stats]{NegBinomial}.
#' 
#' The base upper forecasts `fc_upper` must be a list containing the parameters of 
#' the multivariate Gaussian distribution of the upper forecasts.
#' The list must contain only the named elements `mu` (vector of length n_upper) 
#' and `Sigma` (n_upper x n_upper matrix).
#' 
#' The order of the upper and bottom base forecasts must match the order of (respectively) the rows and the columns of A.
#'  
#' A PMF object is a numerical vector containing the probability mass function of a discrete distribution.
#' Each element corresponds to the probability of the integers from 0 to the last value of the support.
#' See also \link{PMF.get_mean}, \link{PMF.get_var}, \link{PMF.sample}, \link{PMF.get_quantile}, 
#' \link{PMF.summary} for functions that handle PMF objects. 
#' 
#' Warnings are triggered from the Importance Sampling step if:
#' 
#' * weights are all zeros, then the upper forecast is ignored during reconciliation;
#' * the effective sample size is < 200;
#' * the effective sample size is < 1% of the sample size.
#' 
#' Note that warnings are an indication that the base forecasts might have issues. 
#' Please check the base forecasts in case of warnings.
#'
#' @param A Aggregation matrix (n_upper x n_bottom).  
#' @param fc_bottom A list containing the bottom base forecasts, see details.
#' @param fc_upper A list containing the upper base forecasts, see details.
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
#' @param suppress_warnings Logical. If \code{TRUE}, no warnings about samples
#'        are triggered. If \code{FALSE}, warnings are generated. Default is \code{FALSE}. See Details.
#' @param seed Seed for reproducibility.
#'
#' @return A list containing the reconciled forecasts. The list has the following named elements:
#'
#' * `bottom_reconciled`: a list containing the pmf, the samples (matrix n_bottom x `num_samples`) or both, 
#'    depending on the value of `return_type`;
#' * `upper_reconciled`: a list containing the pmf, the samples (matrix n_upper x `num_samples`) or both, 
#'    depending on the value of `return_type`.
#'
#' @examples
#'
#' library(bayesRecon)
#' 
#' # Consider a simple hierarchy with two bottom and one upper
#' A <- matrix(c(1,1),nrow=1)
#' # The bottom forecasts are Poisson with lambda=15
#' lambda <- 15
#' n_tot <- 60
#' fc_bottom <- list()
#' fc_bottom[[1]] <- apply(matrix(seq(0,n_tot)),MARGIN=1,FUN=function(x) dpois(x,lambda=lambda))
#' fc_bottom[[2]] <- apply(matrix(seq(0,n_tot)),MARGIN=1,FUN=function(x) dpois(x,lambda=lambda))
#' 
#' # The upper forecast is a Normal with mean 40 and std 5
#' fc_upper<- list(mu=40, Sigma=matrix(5^2))
#' 
#' # We can reconcile with reconc_MixCond
#' res.mixCond <- reconc_MixCond(A, fc_bottom, fc_upper)
#' 
#' # Note that the bottom distributions are slightly shifted to the right
#' PMF.summary(res.mixCond$bottom_reconciled$pmf[[1]])
#' PMF.summary(fc_bottom[[1]])
#' 
#' PMF.summary(res.mixCond$bottom_reconciled$pmf[[2]])
#' PMF.summary(fc_bottom[[2]])
#' 
#' # The upper distribution is slightly shifted to the left
#' PMF.summary(res.mixCond$upper_reconciled$pmf[[1]])
#' PMF.get_var(res.mixCond$upper_reconciled$pmf[[1]])
#' 
#' @references
#' Zambon, L., Azzimonti, D., Rubattu, N., Corani, G. (2024). 
#' *Probabilistic reconciliation of mixed-type hierarchical time series*. 
#' The 40th Conference on Uncertainty in Artificial Intelligence, accepted.
#'
#' @seealso [reconc_TDcond()], [reconc_BUIS()]
#'
#' @export
reconc_MixCond = function(A, fc_bottom, fc_upper, 
                         bottom_in_type = "pmf", distr = NULL,
                         num_samples = 2e4, return_type = "pmf", 
                         suppress_warnings = FALSE, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Check inputs
  .check_input_TD(A, fc_bottom, fc_upper, 
                  bottom_in_type, distr,
                  return_type)
  
  n_u = nrow(A)
  n_b = ncol(A)
  
  # Prepare samples from the base bottom distribution
  if (bottom_in_type == "pmf") {
    B = lapply(fc_bottom, PMF.sample, N_samples=num_samples)
    B = do.call("cbind", B)    # matrix of bottom samples (N_samples x n_bottom)
  } else if (bottom_in_type == "samples") {
    B = do.call("cbind", fc_bottom)
    num_samples = nrow(B)
  } else if (bottom_in_type == "params") {
    L_pmf = lapply(fc_bottom, PMF.from_params, distr = distr)
    B = lapply(L_pmf, PMF.sample, N_samples=num_samples)
    B = do.call("cbind", B)    # matrix of bottom samples (N_samples x n_bottom)
  } 
  
  # Get mean and covariance matrix of the MVN upper base forecasts
  mu_u    = fc_upper$mu
  Sigma_u = as.matrix(fc_upper$Sigma)
  
  # IS using MVN
  U = B %*% t(A)
  weights = .MVN_density(x=U, mu = mu_u, Sigma = Sigma_u)
  
  
  check_weights.res = .check_weights(weights)
  if (check_weights.res$warning & !suppress_warnings) {
    warning_msg = check_weights.res$warning_msg
    warning(warning_msg)
  }
  if(!(check_weights.res$warning & (1 %in% check_weights.res$warning_code))){
    B = .resample(B, weights, num_samples)
  }
  
  ESS = sum(weights)**2/sum(weights**2)
  
  B = t(B)
  U = A %*% B
  

  # Prepare output: include the marginal pmfs and/or the samples (depending on "return" inputs)
  out = list(bottom_reconciled=list(), upper_reconciled=list(), ESS = ESS)
  if (return_type %in% c('pmf', 'all')) {
    upper_pmf  = lapply(1:n_u, function(i) PMF.from_samples(U[i,]))
    bottom_pmf = lapply(1:n_b, function(i) PMF.from_samples(B[i,]))
    
    out$bottom_reconciled$pmf = bottom_pmf
    out$upper_reconciled$pmf = upper_pmf
    
  }
  if (return_type %in% c('samples','all')) {
    
    out$bottom_reconciled$samples = B
    out$upper_reconciled$samples = U
    
  } 
  
  return(out)
}
