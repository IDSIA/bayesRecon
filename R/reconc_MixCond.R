###############################################################################
# Reconciliation with mixed-conditioning (Mix-Cond)
###############################################################################



#' @title Probabilistic Reconciliation of forecasts via mixed conditioning
#'
#' @description
#'
#' Uses the mixed conditioning algorithm to draw samples from the reconciled
#' forecast distribution. 
#'
#' @details
#' 
#' The base (unreconciled) forecasts are passed with the parameters
#' 
#' * `fc_bottom`: a list of length n_bottom where each element is either a pmf object (`bottom_in_type='pmf'`),
#'    a vector of samples (`bottom_in_type='samples'`) or the parameters (`bottom_in_type='params'`) of the 
#'    parametric distribution specified in `distr`.
#' * `fc_upper`: a list containing the parameters of the multivariate Gaussian distribution of the upper forecasts. 
#'    The list must contain only the named elements `mu` (vector of length n_upper) and `Sigma` (n_upper x n_upper matrix) 
#'
#' 
#' A PMF object is a numerical vector containing the probability mass function for a forecast. 
#' The first element of a PMF vector for the variable X always contains the probability P(X=0) and 
#' there is one element for each integer. 
#' The last element of the PMF corresponds to the probability of the last value in the support of X.
#' See also \link{PMF.get_mean}, \link{PMF.get_var}, \link{PMF.sample}, \link{PMF.get_quantile}, \link{PMF.summary} for functions that handle PMF objects. 
#' 
#' 
#' A warnings is triggered if the intersection of the support for the reconciled uppers 
#' and the support of the bottom-up distribution is too small. In this case only 
#' few samples from the reconciled upper are kept. The warning reports the percentage
#' of samples kept. 
#' 
#' Note that warnings are an indication that the base forecasts might have issues. Please check the base forecasts in case of warnings.
#'
#' @param S Summing matrix (n x n_bottom).   return_pmf = TRUE, return_samples = FALSE, suppress_warnings = FALSE, seed = NULL
#' @param fc_bottom A list containing the bottom base forecasts, see details.
#' @param fc_upper A list containing the bottom base forecasts, see details.
#' @param bottom_in_type A string with three possible values:
#'
#' * 'pmf' if the bottom base forecasts are in the form of pmf, see details;
#' * 'samples' if the bottom base forecasts are in the form of samples;
#' * 'params'  if the bottom base forecasts are in the form of estimated parameters.
#'
#' @param distr A string describing the type of bottom base forecasts ('gaussian', 'poisson' or 'nbinom').
#' 
#' This is only used if `bottom_in_type=='params'`.
#'
#' @param num_samples Number of samples drawn from the reconciled distribution.
#' @param num_resample Number of importance sampling resamples. 
#' @param our_sampler TO BE REMOVED AFTER THE TESTS.
#' @param return_type The return type of the reconciled distributions. A string with three possible values:
#' 
#' * 'pmf' returns a list containing reconciled pmf objects;
#' * 'samples' returns a list containing reconciled samples;
#' * 'all' returns a list with both pmf objects and samples.
#' 
#' @param ... additional parameters to be passed to the smoothing functions for PMF objects.
#' 
#' @param suppress_warnings Logical. If \code{TRUE}, no warnings about samples
#'        are triggered. If \code{FALSE}, warnings are generated. Default is \code{FALSE}. See Details.
#' @param seed Seed for reproducibility.
#'
#' @return A list containing the reconciled forecasts. The list has the following named elements:
#'
#' * `bottom_reconciled`: a list containing the pmf, the samples (matrix n_bottom x `num_samples`) or both, depending on the value of `return_type`;
#' * `upper_reconciled`: a list containing the pmf, the samples (matrix n_upper x `num_samples`) or both, depending on the value of `return_type`.
#'
#' @examples
#'
#' library(bayesRecon)
#' 
#' # Consider a simple hierarchy with two bottom and one upper
#' A <- matrix(c(1,1),nrow=1)
#' S <- rbind(A,diag(nrow=2))
#' # The bottom forecasts are Poisson with lambda=15
#' lambda <- 15
#' n_tot <- 60
#' fc_bottom <- list()
#' fc_bottom[[1]] <- apply(matrix(seq(0,n_tot)),MARGIN=1,FUN=function(x) dpois(x,lambda=lambda))
#' fc_bottom[[2]] <- apply(matrix(seq(0,n_tot)),MARGIN=1,FUN=function(x) dpois(x,lambda=lambda))
#' 
#' # The upper forecast is a Normal with mean 40 and std 5
#' fc_upper<- list(mu=40, Sigma=matrix(c(5^2)))
#' 
#' # We can reconcile with reconc_TDcond
#' res.mixCond <- reconc_MixCond(S, fc_bottom, fc_upper)
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
#' Corani, G., Azzimonti, D., Augusto, J.P.S.C., Zaffalon, M. (2021). *Probabilistic Reconciliation of Hierarchical Forecast via Bayes' Rule*. In: Hutter, F., Kersting, K., Lijffijt, J., Valera, I. (eds) Machine Learning and Knowledge Discovery in Databases. ECML PKDD 2020. Lecture Notes in Computer Science(), vol 12459. Springer, Cham. \doi{10.1007/978-3-030-67664-3_13}.
#'
#' Zambon, L., Agosto, A., Giudici, P., Corani, G. (2023). *Properties of the reconciled distributions for Gaussian and count forecasts*. \doi{10.48550/arXiv.2303.15135}.
#'
#' Zambon, L., Azzimonti, D., Rubattu, N., Corani, G. (2024). *Probabilistic reconciliation of mixed-type hierarchical time series* The 40th Conference on Uncertainty in Artificial Intelligence, accepted.
#'
#' @seealso [reconc_TDcond], [reconc_BUIS()]
#'
#' @export
reconc_MixCond = function(S, fc_bottom, fc_upper, 
                         bottom_in_type = "pmf", distr = NULL,
                         num_samples = 2e4, num_resample = 2e4,
                         return_type = "pmf", 
                         ...,
                         suppress_warnings = FALSE, seed = NULL,our_sampler=TRUE) {
  
  set.seed(seed)
  
  # Parameters for convolution
  # toll=1e-16
  # Rtoll=1e-7
  # smooth_bottom=TRUE
  # al_smooth=NULL
  # lap_smooth=FALSE 
  
  # After testing the convolution parameters:
  # remove dots, remove comment above, and set the "best parameters" as default in 
  # PMF.check_support and .TD_sampling
  
  # Check inputs
  .check_input_TD(S, fc_bottom, fc_upper, 
                  bottom_in_type, distr,
                  return_type)
  
  # Get aggr. matrix A 
  A = .get_A_from_S(S)$A
  n_u = nrow(A)
  n_b = ncol(A)
  
  # Prepare samples from the base bottom distribution
  if (bottom_in_type == "pmf") {
    B = lapply(fc_bottom, PMF.sample, N_samples=num_samples)
    B = do.call("cbind", B)    # matrix of bottom samples (N_samples x n_bottom)
  } else if (bottom_in_type == "samples") {
    B = do.call("cbind", fc_bottom)
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
  if(our_sampler){
    weights = .MVN_density(x=U, mu = mu_u, Sigma = Sigma_u)
  }else{
    weights = emdbook::dmvnorm(U, mu = mu_u, Sigma = Sigma_u)
  }
  
  
  check_weights.res = .check_weigths(weights)
  if (check_weights.res$warning & !suppress_warnings) {
    warning_msg = check_weights.res$warning_msg
    warning(warning_msg)
  }
  if(!(check_weights.res$warning & (1 %in% check_weights.res$warning_code))){
    B = .resample(B, weights, num_resample)
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
