# Sample from the distribution p(B_1, B_2 | B_1 + B_2 = u),
# where B_1 and B_2 are distributed as pmf1 and pmf2.
# u is a vector
.cond_biv_sampling = function(u, pmf1, pmf2) {
  
  # In this way then we iterate over the one with shorter support:
  sw = FALSE
  if (length(pmf1) > length(pmf2)) {
    pmf_ = pmf1
    pmf1 = pmf2
    pmf2 = pmf_
    sw = TRUE
  }
  
  b1 = rep(NA, length(u))  # initialize empty vector
  
  for (u_uniq in unique(u)) {  # loop over different values of u
    
    len_supp1 = length(pmf1)
    supp1 = 0:(len_supp1-1)
    p1 = pmf1
    
    supp2 = u_uniq - supp1
    supp2[supp2<0] = Inf  # trick to get NA when we access pmf2 outside the support
    p2 = pmf2[supp2+1]    # +1 because support starts from 0, but vector indexing from 1
    p2[is.na(p2)] = 0     # set NA to zero
    
    p = p1 * p2
    p = p / sum(p)
    
    u_posit = (u == u_uniq)
    b1[u_posit] = sample(supp1, size = sum(u_posit), replace = TRUE, prob = p)
  }
  
  if (sw) b1 = u - b1   # if we have switched, switch back
  
  return(list(b1, u-b1))
}

# Given a vector u of the upper values and a list of the bottom distr pmfs,
# returns samples (dim: n_bottom x length(u)) from the conditional distr 
# of the bottom given the upper values
.TD_sampling = function(u, bott_pmf, 
                        toll=.TOLL, Rtoll=.RTOLL, smoothing=TRUE, 
                        al_smooth=.ALPHA_SMOOTHING, lap_smooth=.LAP_SMOOTHING) {
  
  l_l_pmf = rev(PMF.bottom_up(bott_pmf, toll = toll, Rtoll = Rtoll, return_all = TRUE, 
                              smoothing=smoothing, al_smooth=al_smooth, lap_smooth=lap_smooth))
  
  b_old = matrix(u, nrow = 1)
  for (l_pmf in l_l_pmf[2:length(l_l_pmf)]) {
    L = length(l_pmf)
    b_new = matrix(ncol = length(u), nrow = L)
    for (j in 1:(L%/%2)) {
      b = .cond_biv_sampling(b_old[j,], l_pmf[[2*j-1]], l_pmf[[2*j]])
      b_new[2*j-1,] = b[[1]]
      b_new[2*j,]   = b[[2]]
    }
    if (L%%2 == 1) b_new[L,] = b_old[L%/%2 + 1,]
    b_old = b_new
  }
  
  return(b_new)
}



#' @title Probabilistic forecast reconciliation of mixed hierarchies via top-down conditioning
#'
#' @description
#'
#' Uses the top-down conditioning algorithm to draw samples from the reconciled
#' forecast distribution. Reconciliation is performed in two steps: 
#' first, the upper base forecasts are reconciled via conditioning, 
#' using only the hierarchical constraints between the upper variables; then,
#' the bottom distributions are updated via a probabilistic top-down procedure.
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
#' and `Sigma` (n_upper x n_upper matrix) 
#'  
#' A PMF object is a numerical vector containing the probability mass function of a discrete distribution.
#' Each element corresponds to the probability of the integers from 0 to the last value of the support.
#' See also \link{PMF.get_mean}, \link{PMF.get_var}, \link{PMF.sample}, \link{PMF.get_quantile}, 
#' \link{PMF.summary} for functions that handle PMF objects. 
#' 
#' If some of the reconciled upper samples lie outside the support of the bottom-up distribution, 
#' those samples are discarded and a warning is triggered.
#' The warning reports the percentage of samples kept. 
#'
#' @param S Summing matrix (n x n_bottom).   
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
#' This is only used if `bottom_in_type=='params'`.
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
#' res.TDcond <- reconc_TDcond(S, fc_bottom, fc_upper)
#' 
#' # Note that the bottom distributions are shifted to the right
#' PMF.summary(res.TDcond$bottom_reconciled$pmf[[1]])
#' PMF.summary(fc_bottom[[1]])
#' 
#' PMF.summary(res.TDcond$bottom_reconciled$pmf[[2]])
#' PMF.summary(fc_bottom[[2]])
#' 
#' # The upper distribution remains similar
#' PMF.summary(res.TDcond$upper_reconciled$pmf[[1]])
#' PMF.get_var(res.TDcond$upper_reconciled$pmf[[1]])
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
#' Zambon, L., Azzimonti, D., Rubattu, N., Corani, G. (2024). 
#' *Probabilistic reconciliation of mixed-type hierarchical time series*. 
#' The 40th Conference on Uncertainty in Artificial Intelligence, accepted.
#' 
#' @seealso [reconc_MixCond()], [reconc_BUIS()]
#'
#' @export
reconc_TDcond = function(S, fc_bottom, fc_upper, 
                         bottom_in_type = "pmf", distr = NULL,
                         num_samples = 2e4, return_type = "pmf", 
                         suppress_warnings = FALSE, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Check inputs
  .check_input_TD(S, fc_bottom, fc_upper, 
                  bottom_in_type, distr,
                  return_type)
  
  # Get aggr. matrix A and find the "lowest upper" 
  A = .get_A_from_S(S)$A
  n_u = nrow(A)
  n_b = ncol(A)
  lowest_rows = .lowest_lev(A)
  n_u_low = length(lowest_rows)  # number of lowest upper
  
  # Get mean and covariance matrix of the MVN upper base forecasts
  mu_u    = fc_upper$mu
  Sigma_u = as.matrix(fc_upper$Sigma)
  
  ### Get upper samples
  if (n_u == n_u_low) {     
    # If all the upper are lowest-upper, just sample from the base distribution
    U = .MVN_sample(num_samples, mu_u, Sigma_u)   # (dim: num_samples x n_u_low)
    U = round(U)                 # round to integer
    U_js = asplit(U, MARGIN = 2) # split into list of column vectors
    
  } else {
    # Else, analytically reconcile the upper and then sample from the lowest-uppers
    
    # Get the aggregation matrix A_u and the summing matrix S_u for the upper sub-hierarchy
    A_u = .get_Au(A, lowest_rows)
    S_u = matrix(nrow = n_u, ncol = n_u_low)
    S_u[-lowest_rows,] = A_u
    S_u[lowest_rows,] = diag(n_u_low)
    
    # Analytically reconcile the upper
    rec_gauss_u = reconc_gaussian(S_u, mu_u, Sigma_u)
    
    # Sample from reconciled MVN on the lowest level of the upper (dim: num_samples x n_u_low)
    U = .MVN_sample(n_samples = num_samples,
                    mu    = rec_gauss_u$bottom_reconciled_mean, 
                    Sigma = rec_gauss_u$bottom_reconciled_covariance)  
    U = round(U)                 # round to integer
    U_js = asplit(U, MARGIN = 2) # split into list of column vectors
  }

  # Prepare list of bottom pmf
  if (bottom_in_type == "pmf") {
    L_pmf = fc_bottom
  } else if (bottom_in_type == "samples") {
    L_pmf = lapply(fc_bottom, PMF.from_samples)
  } else if (bottom_in_type == "params") {
    L_pmf = lapply(fc_bottom, PMF.from_params, distr = distr)
  }
  
  # Prepare list of lists of bottom pmf relative to each lowest upper
  L_pmf_js = list()   
  for (j in lowest_rows) {
    Aj = A[j,]
    L_pmf_js = c(L_pmf_js, list(L_pmf[as.logical(Aj)]))
  }
  
  # Check that each multiv. sample of U is contained in the supp of the bottom-up distr
  samp_ok = mapply(PMF.check_support, U_js, L_pmf_js)
  samp_ok = rowSums(samp_ok) == n_u_low
  # Only keep the "good" upper samples, and throw a warning if some samples are discarded:
  U_js = lapply(U_js, "[", samp_ok) 
  if (sum(samp_ok) != num_samples & !suppress_warnings) {
    # We round down to the nearest decimal
    warning(paste0("Only ", floor(sum(samp_ok)/num_samples*1000)/10, "% of the upper samples ",
                   "are in the support of the bottom-up distribution; ",
                   "the others are discarded."))
  }
  
  # Get bottom samples via the prob top-down
  B = list()
  for (j in 1:n_u_low) {
    B[[j]] = .TD_sampling(U_js[[j]], L_pmf_js[[j]])
  }
  B = do.call("rbind", B)  # dim: n_bottom x num_samples
  U = A %*% B              # dim: n_upper x num_samples
  
  # Prepare output: include the marginal pmfs and/or the samples (depending on "return" inputs)
  out = list(bottom_reconciled=list(), upper_reconciled=list())
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






















