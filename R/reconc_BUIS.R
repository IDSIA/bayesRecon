###############################################################################
# Reconciliation with Bottom-Up Importance Sampling (BUIS)
###############################################################################

# Checks that there is no bottom continuous variable child of a 
# discrete upper variable
.check_hierfamily_rel <- function(sh.res, distr, debug=FALSE) {
  for (bi in seq_along(distr[sh.res$bottom_idxs])) {
    distr_bottom = distr[sh.res$bottom_idxs][[bi]]
    rel_upper_i = sh.res$A[,bi]
    rel_distr_upper = unlist(distr[sh.res$upper_idxs])[rel_upper_i == 1]
    err_message = "A continuous bottom distribution cannot be child of a discrete one."
    if (distr_bottom == "continuous") {              
      if (sum(rel_distr_upper == "discrete") | sum(rel_distr_upper %in% .DISCR_DISTR)) {
        if (debug) { return(-1) } else { stop(err_message) }
      }
    }
    if (distr_bottom %in% .CONT_DISTR) {                
      if (sum(rel_distr_upper == "discrete") | sum(rel_distr_upper %in% .DISCR_DISTR)) {
        if (debug) { return(-1) } else { stop(err_message) }
      }
    }
  }
  if (debug) { return(0) }
}

.emp_pmf <- function(l, density_samples) {
  empirical_pmf = PMF.from_samples(density_samples)
  w = sapply(l, function(i) empirical_pmf[i + 1])
  return(w)
}

.compute_weights <- function(b, u, in_type_, distr_) {
  if (in_type_ == "samples") {
    if (distr_ == "discrete") {
      # Discrete samples
      w = .emp_pmf(b, u)
    } else if (distr_ == "continuous") {
      # KDE
      d = stats::density(u, bw = "SJ", n = 2 ** 16)
      df = stats::approxfun(d)
      w = df(b)
    }
    # be sure no NA are returned, if NA, we want 0:
    # for the discrete branch:   if b_i !in u        --> NA
    # for the continuous branch: if b_i !in range(u) --> NA
    w[is.na(w)] = 0
  } else if (in_type_ == "params") {
    w = .distr_pmf(b, u, distr_)   # this never returns NA
  }
  # be sure not to return all 0 weights, return ones instead
  # if (sum(w) == 0) { w = w + 1 }
  return(w)
}

#' @title BUIS for Probabilistic Reconciliation of forecasts via conditioning
#'
#' @description
#'
#' Uses the Bottom-Up Importance Sampling algorithm to draw samples from the reconciled
#' forecast distribution, obtained via conditioning.
#'
#' @details
#'
#' The parameter `base_forecast` is a list containing n elements where the i-th element depends on
#' the values of `in_type[[i]]` and `distr[[i]]`.
#'
#' If `in_type[[i]]`='samples', then `base_forecast[[i]]` is a vector containing samples from the base forecast distribution.
#'
#' If `in_type[[i]]`='params', then `base_forecast[[i]]` is a list containing the estimated:
#'
#' * mean and sd for the Gaussian base forecast if `distr[[i]]`='gaussian', see \link[stats]{Normal};
#' * lambda for the Poisson base forecast if `distr[[i]]`='poisson', see \link[stats]{Poisson};
#' * size and prob (or mu) for the negative binomial base forecast if `distr[[i]]`='nbinom', see \link[stats]{NegBinomial}.
#' 
#' See the description of the parameters `in_type` and `distr` for more details. 
#'
#' The order of the `base_forecast` list is given by the order of the time series in the summing matrix.
#' 
#' Warnings are triggered from the Importance Sampling step if:
#' 
#' * weights are all zeros, then the upper is ignored during reconciliation;
#' * the effective sample size is < 200;
#' * the effective sample size is < 1% of the sample size (`num_samples` if `in_type` is 'params' or the size of the base forecast if if `in_type` is 'samples').
#' 
#' Note that warnings are an indication that the base forecasts might have issues. 
#' Please check the base forecasts in case of warnings.
#'
#' @param S Summing matrix (n x n_bottom).
#' @param base_forecasts A list containing the base_forecasts, see details.
#' @param in_type A string or a list of length n. If it is a list the i-th element is a string with two possible values:
#'
#' * 'samples' if the i-th base forecasts are in the form of samples;
#' * 'params'  if the i-th base forecasts are in the form of estimated parameters.
#' 
#' If it `in_type` is a string it is assumed that all base forecasts are of the same type. 
#'
#' @param distr A string or a list of length n describing the type of base forecasts. 
#' If it is a list the i-th element is a string with two possible values:
#'
#' * 'continuous' or 'discrete' if `in_type[[i]]`='samples';
#' * 'gaussian', 'poisson' or 'nbinom' if `in_type[[i]]`='params'.
#' 
#' If `distr` is a string it is assumed that all distributions are of the same type.
#'
#' @param num_samples Number of samples drawn from the reconciled distribution. 
#'        This is ignored if `bottom_in_type='samples'`; in this case, the number of reconciled samples is equal to 
#'        the number of samples of the base forecasts. 
#'        
#' @param suppress_warnings Logical. If \code{TRUE}, no warnings about effective sample size
#'        are triggered. If \code{FALSE}, warnings are generated. Default is \code{FALSE}. See Details.
#' @param seed Seed for reproducibility.
#'
#' @return A list containing the reconciled forecasts. The list has the following named elements:
#'
#' * `bottom_reconciled_samples`: a matrix (n_bottom x `num_samples`) containing the reconciled samples for the bottom time series;
#' * `upper_reconciled_samples`: a matrix (n_upper x `num_samples`) containing the reconciled samples for the upper time series;
#' * `reconciled_samples`: a matrix (n x `num_samples`) containing the reconciled samples for all time series.
#'
#' @examples
#'
#'library(bayesRecon)
#'
#'# Create a minimal hierarchy with 2 bottom and 1 upper variable
#'rec_mat <- get_reconc_matrices(agg_levels=c(1,2), h=2)
#'S <- rec_mat$S
#'
#'
#'#1) Gaussian base forecasts
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
#'base_forecasts = list()
#'for (i in 1:nrow(S)) {
#' base_forecasts[[i]] = list(mean = mus[[i]], sd = sigmas[[i]])
#'}
#'
#'
#'#Sample from the reconciled forecast distribution using the BUIS algorithm
#'buis <- reconc_BUIS(S, base_forecasts, in_type="params",
#'                  distr="gaussian", num_samples=100000, seed=42)
#'
#'samples_buis <- buis$reconciled_samples
#'
#'#In the Gaussian case, the reconciled distribution is still Gaussian and can be
#'#computed in closed form
#'Sigma <- diag(sigmas^2)  #transform into covariance matrix
#'analytic_rec <- reconc_gaussian(S, base_forecasts.mu = mus,
#'                                 base_forecasts.Sigma = Sigma)
#'
#'#Compare the reconciled means obtained analytically and via BUIS
#'print(c(S %*% analytic_rec$bottom_reconciled_mean))
#'print(rowMeans(samples_buis))
#'
#'
#'#2) Poisson base forecasts
#'
#'#Set the parameters of the Poisson base forecast distributions
#'lambda1 <- 2
#'lambda2 <- 4
#'lambdaY <- 9
#'lambdas <- c(lambdaY,lambda1,lambda2)
#'
#'base_forecasts <- list()
#'for (i in 1:nrow(S)) {
#'  base_forecasts[[i]] = list(lambda = lambdas[i])
#'}
#'
#'#Sample from the reconciled forecast distribution using the BUIS algorithm
#'buis <- reconc_BUIS(S, base_forecasts, in_type="params",
#'                           distr="poisson", num_samples=100000, seed=42)
#'samples_buis <- buis$reconciled_samples
#'
#'#Print the reconciled means
#'print(rowMeans(samples_buis))
#'
#' @references
#' Zambon, L., Azzimonti, D. & Corani, G. (2024). 
#' *Efficient probabilistic reconciliation of forecasts for real-valued and count time series*. 
#' Statistics and Computing 34 (1), 21.
#' \doi{10.1007/s11222-023-10343-y}.
#'
#'
#' @seealso
#' [reconc_gaussian()]
#'
#' @export
reconc_BUIS <- function(S,
                   base_forecasts,
                   in_type,
                   distr,
                   num_samples = 2e4,
                   suppress_warnings = FALSE,
                   seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)

  # Transform distr and in_type into lists
  if (!is.list(distr)) {
    distr = rep(list(distr), nrow(S))
  }
  if (!is.list(in_type)) {
    in_type = rep(list(in_type), nrow(S))
  }
  
  # Ensure that data inputs are valid
  .check_input_BUIS(S, base_forecasts, in_type, distr)

  # Split bottoms, uppers
  split_hierarchy.res = .split_hierarchy(S, base_forecasts)
  A = split_hierarchy.res$A
  upper_base_forecasts = split_hierarchy.res$upper
  bottom_base_forecasts = split_hierarchy.res$bottom
  
  # Check on continuous/discrete in relationship to the hierarchy
  .check_hierfamily_rel(split_hierarchy.res, distr)

  # H, G
  is.hier = .check_hierarchical(A)  
  # If A is hierarchical we do not solve the integer linear programming problem
  if(is.hier) {
    H = A
    G = NULL
    upper_base_forecasts_H = upper_base_forecasts
    upper_base_forecasts_G = NULL
    in_typeH = in_type[split_hierarchy.res$upper_idxs]
    distr_H  = distr[split_hierarchy.res$upper_idxs]
    in_typeG = NULL
    distr_G  = NULL
  } else {
    get_HG.res = .get_HG(A, upper_base_forecasts, distr[split_hierarchy.res$upper_idxs], in_type[split_hierarchy.res$upper_idxs])
    H = get_HG.res$H
    upper_base_forecasts_H = get_HG.res$Hv
    G = get_HG.res$G
    upper_base_forecasts_G = get_HG.res$Gv
    in_typeH = get_HG.res$Hin_type
    distr_H  = get_HG.res$Hdistr
    in_typeG = get_HG.res$Gin_type
    distr_G  = get_HG.res$Gdistr
  }
  
  # Reconciliation using BUIS
  n_upper = nrow(A)
  n_bottom = ncol(A)
  # 1. Bottom samples
  B = list()
  in_type_bottom = in_type[split_hierarchy.res$bottom_idxs]
  for (bi in 1:n_bottom) {
    if (in_type_bottom[[bi]] == "samples") {
      B[[bi]] = unlist(bottom_base_forecasts[[bi]])
    } else if (in_type_bottom[[bi]] == "params") {
      B[[bi]] = .distr_sample(bottom_base_forecasts[[bi]],
                              distr[split_hierarchy.res$bottom_idxs][[bi]],
                              num_samples)
    }
  }
  B = do.call("cbind", B) # B is a matrix (num_samples x n_bottom)

  # Bottom-Up IS on the hierarchical part
  for (hi in 1:nrow(H)) {
    c = H[hi, ]
    b_mask = (c != 0)
    weights = .compute_weights(
      b = (B %*% c),
      # (num_samples x 1)
      u = upper_base_forecasts_H[[hi]],
      in_type_ = in_typeH[[hi]],
      distr_ = distr_H[[hi]]
    )
    check_weights.res = .check_weigths(weights)
    if (check_weights.res$warning & !suppress_warnings) {
      warning_msg = check_weights.res$warning_msg
      # add information to the warning message
      upper_fromS_i = which(lapply(seq_len(nrow(S)), function(i) sum(abs(S[i,] - c))) == 0)
      for (wmsg in warning_msg) {
        wmsg = paste(wmsg, paste0("Check the upper forecast at index: ", upper_fromS_i,"."))
        warning(wmsg)
      }
    }
    if(check_weights.res$warning & (1 %in% check_weights.res$warning_code)){
      next
    }
    B[, b_mask] = .resample(B[, b_mask], weights)
  }

  if (!is.null(G)) {
    # Plain IS on the additional constraints
    weights = matrix(1, nrow = nrow(B))
    for (gi in 1:nrow(G)) {
      c = G[gi, ]
      weights = weights * .compute_weights(
        b = (B %*% c),
        u = upper_base_forecasts_G[[gi]],
        in_type_ = in_typeG[[gi]],
        distr_ = distr_G[[gi]]
      )
    }
    check_weights.res = .check_weigths(weights)
    if (check_weights.res$warning & !suppress_warnings) {
      warning_msg = check_weights.res$warning_msg
      # add information to the warning message
      upper_fromS_i = c()
      for (gi in 1:nrow(G)) {
        c = G[gi, ]
        upper_fromS_i = c(upper_fromS_i,
                          which(lapply(seq_len(nrow(S)), function(i) sum(abs(S[i,] - c))) == 0))
      }
      for (wmsg in warning_msg) {
        wmsg = paste(wmsg, paste0("Check the upper forecasts at index: ", paste0("{",paste(upper_fromS_i, collapse = ","), "}.")))
        warning(wmsg)
      }
    }
    if(!(check_weights.res$warning & (1 %in% check_weights.res$warning_code))){
      B = .resample(B, weights)
    }
    
  }

  B = t(B)
  U = A %*% B
  Y_reconc = rbind(U, B)

  out = list(
    bottom_reconciled_samples = B,
    upper_reconciled_samples = U,
    reconciled_samples = Y_reconc
  )
  return(out)
}

