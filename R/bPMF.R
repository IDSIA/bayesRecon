# A bivariate pmf is represented as a normalized numeric matrix M: 
# for each (i,j) in {0,...,Mi} x {0,...,Mj}, the probability of (i,j) is the value M[i+1,j+1]

.LAMBDA_SHR      = 1e-2
.NEG_TOLL_BPMF   = 1e-3
.TOLL_IPFP       = 1e-6
.N_MAX_IPFP      = 1000

################################################################################
# Copula CDFs

# Compute Clayton CDF at given points u and v, which can be vectors
.clayton_cdf = function(u, v, theta) {
  z = u^(-theta) + v^(-theta) - 1
  z[z<0] = 0
  return(z^(-1/theta))
}

# Function to compute the Clayton copula CDF, given the marginal CDFs F1 and F2.
# Vectorized implementation
.get_cdf_clayton = function(F1, F2, theta) {
  
  M = length(F1)
  N = length(F2)
  
  F_grid = expand.grid(F1,F2)
  cdf = .clayton_cdf(F_grid[,1], F_grid[,2], theta)
  
  cdf = matrix(cdf, nrow = M, ncol = N)
  
  return(cdf)
}

# Compute Gumbel CDF at given points u and v, which can be vectors
.gumbel_cdf = function(u, v, theta) {
  # Clip u and v between toll and 1
  toll = 1e-9
  u = pmin(pmax(toll, u), 1)
  v = pmin(pmax(toll, v), 1)
  
  z = exp(-((-log(u))^theta + (-log(v))^theta)^(1/theta))
  return(z)
}

# Function to compute the Gumbel copula CDF, given the marginal CDFs F1 and F2.
# Vectorized implementation
.get_cdf_gumbel = function(F1, F2, theta) {
  
  M = length(F1)
  N = length(F2)
  
  F_grid = expand.grid(F1,F2)
  cdf = .gumbel_cdf(F_grid[,1], F_grid[,2], theta)
  
  if (any(is.na(cdf))) browser()
  
  cdf = matrix(cdf, nrow = M, ncol = N)
  
  return(cdf)
}

################################################################################
# Utils

# Transform a CDF into PMF via finite differences
.from_cdf_to_bpmf = function(cdf, neg_toll = .NEG_TOLL_BPMF) {
  
  M = nrow(cdf)
  N = ncol(cdf)
  
  cdf_enlarged = matrix(0, nrow=M+1, ncol=N+1)
  cdf_enlarged[2:(M+1), 2:(N+1)] = cdf
  
  bpmf = diff(cdf_enlarged)   # finite differences on the rows
  bpmf = t(diff(t(bpmf)))     # finite differences on the columns
  
  if (any(bpmf < -neg_toll)) {
    warning("In the computation of the copula pmf same values were negative. 
             They have been set to zero.")
  } 
  bpmf[bpmf<0] = 0
  bpmf = bpmf / sum(bpmf)
  
  return(bpmf)
}

# Shrink bpmf towards the uniform
.shrink_bpmf = function(bpmf, lambda_shr = .LAMBDA_SHR) {
  bpmf = (1 - lambda_shr) * bpmf + 
               lambda_shr * matrix(1 / (nrow(bpmf)*ncol(bpmf)),
                                   nrow = nrow(bpmf), ncol = ncol(bpmf))
  return(bpmf)
}

# Compute L1 distance between the 1st marginal of bpmf and pmf
.dist_marg_rows = function(bpmf, pmf) {
  return(sum(abs(rowSums(bpmf) - pmf)))
}

# Compute L1 distance between the 2nd marginal of bpmf and pmf
.dist_marg_cols = function(bpmf, pmf) {
  return(sum(abs(colSums(bpmf) - pmf)))
}

# Iterative Proportional Fitting Procedure (IPFP)
.IPFP = function(bpmf, pmf1, pmf2, toll = .TOLL_IPFP, N_max = .N_MAX_IPFP) {
  
  d = .dist_marg_rows(bpmf, pmf1) + .dist_marg_cols(bpmf, pmf2)
  
  for (i in 1:N_max) {
    # Exit loop when desired tolerance is reached
    if (d < toll) break()
    # Alternate projection on rows and columns
    if (i %% 2 == 1) {
      bpmf = bpmf / rowSums(bpmf) * pmf1            # divide and multiply each row
      d = .dist_marg_rows(bpmf, pmf1) + .dist_marg_cols(bpmf, pmf2)
    } else {
      bpmf = bpmf %*% diag(pmf2 / colSums(bpmf))    # divide and multiply each column
      d = .dist_marg_rows(bpmf, pmf1) + .dist_marg_cols(bpmf, pmf2)
    }
  }
  
  return(bpmf)
}

################################################################################
# Copula PMFs

# Compute the bivariate Gaussian density at x, y (which can be vectors) 
.biv_gauss_pdf <- function(x, y, rho) {
  z     = (x^2 + y^2 - 2 * rho * x * y) / (2 * (1 - rho^2))
  coeff = 1 / (2 * pi * sqrt(1 - rho^2))
  return(coeff * exp(-z))
}

# Compute Gaussian CDF at given points u and v, which can be vectors
# C(u,v) = F(F1^(-1)(u), F2^(-1)(v))
.gaussC_pmf = function(F1, F2, ro, neg_toll = .NEG_TOLL_BPMF) {
  # Clip u and v between toll and 1-toll
  toll = 1e-9
  F1 = pmax(toll, pmin(F1, 1-toll))
  F2 = pmax(toll, pmin(F2, 1-toll))
  
  F_grid = expand.grid(F1,F2)
  
  x = qnorm(F_grid[,1])
  y = qnorm(F_grid[,2])
  
  numer = .biv_gauss_pdf(x, y, ro)
  denom = dnorm(x) * dnorm(y)
  
  bpmf = matrix(numer/denom, nrow = length(F1), ncol = length(F2))
  
  if (any(bpmf < -neg_toll)) {
    warning("In the computation of the copula pmf same values were negative. 
             They have been set to zero.")
  } 
  bpmf[bpmf<0] = 0
  bpmf = bpmf / sum(bpmf)
  
  return(bpmf)
}

################################################################################

# Given a pair of vectors of samples v1 and v2,
# and given a pair of marginal CDFs F1 and F2,
# fit the specified copula family and return the bPMF
.fit_param_bpmf_from_margins = function(v1, v2,                              
                                  F1, F2,                        
                                  family = "clayton",
                                  lambda_shr = .LAMBDA_SHR,
                                  neg_toll = .NEG_TOLL,
                                  check_in = TRUE) {
  
  # Check input if required
  if (check_in) {
    # Check that samples are discrete and non-negative
    .check_discrete_samples(v1)
    .check_discrete_samples(v2)
    if (any(v1<0) | any(v2<0)) {
      stop("Samples must be non-negative")
    }
    if (length(v1) != length(v2)) {
      stop("Input vectors must have the same length")
    }
  }
  
  ### Clayton ##################################################################
  if (family == "clayton") {                                           
    
    tau = cor(v1, v2, method = "kendall")
    tau = min(tau, 1 - 1e-3)  # set max to tau to avoid numerical problems
    theta = 2 * tau / (1-tau)
    
    if (theta == 0) {
      bpmf = matrix(1/(length(F1)*length(F2)), nrow = length(F1), ncol = length(F2)) 
      return(bpmf)
    } else {
      cdf  = .get_cdf_clayton(F1, F2, theta)
      bpmf = .from_cdf_to_bpmf(cdf)
    }
    
    ### Gumbel #################################################################
  } else if (family == "gumbel") {                                      
    
    tau = cor(v1, v2, method = "kendall")
    tau = min(tau, 1 - 1e-3)  # set max to tau to avoid numerical problems
    theta = 1 / (1-tau)
    
    cdf  = .get_cdf_gumbel(F1, F2, theta)
    bpmf = .from_cdf_to_bpmf(cdf)
    
    ### Gaussian ###############################################################
  } else if (family == "gauss") {                                     
    
    tau = cor(v1, v2, method = "kendall")
    ro = sin(pi * tau / 2)
    
    # Here we directly compute the bpmf, not the cdf
    bpmf = .gaussC_pmf(F1, F2, ro)
    
    ### t Student ##############################################################
  } else if (family == "t") {                                        
    
    # Estimate ro via kendall's tau
    tau = cor(v1, v2, method = "kendall")
    ro = sin(pi * tau / 2)
    
    # Estimate coefficient of tail dependence lambda_t
    # lambda_t = (...)
    
    # Grid search to find nu
    nus = 1:20
    lambdas = 2 * pt(-((nus+1)*(1-ro)/(1+ro))^0.5, nus+1)
    nu = nus[which.min(abs(lambdas - lambda_t))]
    
    # build copula pmf (TODO)
    
    ### Independent ############################################################
  } else if (family == "independent") {                            
    
    bpmf = matrix(1/(length(F1)*length(F2)), nrow = length(F1), ncol = length(F2))
    
    ### Automatic selection ####################################################
  } else if (family == "auto") {            
    
    # TODO: fit different copulas an compute AIC to choose
    stop("Automatic selection of copula not yet implemented")
    
  }
  
  ### Shrinkage towards the uniform ###
  bpmf = .shrink_bpmf(bpmf, .LAMBDA_SHR)
  
  # Iterative Proportional Fitting Procedure (IPFP)
  bpmf = .IPFP(bpmf, PMF.from_cdf(F1), PMF.from_cdf(F2))
  
  return(bpmf)
}

# Nonparametric estimation of the copula bPMF via kernel density estimation.
.fit_copula_kde = function(v1, v2,
                           max1 = NULL, max2 = NULL, 
                           estim_params = NULL,
                           check_in = TRUE) {
  
  # Check input if required
  if (check_in) {
    # Check that samples are discrete and non-negative
    .check_discrete_samples(v1)
    .check_discrete_samples(v2)
    if (any(v1<0) | any(v2<0)) {
      stop("Samples must be non-negative")
    }
    if (length(v1) != length(v2)) {
      stop("Input vectors must have the same length")
    }
  }
  
  # If the maximum of the support is not specified, it is given by the max of v1 and v2
  if (is.null(max1)) {max1 = max(v1)}
  if (is.null(max2)) {max2 = max(v2)}
  
  stop("Copula KDE not yet implemented")
  
}

.empirical_bPMF_from_samples = function(v1_, v2_,
                                        max_supp1 = NULL, max_supp2 = NULL,
                                        smoothing = FALSE,
                                        check_input = TRUE) {
  
  # Remove NA
  v1 = v1_[!is.na(v1_)]
  v2 = v2_[!is.na(v2_)]
  
  if (check_input) {
    # Check that samples are discrete and non-negative
    .check_discrete_samples(v1)
    .check_discrete_samples(v2)
    if (any(v1<0) | any(v2<0)) {
      stop("Samples must be non-negative")
    }
  }
  
  # If the maximum of the support is not specified, it is given by the max of v1 and v2
  if (is.null(max_supp1)) {max_supp1 = max(v1)}
  if (is.null(max_supp2)) {max_supp2 = max(v2)}
  
  tab = table(factor(v1, levels = 0:max_supp1),
              factor(v2, levels = 0:max_supp2))
  bpmf = as.matrix(tab)
  rownames(bpmf) = NULL
  colnames(bpmf) = NULL
  
  bpmf = bpmf / sum(bpmf)       # normalize
  
  if (smoothing) {
    # TODO: implement smoothing
    warning("Smoothing not yet implemented")
  }
   
  return(bpmf)                           
  
}

# Estimate the bivariate pmf from a pair of vectors of samples
# v1 and v2 have to be of the same length;
# they can, however, have NA, corresponding to values that are not observed.
#
# Each marginal is estimated with all the available marginal data.
# There are several possibilities for the estimate; see PMF.from_samples
#
# The dependence structure (i.e., copula) is estimated using only pairs where
# both values are not NA.
# There are several possibilities for the estimate (...)
#
# Returns a bpmf (where v1 is on the rows, v2 on the columns)
bPMF.from_samples = function(v1_, v2_, 
                             estim_type_marginals = "parametric",
                             estim_params_marginals = NULL,
                             estim_type_joint = "parametric",
                             copula_family = "clayton",
                             estim_params_joint = NULL,
                             min_supp1 = NULL, min_supp2 = NULL,
                             L_max_copula = NULL) {
  
  # Remove NA
  v1 = v1_[!is.na(v1_)]
  v2 = v2_[!is.na(v2_)]
  
  # Check that samples are discrete and non-negative
  .check_discrete_samples(v1)
  .check_discrete_samples(v2)
  if (any(v1<0) | any(v2<0)) {
    stop("Samples must be non-negative")
  }
  
  # 1) Estimate the marginals
  pmf1 = PMF.from_samples(v1, 
                          estim_type = estim_type_marginals,
                          estim_params = estim_params_marginals,
                          min_supp = min_supp1,
                          check_in = FALSE)
  pmf2 = PMF.from_samples(v2, 
                          estim_type = estim_type_marginals,
                          estim_params = estim_params_marginals,
                          min_supp = min_supp2,
                          check_in = FALSE)
  
  # 2) Estimate the full bPMF
  
  # reduced vectors: I only keep indices for which both vectors have values
  v1_red = v1_[!is.na(v1_) & !is.na(v2_)]
  v2_red = v2_[!is.na(v1_) & !is.na(v2_)]
  
  # If L_max is specified, only use (at most) the last L_max observations
  L = length(v1_red)
  if (!is.null(L_max_copula) && L_max_copula < L) {
    v1_red = v1_red[(L-L_max_copula+1):L]
    v2_red = v2_red[(L-L_max_copula+1):L]
  }
  
  # If one of the reduced vectors has all elements equal, 
  # just assume independence:
  if (all(v1_red == v1_red[1]) || all(v2_red == v2_red[1])) {
    bpmf = outer(pmf1, pmf2)
    return(bpmf)
  }
  
  # Otherwise, estimate the copula pmf from (v1_red, v2_red):
  if (estim_type_joint == "naive") {
      bpmf = .empirical_bPMF_from_samples(v1_red, v2_red, 
                                          length(pmf1)-1, length(pmf2)-1,
                                          check_input = F)
      
      stop("Not yet implemented")
      
  } else if (estim_type_joint == "parametric") {
    
    F1 = cumsum(pmf1)
    F2 = cumsum(pmf2)
    bpmf = .fit_param_bpmf_from_margins(v1_red, v2_red,
                                        F1, F2,
                                        family = copula_family, 
                                        check_in = FALSE)
    
  }  else if (estim_type_joint == "kde") {
      
    bpmf = .fit_copula_kde(v1_red, v2_red,
                           length(pmf1)-1, length(pmf2)-1, 
                           check_in = FALSE) 
    }

  if (sum(bpmf)!=1) bpmf = bpmf / sum(bpmf)
  
  return(bpmf)
}


#' @title Sample from the bivariate distribution given as a bPMF object. TODO!!
#'
#' @description
#' 
#' Samples (with replacement) from the bivariate probability distribution specified by `bpmf`.
#' 
#' @param bpmf the bPMF object.
#' @param N_samples number of samples. 
#' 
#' @return Matrix N_samples x 2 of bivariate samples drawn from the distribution 
#' specified by `bpmf`.  
#' 
#' @export
bPMF.sample = function(bpmf, N_samples) {
  Mi = nrow(bpmf)
  Mj = ncol(bpmf)
  s = sample(1:(Mi * Mj), 
             prob = as.vector(bpmf), 
             replace = TRUE, size = N_samples)
  si = (s-1) %% Mi
  sj = (s-1) %/% Mi
  S = matrix(c(si,sj), ncol=2)
  return(S)
}

