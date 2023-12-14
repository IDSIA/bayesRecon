###############################################################################
# Reconciliation with Bottom-Up Importance Sampling (BUIS)

.distr_sample <- function(params, distr_, n) {
  switch(
    distr_,
    "gaussian" = {
      samples = stats::rnorm(n=n, mean = params[[1]], sd = params[[2]]) },
    "poisson"  = {
      samples = stats::rpois(n=n, lambda = params[[1]]) },
    "nbinom"   = {
      samples <-if (params[[2]] == 0) {
          stop("Parameter size=0 in nbinom, this is not a valid distribution.")
          stats::rpois(n=n, lambda = params[[1]])
        } else {
          stats::rnbinom(n=n, mu = params[[1]], size = params[[2]])
        } },
  )
  return(samples)
}

.distr_pmf <- function(x, params, distr_) {
  switch(
    distr_,
    "gaussian" = {
      pmf = stats::dnorm(x=x, mean = params[[1]], sd = params[[2]]) },
    "poisson"  = {
      pmf = stats::dpois(x=x, lambda = params[[1]]) },
    "nbinom"   = {
      pmf = stats::dnbinom(x=x, mu = params[[1]], size = params[[2]]) },
  )
  return(pmf)
}

.emp_pmf <- function(l, density_samples) {
  empirical_pmf = sapply(0:max(density_samples), function(i)
    sum(density_samples == i) / length(density_samples))
  w = sapply(l, function(i) empirical_pmf[i + 1])
  return(w)
}

.check_weigths <- function(w, n_eff_min=200, p_n_eff=0.01) {
  warning = FALSE
  warning_code = c()
  warning_msg = c()
  
  n = length(w)
  n_eff = n
  
  
  # 1. w==0
  if (all(w==0)) {
    warning = TRUE
    warning_code = c(warning_code, 1) 
    warning_msg = c(warning_msg, 
                    "Importance Sampling: all the weights are zeros. This is probably caused by a strong incoherence between bottom and upper base forecasts.")
  }else{
    
    # Effective sample size
    n_eff = (sum(w)^2) / sum(w^2)
    
    # 2. n_eff < threshold
    if (n_eff < n_eff_min) {
      warning = TRUE
      warning_code = c(warning_code, 2) 
      warning_msg = c(warning_msg, 
                      paste0("Importance Sampling: effective_sample_size= ", round(n_eff,2), " (< ", n_eff_min,")."))
    }
    
    # 3. n_eff < p*n, e.g. p = 0.05
    if (n_eff < p_n_eff*n) {
      warning = TRUE
      warning_code = c(warning_code, 3) 
      warning_msg = c(warning_msg, 
                      paste0("Importance Sampling: effective_sample_size= ", round(n_eff,2), " (< ", round(p_n_eff * 100, 2),"%)."))
    }
  }
  res = list(warning = warning,
             warning_code = warning_code,
             warning_msg = warning_msg,
             n_eff = n_eff)
  
  return(res)
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

.resample <- function(S_, weights, num_samples = NA) {
  if (is.na(num_samples)) {
    num_samples = length(weights)
  }
  tmp_idx = sample(x = 1:num_samples, num_samples, replace = TRUE, prob = weights)
  return(S_[tmp_idx, ])
}

#' @title BUIS for Probabilistic Reconciliation of forecasts via conditioning
#'
#' @description
#'
#' Uses the Bottom-Up Importance Sampling algorithm to draw samples from the reconciled
#' forecast distribution, which is obtained via conditioning.
#'
#' @details
#'
#' The parameter `base_forecast` is a list containing n elements where the i-th element depends on
#' the values of `in_type[[i]]` and `distr[[i]]`.
#'
#' If `in_type[[i]]`='samples', then `base_forecast[[i]]` is a vector containing samples from the base forecast distribution.
#'
#' If `in_type[[i]]`='params', then `base_forecast[[i]]` is a vector containing the estimated:
#'
#' * mean and sd for the Gaussian base forecast if `distr[[i]]`='gaussian', see \link[stats]{Normal};
#' * lambda for the Poisson base forecast if `distr[[i]]`='poisson', see \link[stats]{Poisson};
#' * mu and size for the negative binomial base forecast if `distr[[i]]`='nbinom', see \link[stats]{NegBinomial}.
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
#' Note that warnings are an indication that the base forecasts might have issues. Please check the base forecasts in case of warnings.
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
#' @param distr A string or a list of length n describing the type of base forecasts. If it is a list the i-th element is a string with two possible values:
#'
#' * 'continuous' or 'discrete' if `in_type[[i]]`='samples';
#' * 'gaussian', 'poisson' or 'nbinom' if `in_type[[i]]`='params'.
#' 
#' If `distr` is a string it is assumed that all distributions are of the same type.
#'
#' @param num_samples Number of samples drawn from the reconciled distribution.
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
#' base_forecasts[[i]] = c(mus[[i]], sigmas[[i]])
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
#'  base_forecasts[[i]] = lambdas[i]
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
#' Zambon, L., Azzimonti, D. & Corani, G. (2024). *Efficient probabilistic reconciliation of forecasts for real-valued and count time series*. \doi{10.1007/s11222-023-10343-y}.
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
  set.seed(seed)

  # Ensure that data inputs are valid
  .check_input(S, base_forecasts, in_type, distr)
  if (!is.list(distr)) {
    distr = rep(list(distr), nrow(S))
  }
  
  if (!is.list(in_type)) {
    in_type = rep(list(in_type), nrow(S))
  }

  # Split bottoms, uppers
  split_hierarchy.res = .split_hierarchy(S, base_forecasts)
  A = split_hierarchy.res$A
  upper_base_forecasts = split_hierarchy.res$upper
  bottom_base_forecasts = split_hierarchy.res$bottom
  
  # Check on continuous/discrete in relationship to the hierarchy
  .check_hierfamily_rel(split_hierarchy.res, distr)

  # H, G
  is.hier = .check_hierarchical(A)
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
      u = unlist(upper_base_forecasts_H[[hi]]),
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
        u = unlist(upper_base_forecasts_G[[gi]]),
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
#'Z = matrix(rnorm(n = 2000), nrow = 2) # Sample from standard normal
#'B = chol_decomp %*% Z + matrix(rep(bottom_mu_reconc, 1000), nrow=2) # Apply the transformation
#'
#'U = S %*% B
#'Y_reconc = rbind(U, B)
#'
#' @references
#' Corani, G., Azzimonti, D., Augusto, J.P.S.C., Zaffalon, M. (2021). *Probabilistic Reconciliation of Hierarchical Forecast via Bayes' Rule*. In: Hutter, F., Kersting, K., Lijffijt, J., Valera, I. (eds) Machine Learning and Knowledge Discovery in Databases. ECML PKDD 2020. Lecture Notes in Computer Science(), vol 12459. Springer, Cham. \doi{10.1007/978-3-030-67664-3_13}.
#'
#' Zambon, L., Agosto, A., Giudici, P., Corani, G. (2023). *Properties of the reconciled distributions for Gaussian and count forecasts*. \doi{10.48550/arXiv.2303.15135}.
#'
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
  .check_cov(base_forecasts.Sigma)
  if (!(nrow(base_forecasts.Sigma) == n)) {
    stop("Input error: nrow(base_forecasts.Sigma) != length(base_forecasts.mu)")
  }
  if (!(k + m == n)) {
    stop("Input error: the shape of S is not correct")
  }

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
  invQ = solve(Q)
  mu_b_tilde = mu_b + (t(Sigma_ub) - Sigma_b %*% t(A)) %*% invQ %*% (A %*% mu_b - mu_u)
  Sigma_b_tilde = Sigma_b - (t(Sigma_ub) - Sigma_b %*% t(A)) %*% invQ %*% t(t(Sigma_ub) - Sigma_b %*% t(A))

  out = list(
    bottom_reconciled_mean = mu_b_tilde,
    bottom_reconciled_covariance = Sigma_b_tilde
  )
  return(out)
}
