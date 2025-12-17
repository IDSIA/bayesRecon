# Function for estimating the covariance matrix of the residuals of the naive or seasonal naive forecasts
# For each series, choose by using some criterion (RSS or test of seasonality)
# TODO: implement statistical test; leave dependence as suggested; default: choose using RSS
compute_naive_cov = function(y_train, freq = 1, criterion = "RSS") {
  
  n = ncol(y_train)
  L = nrow(y_train)
  
  # Compute residuals of naive
  res_n = y_train[2:L, ] - y_train[1:(L - 1), ]
  
  if (freq == 1) {
    residuals = res_n
    
  } else {
    # compute residuals of seasonal naive
    res_seas = y_train[(freq + 1):L, ] - y_train[1:(L - freq), ]
    
    # Choose which residuals to use based on the criterion
    if (criterion == "seas-test") {
      stop("seas-test criterion not yet implemented")
      
    } else if (criterion == "RSS") {
      # Compute RSS for both methods, for each series
      RSS_n = colSums(res_n^2, na.rm = TRUE)
      RSS_seas = colSums(res_seas^2, na.rm = TRUE)
      is_seas = RSS_seas < RSS_n
      
    } else {
      stop("Input error: criterion must be either 'RSS' or 'seas-test'")
    }
    
    # No series is seasonal
    if (sum(is_seas) == 0) {
      residuals = res_n
      
    # All series are seasonal
    } else if (sum(is_seas) == n) {
      residuals = res_seas
      
    # Some series are seasonal, some are not (use criterion)
    } else  {
      residuals = matrix(NA, nrow = L - freq, ncol = n)
      residuals[, is_seas] = res_seas[, is_seas]
      residuals[, !is_seas] = res_n[freq:(L - 1), !is_seas]
    }
  }
  
  Sigma_naive = schaferStrimmer_cov(residuals)$shrink_cov
  
  return(Sigma_naive)
}


#' Optimize Degrees of Freedom (nu) via LOO Cross-Validation
#'
#' @param res Matrix of residuals (n_obs x n_var).
#' @param prior_mean The prior mean covariance matrix (n_var x n_var).
#' @param trim Fraction of observations to trim (0 to 1). Default is 0.1 (10%).
#'
#' @details
#' TODO: check these details
#' \strong{1. Model Parameterization:}
#' The function uses the "Bayesian" parameterization of the Multivariate Student-t. 
#' Instead of the standard scale matrix $\Sigma$, it operates on the posterior sum-of-squares 
#' matrix $\Psi$.
#' \itemize{
#'   \item Standard t-density: $f(x) \propto |\Sigma|^{-1/2} (1 + \frac{1}{\nu} x^T \Sigma^{-1} x)^{-\frac{\nu+p}{2}}$
#'   \item Implemented t-density: $f(x) \propto |\Psi|^{-1/2} (1 + x^T \Psi^{-1} x)^{-\frac{\nu+p}{2}}$
#' }
#' These are mathematically equivalent because the code sets $\Psi \approx \nu \Sigma$.
#' This absorbs the terms $\nu^{-p/2}$ (into the determinant) and $1/\nu$ (into the inverse).
#'
#' \strong{2. Efficient LOO Computation:}
#' Computing the LOO score normally requires inverting the matrix $\Psi_{-i}$ for every 
#' observation $i=1 \dots T$. This function uses algebraic shortcuts to compute it 
#' using only the full matrix $\Psi$:
#'
#' \emph{Step A: Determinant Adjustment (Matrix Determinant Lemma)}
#' \[ |\Psi_{-i}| = |\Psi| (1 - h_i) \]
#'
#' \emph{Step B: Quadratic Form Adjustment (Sherman-Morrison Formula)}
#' The LOO quadratic form $Q_{-i} = r_i^T \Psi_{-i}^{-1} r_i$ simplifies analytically to:
#' \[ 1 + Q_{-i} = \frac{1}{1 - h_i} \]
#' where $h_i = r_i^T \Psi^{-1} r_i$ is the leverage of the $i$-th observation.
#'
#' \strong{3. Final Log-Likelihood:}
#' Combining these, the LOO log-likelihood for observation $i$ becomes:
#' \[ \text{LL}_i = \text{Const} + \frac{\nu + T - 1}{2} \log(1 - h_i) \]
#' This formula (implemented in step 4 of the code) captures both the change in the 
#' determinant volume and the change in the Mahalanobis distance simultaneously.
#'
#' @return A list containing the optimization results:
#' * `optimal_nu`: The optimal degrees of freedom found.
#' * `min_neg_log_score`: The minimum negative log score achieved.
#' * `convergence`: The convergence status of the optimizer.
#' * `time_elapsed`: The time taken for the optimization.
#'
multi_log_score_optimization <- function(res, prior_mean, trim = 0.1) {
  
  n_obs <- nrow(res)
  n_var <- ncol(res)
  
  # Pre-compute cross-product of the residuals (to avoid recomputing inside the loop)
  RRt <- crossprod(res) 
  
  # Define the Objective Function (Negative Log Score to minimize)
  objective_function <- function(nu) {
    
    # Compute posterior Psi
    Psi <- (prior_mean * (nu - n_var - 1)) + RRt

    # Compute Cholesky decomposition 
    Psi_chol <- tryCatch(chol(Psi), error = function(e) return(NULL))
    if (is.null(Psi_chol)) return(Inf) # Return Inf if Psi is not positive definite
    
    # Invert (using Cholesky) 
    inv_Psi <- chol2inv(Psi_chol)
    
    # Calculate log-determinant
    log_det_Psi <- 2 * sum(log(diag(Psi_chol))) 
    
    # Log Gamma terms and determinant
    const_term <- lgamma((nu + n_obs) / 2) - lgamma((nu + n_obs - n_var) / 2) - log_det_Psi / 2
    
    # Compute LOO Terms (using Sherman-Morrison)
    h_i <- rowSums((res %*% inv_Psi) * res)
    if (any(h_i >= 1)) return(Inf)   # If h_i >= 1, the log is undefined
    log_adjustments <- log(1 - h_i) 
    
    # Combine to get vector of LOO Log-Likelihoods per observation
    loo_log_liks <- const_term + ((nu + n_obs - 1) / 2) * log_adjustments
    
    # Calculate how many observations to remove (based on 'trim')
    n_trim <- round(trim * n_obs)
    
    # Sort likelihoods (ascending) and remove the lowest 'n_trim' values
    # Return negative sum because optimizer minimizes 
    return(-sum(sort(loo_log_liks)[(n_trim + 1):n_obs]))
  }
  
  # Optimization Setup
  initial_guess <- n_var + 3
  lower_bound   <- n_var + 2
  upper_bound   <- max(5 * n_var, n_obs)
  
  opts <- list("algorithm" = "NLOPT_LN_BOBYQA", 
               "xtol_rel"  = 1e-5, 
               "maxeval"   = 1000)
  
  # Run Optimization
  start_time <- Sys.time()
  results <- nloptr::nloptr(
    x0     = initial_guess, 
    eval_f = objective_function, 
    lb     = lower_bound, 
    ub     = upper_bound, 
    opts   = opts
  )
  end_time <- Sys.time()
  
  return(list(
    optimal_nu = results$solution,
    min_neg_log_score = results$objective,
    convergence = results$status,
    time_elapsed = end_time - start_time
  ))
}


