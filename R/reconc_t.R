# Function for estimating the covariance matrix of the residuals of the naive or seasonal naive forecasts
# For each series, choose by using some criterion (RSS or test of seasonality)
# TODO: implement statistical test; leave dependence as suggested; default: choose using RSS
compute_naive_cov <- function(y_train, freq = 1, criterion = "RSS") {
  n <- ncol(y_train)
  L <- nrow(y_train)

  # Compute residuals of naive
  res_n <- y_train[2:L, ] - y_train[1:(L - 1), ]

  if (freq == 1) {
    residuals <- res_n
  } else {
    # compute residuals of seasonal naive
    res_seas <- y_train[(freq + 1):L, ] - y_train[1:(L - freq), ]

    # Choose which residuals to use based on the criterion
    if (criterion == "seas-test") {
      stop("seas-test criterion not yet implemented")
    } else if (criterion == "RSS") {
      # Compute RSS for both methods, for each series
      RSS_n <- colSums(res_n^2, na.rm = TRUE)
      RSS_seas <- colSums(res_seas^2, na.rm = TRUE)
      is_seas <- RSS_seas < RSS_n
    } else {
      stop("Input error: criterion must be either 'RSS' or 'seas-test'")
    }

    # No series is seasonal
    if (sum(is_seas) == 0) {
      residuals <- res_n

      # All series are seasonal
    } else if (sum(is_seas) == n) {
      residuals <- res_seas

      # Some series are seasonal, some are not (use criterion)
    } else {
      residuals <- matrix(NA, nrow = L - freq, ncol = n)
      residuals[, is_seas] <- res_seas[, is_seas]
      residuals[, !is_seas] <- res_n[freq:(L - 1), !is_seas]
    }
  }

  Sigma_naive <- schaferStrimmer_cov(residuals)$shrink_cov

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
#'
#' \strong{Leave-One-Out (LOO) Cross-Validation:}
#' This function estimates the optimal degrees of freedom \eqn{\nu} by maximizing the
#' out-of-sample predictive performance. This is achieved by computing the
#' log-density of each held-out observation \eqn{\mathbf{r}_i} given the remaining
#' data \eqn{\mathbf{R}_{-i}}. The total objective function is the sum of these
#' predictive log-densities:
#' \deqn{\mathcal{L}(\nu) = \sum_{i=1}^T \log f(\mathbf{r}_i | \mathbf{R}_{-i}, \nu)}
#'
#' \strong{The Log-Density Function:}
#' For each LOO step, the residuals are assumed to follow a Multivariate
#' Student-t distribution. The density is expressed directly as a function of the
#' posterior sum-of-squares matrix \eqn{\Psi}, where \eqn{\Psi} scales implicitly with \eqn{\nu}:
#' \deqn{f(\mathbf{r}_i | \Psi, \nu) = \frac{\Gamma(\frac{\nu + T}{2})}{\Gamma(\frac{\nu + T - p}{2}) \pi^{p/2}} |\Psi|^{-1/2} \left( 1 + \mathbf{r}_i^\top \Psi^{-1} \mathbf{r}_i \right)^{-\frac{\nu + T}{2}}}
#' In the code, \eqn{\Psi} is constructed as:
#' \deqn{\Psi = (\nu - p - 1)\bar{\Sigma}_{prior} + \mathbf{R}^\top\mathbf{R}}
#' By using this formulation, the standard scaling factors \eqn{1/\nu} and \eqn{\nu^{-p/2}}
#' are absorbed into the matrix inverse and determinant, respectively.
#'
#' \strong{Efficient Computation via Sherman-Morrison:}
#' Rather than recomputing \eqn{\Psi_{-i}} and its inverse \eqn{T} times, the function uses
#' the full-sample matrix \eqn{\Psi} and adjusts it using the leverage
#' \eqn{h_i = \mathbf{r}_i^\top \Psi^{-1} \mathbf{r}_i}.
#'
#' Through the Matrix Determinant Lemma and the Sherman-Morrison formula, the
#' internal term \eqn{(1 + \mathbf{r}_i^\top \Psi_{-i}^{-1} \mathbf{r}_i)} simplifies to \eqn{(1 - h_i)^{-1}}.
#' The final log-density contribution used in the code is:
#' \deqn{\log f_i \propto \textrm{const} - \frac{1}{2}\log|\Psi| + \frac{\nu + T - 1}{2} \log(1 - h_i)}
#'
#' @return A list containing the optimization results:
#' * `optimal_nu`: The optimal degrees of freedom found.
#' * `min_neg_log_score`: The minimum negative log score achieved.
#' * `convergence`: The convergence status of the optimizer.
#' * `time_elapsed`: The time taken for the optimization.
#'
#' @import nloptr
#'
#' @export
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
    Psi_chol <- tryCatch(chol(Psi), error = function(e) {
      return(NULL)
    })
    if (is.null(Psi_chol)) {
      return(Inf)
    } # Return Inf if Psi is not positive definite

    # Invert (using Cholesky)
    inv_Psi <- chol2inv(Psi_chol)

    # Calculate log-determinant
    log_det_Psi <- 2 * sum(log(diag(Psi_chol)))

    # Log Gamma terms and determinant
    const_term <- lgamma((nu + n_obs) / 2) - lgamma((nu + n_obs - n_var) / 2) - log_det_Psi / 2

    # Compute LOO Terms (using Sherman-Morrison)
    h_i <- rowSums((res %*% inv_Psi) * res)
    if (any(h_i >= 1)) {
      return(Inf)
    } # If h_i >= 1, the log is undefined
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
  lower_bound <- n_var + 2
  upper_bound <- max(5 * n_var, n_obs)

  opts <- list(
    "algorithm" = "NLOPT_LN_BOBYQA",
    "xtol_rel" = 1e-5,
    "maxeval" = 1000
  )

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


#' t-Rec: Reconciliation via Conditioning with uncertain covariance via Multivariate Student-t
#'
#' Reconciles base forecasts in a hierarchy by conditioning on the hierarchical
#' constraints, specified by the aggregation matrix A.
#' The base forecasts are assumed to be jointly Gaussian, conditionally on the
#  covariance matrix of the forecast errors. To account for uncertainty in the
#' covariance matrix, a Bayesian approach is adopted using an Inverse-Wishart prior,
#' leading to a Multivariate Student-t distribution for the base forecasts.
#' The reconciliation is in closed-form, yielding a multivariate Student-t reconciled distribution.
#'
#' @param A Matrix (n_upp x n_bott) defining the hierarchy (u = Ab).
#' @param point_fc Vector of base forecasts (length n = n_upp + n_bott).
#' @param y_train Matrix of historical training data (T x n) used for setting prior parameters.
#' @param residuals Optional matrix (T x n) of base forecast residuals.
#' @param freq Seasonal frequency for naive covariance estimation (default is 1).
#' @param prior Optional list containing 'nu' and 'Psi' (prior parameters).
#' @param posterior Optional list containing 'nu' and 'Psi' (posterior parameters).
#' @param l_shr Shrinkage intensity (0 to 1) for stabilizing the sample covariance matrix (default 1e-4).
#' @param return_uppers Logical; if TRUE, also returns parameters for the upper level reconciled distribution.
#' @param return_parameters Logical; if TRUE, returns internal parameters like C and posterior nu.
#'
#' @details
#' \strong{Standard Usage and Parameter Estimation:}
#' The standard workflow for this function is to provide the in-sample \code{residuals}
#' and the historical training data \code{y_train}.
#' \itemize{
#'   \item \strong{Prior Scale (\eqn{\Psi_0}):} Set as the covariance of the residuals of naive (or seasonal naive,
#'         a criterion is used to choose between the 2) forecasts computed on \code{y_train}.
#'   \item \strong{Prior Degrees of Freedom (\eqn{\nu_0}):} Estimated via Bayesian Leave-One-Out
#'         Cross-Validation (LOOCV) to maximize out-of-sample performance.
#' }
#'
#' \strong{Advanced Options:}
#' Users can bypass the automated estimation by:
#' \enumerate{
#'   \item Directly passing the \code{prior} parameters (requires \code{residuals} to compute the posterior).
#'   \item Directly passing the \code{posterior} parameters, which skips all internal estimation and updating logic.
#' }
#'
#' \strong{The Reconciled Bottom Distribution:}
#' The reconciliation yields a distribution:
#' \deqn{\tilde{\mathbf{b}} \sim t(\hat{\mathbf{b}}_{tilde}, \tilde{\Sigma}_B, \tilde{\nu})}
#' where the reconciled mean is:
#' \deqn{\hat{\mathbf{b}}_{tilde} = \hat{\mathbf{b}} + (\Psi'_{UB}^\top - \Psi'_B A^\top) Q^{-1} (A\hat{\mathbf{b}} - \hat{\mathbf{u}})}
#' and the scale matrix is:
#' \deqn{\tilde{\Sigma}_B = C [\Psi'_B - (\Psi'_{UB}^\top - \Psi'_B A^\top) Q^{-1} (\Psi'_{UB}^\top - \Psi'_B A^\top)^\top]}
#' with scalar \deqn{C = \frac{1 + (A\hat{\mathbf{b}} - \hat{\mathbf{u}})^\top Q^{-1} (A\hat{\mathbf{b}} - \hat{\mathbf{u}})}{\tilde{\nu}}.}
#'
#' @references
#' Carrara, C., Corani, G., Azzimonti, D., & Zambon, L. (2025). Modeling the uncertainty on the covariance
#' matrix for probabilistic forecast reconciliation. arXiv preprint arXiv:2506.19554.
#' \url{https://arxiv.org/abs/2506.19554}
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{bottom_mean}: Reconciled bottom-level mean forecasts.
#'   \item \code{bottom_scale_matrix}: Reconciled bottom-level scale matrix.
#'   \item \code{bottom_df}: Reconciled degrees of freedom.
#' }
#' If \code{return_parameters} is TRUE, also returns:
#' \itemize{
#'   \item \code{prior_nu}: Prior degrees of freedom.
#'   \item \code{posterior_nu}: Posterior degrees of freedom.
#'   \item \code{posterior_Psi}: Posterior scale matrix.
#'   \item \code{C}: Scaling factor for the scale matrix.
#' }
#'
#' @examples
#' 
#' \donttest{
#' library(bayesRecon)
#'
#' if (requireNamespace("forecast", quietly = TRUE)) {
#' 
#'   set.seed(1234)
#'   n_obs <- 100
#' 
#'   # Simulate 2 bottom series from AR(1) processes
#'   y1 <- arima.sim(model = list(ar = 0.8), n = n_obs)
#'   y2 <- arima.sim(model = list(ar = 0.5), n = n_obs)
#'
#'   y_upper <- y1 + y2   # upper series is the sum of the two bottoms
#'   A <- matrix(c(1, 1), nrow = 1)  # Aggregation matrix
#'
#'   # Fit additive ETS models
#'   fit1 <- forecast::ets(y1, additive.only = TRUE)
#'   fit2 <- forecast::ets(y2, additive.only = TRUE)
#'   fit_upper <- forecast::ets(y_upper, additive.only = TRUE)
#'
#'   # Point forecasts (h = 1)
#'   fc_upper <- as.numeric(forecast::forecast(fit_upper, h = 1)$mean)
#'   fc1 <- as.numeric(forecast::forecast(fit1, h = 1)$mean)
#'   fc2 <- as.numeric(forecast::forecast(fit2, h = 1)$mean)
#'   point_fc <- c(fc_upper, fc1, fc2)
#'
#'   # Residuals and training data (n_obs x n matrices, columns in same order as point_fc)
#'   res <- cbind(residuals(fit_upper), residuals(fit1), residuals(fit2))
#'   y_train <- cbind(y_upper, y1, y2)
#'
#'   # --- 1) Generate joint reconciled samples ---
#'   result <- reconc_t(A, point_fc, y_train = y_train, residuals = res)
#'
#'   # Sample from the reconciled bottom-level Student-t distribution
#'   n_samples <- 2000
#'   L_chol <- t(chol(result$bottom_scale_matrix))
#'   z <- matrix(rt(ncol(A) * n_samples, df = result$bottom_df), nrow = ncol(A))
#'   bottom_samples <- result$bottom_mean + L_chol %*% z  # 2 x n_samples
#'
#'   # Aggregate bottom samples to get upper samples
#'   upper_samples <- A %*% bottom_samples              
#'   joint_samples <- rbind(upper_samples, bottom_samples)
#'   rownames(joint_samples) <- c("upper", "bottom_1", "bottom_2")
#'
#'   cat("Reconciled means (from samples):\n")
#'   print(round(rowMeans(joint_samples), 3))
#' 
#'   cat("Reconciled standard deviations (from samples):\n")
#'   print(round(apply(joint_samples, 1, sd), 3))
#'
#'   # --- 2) 95% prediction intervals via t-distribution quantiles ---
#'   result2 <- reconc_t(A, point_fc, y_train = y_train,
#'                       residuals = res, return_uppers = TRUE)
#'
#'   alpha <- 0.05
#'   # Bottom series intervals
#'   for (i in seq_len(ncol(A))) {
#'     s_i <- sqrt(result2$bottom_scale_matrix[i, i])
#'     lo  <- result2$bottom_mean[i] + s_i * qt(alpha / 2,     df = result2$bottom_df)
#'     hi  <- result2$bottom_mean[i] + s_i * qt(1 - alpha / 2, df = result2$bottom_df)
#'     cat(sprintf("Bottom %d: 95%% PI = [%.3f, %.3f]\n", i, lo, hi))
#'   }
#'   # Upper series interval
#'   s_u <- sqrt(result2$upper_scale_matrix[1, 1])
#'   lo  <- result2$upper_mean[1] + s_u * qt(alpha / 2,     df = result2$upper_df)
#'   hi  <- result2$upper_mean[1] + s_u * qt(1 - alpha / 2, df = result2$upper_df)
#'   cat(sprintf("Upper:    95%% PI = [%.3f, %.3f]\n", lo, hi))
#' }
#' }
#'
#' @export
reconc_t <- function(A,
                     point_fc,
                     y_train = NULL,
                     residuals = NULL,
                     freq = 1,
                     prior = NULL,
                     posterior = NULL,
                     l_shr = 1e-4,
                     return_uppers = FALSE,
                     return_parameters = FALSE) {
  
  .check_input_t(A, point_fc, y_train, residuals, freq, prior, posterior, l_shr)

  ##############################################################################
  ### CASE 1 ###
  # If posterior is provided, check if is a list with entries nu and Psi and extract values
  if (!is.null(posterior)) {
    nu_post <- posterior$nu
    Psi_post <- posterior$Psi
    ### CASE 2 ###
    # If posterior not provided, first check that residuals are provided
  } else {
    L <- nrow(residuals) # number of residual samples (i.e., training length)
    n <- length(point_fc) # number of series
    # TODO: implement fallback
    Samp_cov <- crossprod(residuals) / nrow(residuals) # sample covariance of the residuals
    Samp_cov <- (1 - l_shr) * Samp_cov + l_shr * diag(diag(Samp_cov)) # apply shrinkage to stabilize

    ### CASE 2a ###
    # If prior is provided, check if is a list with entries nu and Psi and extract values
    if (!is.null(prior)) {
      nu_prior <- prior$nu
      Psi_prior <- prior$Psi
      ### CASE 2b ###
      # If prior not provided:
      # - compute Psi using the (shrinked) covariance matrix of the residuals of the naive
      #   or seasonal naive forecasts
      # - set nu using LOOCV
    } else {
      # Compute the covariance residuals of the (seasonal) naive forecasts on training data
      cov_naive <- compute_naive_cov(y_train, freq)
      # Use it to set the prior
      bayesian_LOO <- multi_log_score_optimization(residuals, cov_naive)
      nu_prior <- bayesian_LOO$optimal_nu
      Psi_prior <- (nu_prior - n - 1) * cov_naive
    }
    # Compute posterior parameters
    Psi_post <- Psi_prior + L * Samp_cov
    nu_post <- nu_prior + nrow(residuals)
  }

  ##############################################################################
  # Reconcile via conditioning the t-distribution in closed form

  out <- .core_reconc_t(
    A = A, point_fc = point_fc, Psi_post = Psi_post, nu_post = nu_post,
    return_uppers = return_uppers, return_parameters = return_parameters, suppress_warnings = FALSE
  )

  return(out)
}

#' Core Reconciliation via Multivariate Student-t Distribution.
#'
#' Internal function that performs the core reconciliation logic for the t-distribution based
#' reconciliation method. This function assumes an uncertain covariance matrix with an Inverse-Wishart prior.
#'
#' @param A Matrix (n_upper x n_bottom) defining the hierarchy where upper = A %*% bottom.
#' @param point_fc Vector of length (n_upper + n_bottom) containing the base forecast means
#'   for both upper and bottom levels (upper first, then bottom).
#' @param Psi_post Scale matrix (n_upper + n_bottom x n_upper + n_bottom) of the posterior
#'   Student-t distribution.
#' @param nu_post Degrees of freedom of the posterior Student-t distribution.
#' @param return_uppers Logical. If TRUE, also returns parameters for the upper level reconciled
#'   distribution. Default is FALSE.
#' @param return_parameters Logical. If TRUE, returns internal parameters (C matrix, posterior nu, etc.)
#'   for debugging or advanced use. Default is FALSE.
#' @param suppress_warnings Logical. If TRUE, suppresses warnings about numerical issues. Default is FALSE.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `bottom_mean`: Reconciled mean vector for bottom level.
#'     \item `bottom_scale_matrix`: Reconciled scale matrix for bottom level.
#'     \item `bottom_df`: Reconciled degrees of freedom for bottom level.
#'     \item `upper_mean`: (optional) Reconciled mean vector for upper level.
#'     \item `upper_scale_matrix`: (optional) Reconciled scale matrix for upper level.
#'     \item `upper_df`: (optional) Reconciled degrees of freedom for upper level.
#'     \item `posterior_nu`: (optional) Posterior degrees of freedom.
#'     \item `posterior_Psi`: (optional) Posterior scale matrix.
#'     \item `C`: (optional) Scaling factor used in reconciliation.
#'   }
#'
#' @keywords internal
#' @noRd
#' @export
.core_reconc_t <- function(A, point_fc, Psi_post, nu_post, return_uppers = FALSE,
                           return_parameters = FALSE, suppress_warnings = FALSE) {
  # Indices for Upper and Bottom
  k <- nrow(A)
  m <- ncol(A)
  idx_u <- 1:k
  idx_b <- (k + 1):(k + m)

  # Extract Psi blocks
  Psi_U <- Psi_post[idx_u, idx_u, drop = FALSE]
  Psi_B <- Psi_post[idx_b, idx_b, drop = FALSE]
  Psi_UB <- Psi_post[idx_u, idx_b, drop = FALSE]

  # Extract means
  u_hat <- point_fc[idx_u]
  b_hat <- point_fc[idx_b]

  # Compute Q = Psi_U - (Psi_UB %*% t(A)) - (A %*% t(Psi_UB)) + (A %*% Psi_B %*% t(A))
  Psi_UB_At <- tcrossprod(Psi_UB, A)
  A_Psi_B <- A %*% Psi_B
  Q <- Psi_U - Psi_UB_At - t(Psi_UB_At) + tcrossprod(A_Psi_B, A)

  # Invert Q using Cholesky (should be p.d.)
  Q_chol <- tryCatch(chol(Q), error = function(e) NULL)
  if (is.null(Q_chol)) {
    # Fallback to standard solve if Cholesky fails (numerical issues)
    if (!suppress_warnings) {
      warning("Cholesky decomposition of Q failed; using standard inversion.")
    }
    inv_Q <- solve(Q)
  } else {
    inv_Q <- chol2inv(Q_chol)
  }

  # Incoherence
  delta <- (A %*% b_hat) - u_hat

  # Lambda = Psi_UB^T - Psi_B A^T
  Lambda <- t(Psi_UB) - tcrossprod(Psi_B, A)
  Lambda_invQ <- Lambda %*% inv_Q

  # Compute b_tilde = b_hat + Lambda * Q^{-1} * delta
  b_tilde <- b_hat + (Lambda_invQ %*% delta)

  # Compute nu_tilde = nu' - n_b + 1
  nu_tilde <- nu_post - m + 1

  # Compute C = (1 + delta^T Q^{-1} delta) / nu_tilde
  mahalanobis_term <- sum(delta * (inv_Q %*% delta)) # efficient x^T A x
  C <- (1 + mahalanobis_term) / nu_tilde

  # Compute Sigma_B_tilde = C * [ Psi_B - Lambda * Q^{-1} * Lambda^T ]
  Sigma_tilde_B <- as.numeric(C) * (Psi_B - tcrossprod(Lambda_invQ, Lambda))

  # Prepare output
  out <- list(
    bottom_mean = as.vector(b_tilde),
    bottom_scale_matrix = Sigma_tilde_B,
    bottom_df = nu_tilde
  )
  if (return_uppers) {
    # Compute the parameters of the uppers using closure property
    u_tilde <- A %*% b_tilde
    Sigma_tilde_U <- A %*% Sigma_tilde_B %*% t(A)
    out$upper_mean <- as.vector(u_tilde)
    out$upper_scale_matrix <- Sigma_tilde_U
    out$upper_df <- nu_tilde
  }
  if (return_parameters) {
    out$posterior_nu <- nu_post
    out$posterior_Psi <- Psi_post
    out$C <- C
  }

  return(out)
}
