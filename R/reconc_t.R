#' Estimate covariance from naive or seasonal-naive residuals
#'
#' Estimates via shrinkage the covariance matrix of the residuals of the naive 
#' or seasonal naive forecasts.
#' If the frequency of the time series is > 1, the function chooses between the
#' two methods series-by-series according to a selection criterion.
#' If the frequency is 1 or is not provided, only the naive residuals are used.
#'
#' @param y_train Multivariate time series object or numeric matrix of historical 
#'   observations, with dimensions `T x n` (rows are time points, columns are series).
#' @param freq Positive integer seasonal frequency (optional).
#'   If not provided and `y_train` is a multivariate time series, the frequency 
#'   of the data is used. 
#' @param criterion Character string used when `freq > 1` to choose residuals.
#'   Supported values are `"RSS"` (default) and `"seas-test"`.
#'   `"RSS"` chooses the method with the lower residual sum of squares (RSS), 
#'   while `"seas-test"` uses a statistical test for seasonality (requires `forecast` package).
#'
#' @return A numeric `n x n` shrinkage covariance matrix estimated with
#'   [schaferStrimmer_cov()].
#'
#' @seealso [schaferStrimmer_cov()], [reconc_t()]
#'
#' @keywords internal
#' @export
.compute_naive_cov = function(y_train, freq = NULL, criterion = "RSS") {
  
  n <- ncol(y_train)
  L <- nrow(y_train)
  
  # Compute residuals of naive
  res_n <- y_train[2:L, ] - y_train[1:(L - 1), ]
  
  # If freq is not provided and the data is a multivariate time series, 
  # set freq to the frequency of the data
  if (stats::is.mts(y_train)) {
    if (is.null(freq)) {
      freq <- stats::frequency(y_train)
    }
  }

  if (is.null(freq) || freq == 1) {
    residuals <- res_n
  } else {
    # compute residuals of seasonal naive
    res_seas <- y_train[(freq + 1):L, ] - y_train[1:(L - freq), ]

    if (is.null(criterion)) {
      criterion <- "RSS"
    }

    # Choose which residuals to use based on the criterion
    if (criterion == "seas-test") {
      # check if forecast package is installed
      if (!requireNamespace("forecast", quietly = TRUE)) {
        stop("Package 'forecast' is required for criterion = 'seas-test'. 
              Please install it using install.packages('forecast') or use criterion = 'RSS'."
        )
      }
      # Seasonality test for each time series
      is_seas = as.logical(apply(stats::ts(y_train, frequency = freq), 2, 
                                 function(col) forecast::nsdiffs(stats::ts(col,frequency = freq))))
      
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


#' Optimize degrees of freedom (nu) via LOO Cross-Validation
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
#' For each LOO step, the residuals are assumed to follow a
#' multivariate t-distribution. The density is expressed directly as a function of the
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


#' t-Rec: reconciliation via conditioning with uncertain covariance via multivariate t-distribution
#'
#' Reconciles base forecasts in a hierarchy by conditioning on the hierarchical
#' constraints, specified by the aggregation matrix A.
#' The base forecasts are assumed to be jointly Gaussian, conditionally on the
#' covariance matrix of the forecast errors. 
#' A Bayesian approach is adopted to account for the uncertainty of the covariance matrix. 
#' An Inverse-Wishart prior is specified on the covariance matrix, 
#' leading to a multivariate t-distribution for the base forecasts.
#' The reconciliation via conditioning is in closed-form, yielding a multivariate t 
#' reconciled distribution.
#'
#' @param A Matrix (n_upp x n_bott) defining the hierarchy (u = Ab).
#' @param base_fc_mean Vector of base forecasts (length n = n_upp + n_bott).
#' @param y_train mts (or matrix) of historical training data (T x n) used for setting prior parameters.
#' @param residuals Matrix (T x n) of base forecast residuals.
#' @param ... Additional arguments for advanced usage: see details.
#' @param return_upper Logical; if TRUE, also returns parameters for the upper-level reconciled distribution.
#' @param return_parameters Logical; if TRUE, also returns prior and posterior parameters.
#'
#' @details
#' \strong{Standard usage.}
#' 
#' The standard workflow for this function is to provide the in-sample \code{residuals}
#' and the historical training data \code{y_train}.
#' The parameters of the Inverse-Wishart prior distribution of the covariance matrix 
#' are set as follows: 
#' \itemize{
#'   \item Prior scale matrix (Psi): set as the covariance of the residuals of naive (or seasonal naive,
#'         a criterion is used to choose between the two) forecasts computed on \code{y_train}.
#'   \item Prior degrees of freedom (nu): estimated via Bayesian Leave-One-Out
#'         Cross-Validation (LOOCV) to maximize out-of-sample performance.
#' }
#' The posterior distribution of the covariance matrix is still Inverse-Wishart.
#' The parameters of the posterior are computed in closed-form using the sample 
#' covariance of the provided \code{residuals}.
#'
#' \strong{Advanced Options.}
#' 
#' Users can bypass the automated estimation by specifying:
#' \enumerate{
#'   \item \code{prior}: a list with entries 'nu' and 'Psi'.
#'         This skips the LOOCV step for \eqn{\nu_0} and the covariance estimation from \code{y_train}.
#'         It requires \code{residuals} to compute the posterior.
#'   \item \code{posterior}: a list with entries 'nu' and 'Psi'.
#'         This skips all internal estimation and updating logic.
#' }
#' Moreover, users can specify:
#' \itemize{
#'   \item \code{freq}: positive integer, used as frequency of data for the seasonal naive forecast in the specification of the prior scale matrix.
#'         By default, if \code{y_train} is a multivariate time series, the frequency of the data is used; otherwise, it is set to 1 (no seasonality).
#'   \item \code{criterion}: either 'RSS' (default) or 'seas-test', specifying which criterior is used to choose between 
#'                           the naive and seasonal naive forecasts for the specification of the prior scale matrix. 
#'                           'RSS' computes the residual sum of squares for both methods and chooses the one with lower RSS, 
#'                           while 'seas-test' uses a statistical test for seasonality 
#'                           (currently implemented using the number of seasonal differences suggested by the `forecast` package, 
#'                           which must be installed). 
#' }
#' 
#' \strong{Reconciled distribution.}
#' 
#' The reconciled distribution is a multivariate t-distribution, 
#' specified by a vector of means, a scale matrix, and a number of degrees of freedom.
#' These parameters are computed in closed-form.
#' By default, only the parameters of the reconciled distribution for the bottom-level 
#' series are returned. See examples.
#'
#'
#' @references
#' Carrara, C., Corani, G., Azzimonti, D., & Zambon, L. (2025). Modeling the uncertainty on the covariance
#' matrix for probabilistic forecast reconciliation. arXiv preprint arXiv:2506.19554.
#' \url{https://arxiv.org/abs/2506.19554}
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{bottom_rec_mean}: reconciled bottom-level mean.
#'   \item \code{bottom_rec_scale_matrix}: reconciled bottom-level scale matrix.
#'   \item \code{bottom_rec_df}: reconciled degrees of freedom.
#' }
#' If \code{return_upper} is TRUE, also returns:
#' \itemize{
#'  \item \code{upper_rec_mean}: reconciled upper-level mean.
#'  \item \code{upper_rec_scale_matrix}: reconciled upper-level scale matrix.
#'  \item \code{upper_rec_df}: reconciled upper-level degrees of freedom.
#'  }
#' If \code{return_parameters} is TRUE, also returns:
#' \itemize{
#'   \item \code{prior_nu}: prior degrees of freedom.
#'   \item \code{prior_Psi}: prior scale matrix.
#'   \item \code{posterior_nu}: posterior degrees of freedom.
#'   \item \code{posterior_Psi}: posterior scale matrix.
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
#'   base_fc_mean <- c(fc_upper, fc1, fc2)
#'
#'   # Residuals and training data (n_obs x n matrices, columns in same order as base_fc_mean)
#'   res <- cbind(residuals(fit_upper), residuals(fit1), residuals(fit2))
#'   y_train <- cbind(y_upper, y1, y2)
#'
#'   # --- 1) Generate joint reconciled samples ---
#'   result <- reconc_t(A, base_fc_mean, y_train = y_train, residuals = res)
#'
#'   # Sample from the reconciled bottom-level t-distribution
#'   n_samples <- 2000
#'   L_chol <- t(chol(result$bottom_rec_scale_matrix))
#'   z <- matrix(rt(ncol(A) * n_samples, df = result$bottom_rec_df), nrow = ncol(A))
#'   bottom_samples <- result$bottom_rec_mean + L_chol %*% z  # 2 x n_samples
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
#'   result2 <- reconc_t(A, base_fc_mean, y_train = y_train,
#'                       residuals = res, return_upper = TRUE)
#'
#'   alpha <- 0.05
#'   # Bottom series intervals
#'   for (i in seq_len(ncol(A))) {
#'     s_i <- sqrt(result2$bottom_rec_scale_matrix[i, i])
#'     lo  <- result2$bottom_rec_mean[i] + s_i * qt(alpha / 2,     df = result2$bottom_rec_df)
#'     hi  <- result2$bottom_rec_mean[i] + s_i * qt(1 - alpha / 2, df = result2$bottom_rec_df)
#'     cat(sprintf("Bottom %d: 95%% PI = [%.3f, %.3f]\n", i, lo, hi))
#'   }
#'   # Upper series interval
#'   s_u <- sqrt(result2$upper_rec_scale_matrix[1, 1])
#'   lo  <- result2$upper_rec_mean[1] + s_u * qt(alpha / 2,     df = result2$upper_rec_df)
#'   hi  <- result2$upper_rec_mean[1] + s_u * qt(1 - alpha / 2, df = result2$upper_rec_df)
#'   cat(sprintf("Upper:    95%% PI = [%.3f, %.3f]\n", lo, hi))
#' }
#' }
#' 
#' @seealso 
#' [reconc_gaussian()]
#'
#' @export
reconc_t <- function(A,
                     base_fc_mean,
                     y_train = NULL,
                     residuals = NULL,
                     ...,
                     return_upper = FALSE,
                     return_parameters = FALSE) {
  
  add_args <- list(...)
  unused_names <- setdiff(names(add_args), c("prior", "posterior", "freq", "criterion"))
  if (length(unused_names) > 0) {
    warning(paste("The following additional arguments are not used:", paste(unused_names, collapse = ", ")))
  }
  
  prior <- add_args$prior
  posterior <- add_args$posterior
  freq <- add_args$freq
  criterion <- add_args$criterion
  .check_input_t(A, base_fc_mean, y_train, residuals, ...)  

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
    n <- length(base_fc_mean) # number of series
    
    # Compute sample covariance of the residuals
    Samp_cov <- crossprod(residuals) / nrow(residuals) 
    # Shrink to diagonal for numerical reasons
    Samp_cov <- (1 - .L_SHRINK_RECONC_T) * Samp_cov + .L_SHRINK_RECONC_T * diag(diag(Samp_cov)) 
    
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
      # Compute the covariance residuals of the naive or seasonal naive forecasts on training data
      cov_naive <- .compute_naive_cov(y_train, freq = freq, criterion = criterion)
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
    A = A, base_fc_mean = base_fc_mean, Psi_post = Psi_post, nu_post = nu_post,
    return_upper = return_upper)

  if (return_parameters) {
    if(!is.null(posterior)) {
      warning("Prior parameters are not returned when 'posterior' is provided, 
              as the prior is not used in this case.")
    } else {
      out$prior_nu <- nu_prior
      out$prior_Psi <- Psi_prior
    }
    out$posterior_nu <- nu_post
    out$posterior_Psi <- Psi_post
  }

  return(out)
}

#' Core reconciliation via multivariate t-distribution.
#'
#' Internal function that performs the core reconciliation logic for the t-distribution based
#' reconciliation method. This function assumes an uncertain covariance matrix with an Inverse-Wishart prior.
#'
#' @param A Matrix (n_upper x n_bottom) defining the hierarchy where upper = A %*% bottom.
#' @param base_fc_mean Vector of length (n_upper + n_bottom) containing the base forecast means
#'   for both upper and bottom levels (upper first, then bottom).
#' @param Psi_post Scale matrix (n_upper + n_bottom x n_upper + n_bottom) of the posterior
#'   multivariate t-distribution.
#' @param nu_post Degrees of freedom of the posterior multivariate t-distribution.
#' @param return_upper Logical, whether to return the reconciled parameters for the upper variables (default is FALSE).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{bottom_rec_mean}: reconciled bottom-level mean.
#'   \item \code{bottom_rec_scale_matrix}: reconciled bottom-level scale matrix.
#'   \item \code{bottom_rec_df}: reconciled degrees of freedom.
#' }
#' If \code{return_upper} is TRUE, also returns:
#' \itemize{
#'  \item \code{upper_rec_mean}: reconciled upper-level mean.
#'  \item \code{upper_rec_scale_matrix}: reconciled upper-level scale matrix.
#'  \item \code{upper_rec_df}: reconciled upper-level degrees of freedom.
#'  }
#'
#' @keywords internal
#' @export
.core_reconc_t <- function(A, base_fc_mean, Psi_post, nu_post, return_upper = FALSE,
                           return_parameters = FALSE) {
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
  u_hat <- base_fc_mean[idx_u]
  b_hat <- base_fc_mean[idx_b]

  # Compute Q = Psi_U - (Psi_UB %*% t(A)) - (A %*% t(Psi_UB)) + (A %*% Psi_B %*% t(A))
  Psi_UB_At <- tcrossprod(Psi_UB, A)
  A_Psi_B <- A %*% Psi_B
  Q <- Psi_U - Psi_UB_At - t(Psi_UB_At) + tcrossprod(A_Psi_B, A)

  # Invert Q using Cholesky (should be p.d.)
  Q_chol <- tryCatch(chol(Q), error = function(e) NULL)
  # If Cholesky fails, stop
  if (is.null(Q_chol)) {
    stop("A numerical error occured during a Cholesky decomposition. 
         This may be due to redundant aggregated series.")
  }
  inv_Q <- chol2inv(Q_chol)

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
    bottom_rec_mean = as.vector(b_tilde),
    bottom_rec_scale_matrix = Sigma_tilde_B,
    bottom_rec_df = nu_tilde
  )
  if (return_upper) {
    # Compute the parameters of the uppers using closure property
    u_tilde <- A %*% b_tilde
    Sigma_tilde_U <- A %*% Sigma_tilde_B %*% t(A)
    out$upper_rec_mean <- as.vector(u_tilde)
    out$upper_rec_scale_matrix <- Sigma_tilde_U
    out$upper_rec_df <- nu_tilde
  }

  return(out)
}
