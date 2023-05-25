###############################################################################
# Reconciliation with Bottom-Up Importance Sampling (BUIS)

.distr_sample <- function(params, distr_, n) {
  switch(
    distr_,
    "gaussian" = {
      samples = stats::rnorm(n=n, mean = params[[1]], sd = params[[2]]) },
    "poisson"  = {
      samples = stats::rpois(n=n, lambda = params[[1]]) },
    "negbin"   = {
      samples = stats::rnbinom(n=n, size = params[[1]], prob = params[[2]]) },
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
    "negbin"   = {
      pmf = stats::dnbinom(x=x, size = params[[1]], prob = params[[2]]) },
  )
  return(pmf)
}
.emp_pmf <- function(l, density_samples) {
  empirical_pmf = sapply(0:max(density_samples), function(i)
    sum(density_samples == i) / length(density_samples))
  w = sapply(l, function(i) empirical_pmf[i + 1])
  return(w)
}
.fix_weights <- function(w) {
  w[is.na(w)] = 0
  if (sum(w) == 0) {
    w = w + 1
    warning("WARNING: all IS weights are zero, increase sample size or check your forecasts.")
  }
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
  } else if (in_type_ == "params") {
    w = .distr_pmf(b, u, distr_)
  }
  w = .fix_weights(w)
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
#' The parameter `base_forecast` is a list containing n elements that depend on
#' the options `in_type` and `distr`.
#'
#' If `in_type`='samples', each element of `base_forecast` is a vector containing samples from the base forecast distribution.
#'
#' If `in_type`='params', each element of `base_forecast` is a vector containing the estimated:
#'
#' * mean and sd for the Gaussian base forecast, see \link[stats]{Normal}, if `distr`='gaussian';
#' * lambda for the Poisson base forecast, see \link[stats]{Poisson}, if `distr`='poisson';
#' * size and probability of success for the negative binomial base forecast, see \link[stats]{NegBinomial}, if `distr`='nbinom'.
#'
#' The order of the `base_forecast` list is given by the order of the time series in the summing matrix.
#'
#' @param S summing matrix (n x n_bottom).
#' @param base_forecasts a list containing the base_forecasts, see details.
#' @param in_type a string with two possible values:
#'
#' * 'samples' if the base forecasts are in the form of samples;
#' * 'params'  if the base forecasts are in the form of estimated parameters.
#'
#' @param distr a string describing the type of base forecasts:
#'
#' * 'continuous' or 'discrete' if `in_type`='samples';
#' * 'gaussian', 'poisson' or 'nbinom' if `in_type`='params'.
#'
#' @param num_samples number of samples drawn from the reconciled distribution.
#' @param seed seed for reproducibility.
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
#'print(c(analytic_rec$upper_reconciled_mean, analytic_rec$bottom_reconciled_mean))
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
#' Zambon, L., Azzimonti, D. & Corani, G. (2022). *Efficient probabilistic reconciliation of forecasts for real-valued and count time series*. \doi{10.48550/arXiv.2210.02286}.
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
                   seed = NULL) {
  set.seed(seed)

  # Ensure that data inputs are valid
  .check_input(S, base_forecasts, in_type, distr)
  if (!is.list(distr)) {
    distr = rep(list(distr), nrow(S))
  }

  # Split bottoms, uppers
  split_hierarchy.res = .split_hierarchy(S, base_forecasts)
  A = split_hierarchy.res$A
  upper_base_forecasts = split_hierarchy.res$upper
  bottom_base_forecasts = split_hierarchy.res$bottom

  # H, G
  get_HG.res = .get_HG(A, upper_base_forecasts, distr[split_hierarchy.res$upper_idxs])
  H = get_HG.res$H
  upper_base_forecasts_H = get_HG.res$Hv
  G = get_HG.res$G
  upper_base_forecasts_G = get_HG.res$Gv

  # Reconciliation using BUIS
  n_upper = nrow(A)
  n_bottom = ncol(A)
  # 1. Bottom samples
  B = list()
  for (bi in 1:n_bottom) {
    if (in_type == "samples") {
      B[[bi]] = unlist(bottom_base_forecasts[[bi]])
    } else if (in_type == "params") {
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
      in_type_ = in_type,
      distr_ = get_HG.res$Hdistr[[hi]]
    )
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
        in_type_ = in_type,
        distr_ = get_HG.res$Gdistr[[gi]]
      )
    }
    B = .resample(B, weights)
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

###############################################################################
.check_cov <- function(cov_matrix) {
  # Check if the matrix is square
  if (!is.matrix(cov_matrix) || nrow(cov_matrix) != ncol(cov_matrix)) {
    stop("base_forecasts.Sigma not square")
  }
  # Check if the matrix is positive semi-definite
  eigen_values <- eigen(cov_matrix, symmetric = TRUE)$values
  if (any(eigen_values <= 0)) {
    stop("base_forecasts.Sigma not positive semi-definite")
  }
  # Check if the matrix is symmetric
  if (!isSymmetric(cov_matrix)) {
    stop("base_forecasts.Sigma not symmetric")
  }
  # Check if the diagonal elements are non-negative
  if (any(diag(cov_matrix) < 0)) {
    stop("base_forecasts.Sigma, diagonal elements are non-positive")
  }
  # If all checks pass, return TRUE
  return(TRUE)
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
#'
#' @return A list containing the reconciled forecasts. The list has the following named elements:
#'
#' * `bottom_reconciled_mean`: reconciled mean for the bottom forecasts;
#' * `bottom_reconciled_covariance`: reconciled covariance for the bottom forecasts;
#' * `upper_reconciled_mean`: reconciled mean for the upper forecasts;
#' * `upper_reconciled_covariance`: reconciled covariance for the upper forecasts.
#'
#' @examples
#'
#'library(bayesRecon)
#'
#'# Create a minimal hierarchy with 2 bottom and 1 upper variable
#'rec_mat <- get_reconc_matrices(agg_levels=c(1,2), h=2)
#'S <- rec_mat$S
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
#'bottom_means <- analytic_rec$bottom_reconciled_mean
#'upper_means  <- analytic_rec$upper_reconciled_mean
#'bottom_cov   <- analytic_rec$bottom_reconciled_covariance
#'upper_cov    <- analytic_rec$upper_reconciled_covariance
#'
#' @references
#' Corani, G., Azzimonti, D., Augusto, J.P.S.C., Zaffalon, M. (2021). *Probabilistic Reconciliation of Hierarchical Forecast via Bayes’ Rule*. In: Hutter, F., Kersting, K., Lijffijt, J., Valera, I. (eds) Machine Learning and Knowledge Discovery in Databases. ECML PKDD 2020. Lecture Notes in Computer Science(), vol 12459. Springer, Cham. \doi{10.1007/978-3-030-67664-3_13}.
#'
#' Zambon, L., Agosto, A., Giudici, P., Corani, G. (2023). *Properties of the reconciled distributions for Gaussian and count forecasts*. \doi{10.48550/arXiv.2303.15135}.
#'
#'
#' @seealso [reconc_BUIS()]
#'
#' @export
reconc_gaussian <- function(S, base_forecasts.mu,
                            base_forecasts.Sigma) {
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
  mu_u_tilde = mu_u + (Sigma_u - Sigma_ub %*% t(A)) %*% invQ %*% (A %*% mu_b - mu_u)
  Sigma_b_tilde = Sigma_b - (t(Sigma_ub) - Sigma_b %*% t(A)) %*% invQ %*% t(t(Sigma_ub) - Sigma_b %*% t(A))
  Sigma_u_tilde = Sigma_u - (Sigma_u - Sigma_ub %*% t(A)) %*% invQ %*% t(Sigma_u - Sigma_ub %*% t(A))

  out = list(
    bottom_reconciled_mean = mu_b_tilde,
    bottom_reconciled_covariance = Sigma_b_tilde,
    upper_reconciled_mean = mu_u_tilde,
    upper_reconciled_covariance = Sigma_u_tilde
  )
  return(out)
}
