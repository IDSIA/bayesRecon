################################################################################
# IMPLEMENTED DISTRIBUTIONS

.DISTR_TYPES <- c("continuous", "discrete")
.DISCR_DISTR <- c("poisson", "nbinom")
.CONT_DISTR <- c("gaussian")

################################################################################
# PARAMETERS FOR PMF CONVOLUTION AND SMOOTHING

.TOLL <- 1e-15
.RTOLL <- 1e-9
.ALPHA_SMOOTHING <- 1e-9
.LAP_SMOOTHING <- FALSE

################################################################################
# OTHER PARAMETERS

.NEGBIN_TOLL <- 1e-6 # used when fitting a Negative Binomial distribution
.L_SHRINK_RECONC_T <- 1e-4  # used for shrinking the empirical covariance matrix in reconc_t

################################################################################
# CHECK INPUT

# Function to check values allowed in S.
.check_S <- function(S) {
  if (!identical(sort(unique(as.vector(S))), c(0, 1))) {
    stop("Input error in S: S must be a matrix containing only 0s and 1s.")
  }

  if (!all(colSums(S) > 1)) {
    stop("Input error in S: all bottom level forecasts must aggregate into an upper.")
  }

  if (nrow(unique(S)) != nrow(S)) {
    warning("S has some repeated rows.")
  }

  # Check that each bottom has a corresponding row with with one 1 and the rest 0s.
  if (nrow(unique(S[rowSums(S) == 1, ])) < ncol(S)) {
    stop("Input error in S: there is at least one bottom that does not have a row with one 1 and the rest 0s.")
  }
}

# Function to check aggregation matrix A
.check_A <- function(A) {
  if (!all(A %in% c(0, 1))) {
    stop("Input error in A: A must be a matrix containing only 0s and 1s.")
  }

  if (any(colSums(A) == 0)) {
    stop("Input error in A: some columns do not have any 1.
          All bottom level forecasts must aggregate into an upper.")
  }

  if (nrow(unique(A)) != nrow(A)) {
    warning("A has some repeated rows.")
  }
}

# Check if it is a covariance matrix (i.e. symmetric p.d.)
.check_cov <- function(cov_matrix, Sigma_str, pd_check = FALSE, symm_check = FALSE) {
  # Check if the matrix is square
  if (!is.matrix(cov_matrix) || nrow(cov_matrix) != ncol(cov_matrix)) {
    stop(paste0(Sigma_str, " is not square"))
  }

  # Check if the matrix is positive semi-definite
  if (pd_check) {
    eigen_values <- eigen(cov_matrix, symmetric = TRUE)$values
    if (any(eigen_values <= 0)) {
      stop(paste0(Sigma_str, " is not positive semi-definite"))
    }
  }
  if (symm_check) {
    # Check if the matrix is symmetric
    if (!isSymmetric(cov_matrix)) {
      stop(paste0(Sigma_str, " is not symmetric"))
    }
  }
  # Check if the diagonal elements are non-negative
  if (any(diag(cov_matrix) < 0)) {
    stop(paste0(Sigma_str, ": some elements on the diagonal are negative"))
  }
  # If all checks pass, return TRUE
  return(TRUE)
}

# Checks if the input is a real number
.check_real_number <- function(x) {
  return(length(x) == 1 & is.numeric(x))
}

# Checks if the input is a positive number
.check_positive_number <- function(x) {
  return(length(x) == 1 && is.numeric(x) && x > 0)
}

# Check that the distr is implemented
.check_implemented_distr <- function(distr) {
  if (!(distr %in% c(.DISCR_DISTR, .CONT_DISTR))) {
    stop(paste(
      "Input error: the distribution must be one of {",
      paste(c(.DISCR_DISTR, .CONT_DISTR), collapse = ", "), "}"
    ))
  }
}

# Check the parameters of distr
.check_distr_params <- function(distr, params) {
  .check_implemented_distr(distr)
  if (!is.list(params)) {
    stop("Input error: the parameters of the distribution must be given as a list.")
  }
  switch(distr,
    "gaussian" = {
      mean <- params$mean
      sd <- params$sd
      if (!.check_real_number(mean)) {
        stop("Input error: mean of Gaussian must be a real number")
      }
      if (!.check_positive_number(sd)) {
        stop("Input error: sd of Gaussian must be a positive number")
      }
    },
    "poisson" = {
      lambda <- params$lambda
      if (!.check_positive_number(lambda)) {
        stop("Input error: lambda of Poisson must be a positive number")
      }
    },
    "nbinom" = {
      size <- params$size
      prob <- params$prob
      mu <- params$mu
      # Check that size is specified, and that is a positive number
      if (is.null(size)) {
        stop("Input error: size parameter for the nbinom distribution must be specified")
      }
      if (!.check_positive_number(size)) {
        stop("Input error: size of nbinom must be a positive number")
      }
      # Check that exactly one of prob, mu is specified
      if (!is.null(prob) & !is.null(mu)) {
        stop("Input error: prob and mu for the nbinom distribution are both specified ")
      } else if (is.null(prob) & is.null(mu)) {
        stop("Input error: either prob or mu must be specified")
      } else {
        if (!is.null(prob)) {
          if (!.check_positive_number(prob) | prob > 1) {
            stop("Input error: prob of nbinom must be positive and <= 1")
          }
        } else if (!is.null(mu)) {
          if (!.check_positive_number(mu)) {
            stop("Input error: mu of nbinom must be positive")
          }
        }
      }
    },
  )
}

# Check that the samples are discrete
.check_discrete_samples <- function(samples) {
  if (!isTRUE(as.vector(all.equal(as.vector(samples), as.integer(samples))))) {
    stop("Input error: samples are not all discrete")
  }
}

# Check input for BUIS (and for MH)
# base_fc, in_type, and distr must be list
.check_input_BUIS <- function(A, base_fc, in_type, distr) {
  .check_A(A)

  n_tot_A <- ncol(A) + nrow(A)

  # Check in_type
  if (!is.list(in_type)) {
    stop("Input error: in_type must be a list")
  }
  if (!(n_tot_A == length(in_type))) {
    stop("Input error: ncol(A)+nrow(A) != length(in_type)")
  }
  for (i in 1:n_tot_A) {
    if (!(in_type[[i]] %in% c("params", "samples"))) {
      stop("Input error: in_type[[", i, "]] must be either 'samples' or 'params'")
    }
  }

  # Check distr and base forecasts
  if (!is.list(distr)) {
    stop("Input error: distr must be a list")
  }
  if (!(n_tot_A == length(distr))) {
    stop("Input error: ncol(A)+nrow(A) != length(distr)")
  }
  if (!is.list(base_fc)) {
    stop("Input error: base_fc must be a list")
  }
  if (!(n_tot_A == length(base_fc))) {
    stop("Input error: ncol(A)+nrow(A) != length(base_fc)")
  }
  for (i in 1:n_tot_A) {
    if (in_type[[i]] == "params") {
      .check_distr_params(distr[[i]], base_fc[[i]])
    } else if (in_type[[i]] == "samples") {
      if (!(distr[[i]] %in% .DISTR_TYPES)) {
        stop(paste(
          "Input error: the distribution must be one of {",
          paste(.DISTR_TYPES, collapse = ", "), "}"
        ))
      }
      if (distr[[i]] == "discrete") {
        .check_discrete_samples(base_fc[[i]])
      }
      # TODO: check sample size?
    } else {
      stop("Input error: in_type[[", i, "]] must be either 'samples' or 'params'")
    }
  }
}

# Check input for TDcond
.check_input_TD <- function(A, base_fc_bottom, base_fc_upper,
                            bottom_in_type, distr,
                            return_type) {
  .check_A(A)

  n_b <- ncol(A) # number of bottom TS
  n_u <- nrow(A) # number of upper TS

  if (!(bottom_in_type %in% c("pmf", "samples", "params"))) {
    stop("Input error: bottom_in_type must be either 'pmf', 'samples', or 'params'")
  }
  if (!(return_type %in% c("pmf", "samples", "all"))) {
    stop("Input error: return_type must be either 'pmf', 'samples', or 'all'")
  }
  if (length(base_fc_bottom) != n_b) {
    stop("Input error: length of base_fc_bottom does not match with A")
  }
  # If cov is a number, transform into a matrix
  if (length(base_fc_upper$cov) == 1) {
    base_fc_upper$cov <- as.matrix(base_fc_upper$cov)
  }
  # Check the dimensions of mean and cov
  if (length(base_fc_upper$mean) != n_u | any(dim(base_fc_upper$cov) != c(n_u, n_u))) {
    stop("Input error: the dimensions of the upper parameters do not match with A")
  }
  # Check that cov is a covariance matrix (symmetric positive semi-definite)
  .check_cov(base_fc_upper$cov, "Upper covariance matrix", symm_check = TRUE)

  # If bottom_in_type is not "params" but distr is specified, throw a warning
  if (bottom_in_type %in% c("pmf", "samples") & !is.null(distr)) {
    warning(paste0("Since bottom_in_type = '", bottom_in_type, "', the input distr is ignored"))
  }
  # If bottom_in_type is params, distr must be one of the implemented discrete distr.
  # Also, check the parameters
  if (bottom_in_type == "params") {
    if (is.null(distr)) {
      stop("Input error: if bottom_in_type = 'params', distr must be specified")
    }
    if (!(distr %in% .DISCR_DISTR)) {
      stop(paste0(
        "Input error: distr must be one of {",
        paste(.DISCR_DISTR, collapse = ", "), "}"
      ))
    }
    for (i in 1:n_b) {
      .check_distr_params(distr, base_fc_bottom[[i]])
    }
  }
}

.check_input_t <- function(A, base_fc_mean, y_train, residuals, ...) {
  .check_A(A)

  n_b <- ncol(A) # number of bottom TS
  n_u <- nrow(A) # number of upper TS

  if (!is.vector(base_fc_mean)) {
    stop("Input error: base_fc_mean must be a vector")
  }
  n <- length(base_fc_mean)
  if (n_u + n_b != n) {
    stop("Input error: the length of base_fc_mean must be equal to nrow(A) + ncol(A)")
  }
  
  add_args <- list(...)
  prior <- add_args$prior
  posterior <- add_args$posterior
  freq <- add_args$freq
  criterion <- add_args$criterion

  ##############################################################################
  ### CASE 1 ###
  # If posterior is provided, check if is a list with entries nu and Psi and extract values
  if (!is.null(posterior)) {
    if (is.list(posterior)) {
      nu_post <- posterior$nu
      Psi_post <- posterior$Psi
      if (is.null(nu_post) | is.null(Psi_post)) {
        stop("Input error: posterior must be a list with entries nu and Psi")
      } else if (!is.numeric(nu_post) | length(nu_post) != 1 | nu_post <= n-1) {
        stop("Input error: nu in posterior must be a number greater than n. of series - 1")
      } else if (!is.matrix(Psi_post) | any(dim(Psi_post) != c(n, n))) {
        stop("Input error: Psi in posterior must be a matrix with dimensions compatible with base_fc_mean")
      }
      # If posterior is provided, then y_train, residuals, and prior are ignored:
      # if they are provided, throw a warning
      if (!is.null(y_train)) {
        warning("Input warning: posterior is provided, ignoring y_train")
      }
      if (!is.null(residuals)) {
        warning("Input warning: posterior is provided, ignoring residuals")
      }
      if (!is.null(prior)) {
        warning("Input warning: posterior is provided, ignoring prior")
      }
    } else {
      stop("Input error: posterior must be a list with entries nu and Psi")
    }

    ### CASE 2 ###
    # If posterior not provided, first check that residuals are provided
  } else {
    if (is.null(residuals)) {
      stop("Input error: either posterior or residuals must be provided")
    }
    if (!is.matrix(residuals)) {
      stop("Input error: residuals must be a matrix")
    }
    if (ncol(residuals) != n) {
      stop("Input error: number of columns of residuals must be equal to length of base_fc_mean")
    }

    L <- nrow(residuals) # number of residual samples (i.e., training length)
    if (L < 10) {
      warning("Warning: number of rows of residuals is less than 10, covariance estimation may be inaccurate")
    }
    # TODO: implement fallback

    ### CASE 2a ###
    # If prior is provided, check if is a list with entries nu and Psi and extract values
    if (!is.null(prior)) {
      if (is.list(prior)) {
        nu_prior <- prior$nu
        Psi_prior <- prior$Psi
        if (is.null(nu_prior) | is.null(Psi_prior)) {
          stop("Input error: prior must be a list with entries nu and Psi")
        } else if (!is.numeric(nu_prior) | length(nu_prior) != 1 | nu_prior <= n-1) {
          stop("Input error: nu in prior must be a number greater than n. of series - 1")
        } else if (!is.matrix(Psi_prior) | any(dim(Psi_prior) != c(n, n))) {
          stop("Input error: Psi in prior must be a matrix with dimensions compatible with base_fc_mean")
        }
        # If prior is provided, then y_train is ignored: if it is provided, throw a warning
        if (!is.null(y_train)) {
          warning("Input warning: prior is provided, ignoring y_train")
        }
      } else {
        stop("Input error: prior must be a list with entries nu and Psi")
      }
    }

    ### CASE 2b ###
    # If prior not provided:
    # - compute Psi using the (shrinked) covariance matrix of the residuals of the naive
    #   or seasonal naive forecasts
    # - set nu using LOOCV
    else {
      if (is.null(y_train)) {
        stop("Input error: y_train must be provided when neither prior nor posterior are given")
      }
      if (!is.matrix(y_train)) {  
        stop("Input error: y_train must be a matrix")
      }
      # this works also if y_train is a mts; TODO: allow for other types, e.g. data.frame
      if (ncol(y_train) != n) {
        stop("Input error: number of columns of y_train must be equal to length of base_fc_mean")
      }
      if (nrow(y_train) != L) {
        warning("Numbers of rows of y_train and of residuals are different!")
      }
      
      if (!is.null(freq)) {
        if (!is.numeric(freq) | length(freq) != 1 | freq < 1 | (freq %% 1) != 0) {
          stop("Input error: freq must be a positive integer")
        }
      }
      
      # Current logic: 
      # * if y_train is a mts object: check if freq is provided and is equal to the frequency of y_train; 
      #                               if they are different, throw a warning and ignore the frequency of y_train (see compute_naive_cov function)
      # * if y_train is not a mts object: check if freq is provided and is a positive integer; 
      #                                   if not, throw a warning  
      if (stats::is.mts(y_train)) {
        if (!is.null(freq) && freq != stats::frequency(y_train)) {
          warning("Input warning: the provided freq is different from the frequency of y_train.
                   The frequency of y_train will be ignored.")
        } 
      } else {
        if (is.null(freq)) {
          warning("Input warning: y_train is not a mts object and freq is not provided. 
                   The prior will be set assuming no seasonality.")
        }
      }
      
      if (!is.null(criterion)) {
        if (!(criterion %in% c("RSS", "seas-test"))) {
          stop("Input error: criterion must be either 'RSS' or 'seas-test'")
        }
      }
    }
  }
}

# Check importance sampling weights
.check_weights <- function(w, n_eff_min = 200, p_n_eff = 0.01) {
  warning <- FALSE
  warning_code <- c()
  warning_msg <- c()

  n <- length(w)
  n_eff <- n

  # 1. w==0
  if (all(w == 0)) {
    warning <- TRUE
    warning_code <- c(warning_code, 1)
    warning_msg <- c(
      warning_msg,
      "Importance Sampling: all the weights are zeros. This is probably caused by a strong incoherence between bottom and upper base forecasts."
    )
  } else {
    # Effective sample size
    w <- w / sum(w)
    n_eff <- 1 / sum(w^2)

    # 2. n_eff < threshold
    if (n_eff < n_eff_min) {
      warning <- TRUE
      warning_code <- c(warning_code, 2)
      warning_msg <- c(
        warning_msg,
        paste0("Importance Sampling: effective_sample_size= ", round(n_eff, 2), " (< ", n_eff_min, ").")
      )
    }

    # 3. n_eff < p*n, e.g. p = 0.05
    if (n_eff < p_n_eff * n) {
      warning <- TRUE
      warning_code <- c(warning_code, 3)
      warning_msg <- c(
        warning_msg,
        paste0("Importance Sampling: effective_sample_size= ", round(n_eff, 2), " (< ", round(p_n_eff * 100, 2), "%).")
      )
    }
  }
  res <- list(
    warning = warning,
    warning_code = warning_code,
    warning_msg = warning_msg,
    n_eff = n_eff
  )

  return(res)
}

################################################################################
# SAMPLE

# Sample from one of the implemented distributions
.distr_sample <- function(params, distr, n) {
  .check_distr_params(distr, params)
  switch(distr,
    "gaussian" = {
      mean <- params$mean
      sd <- params$sd
      samples <- stats::rnorm(n = n, mean = mean, sd = sd)
    },
    "poisson" = {
      lambda <- params$lambda
      samples <- stats::rpois(n = n, lambda = lambda)
    },
    "nbinom" = {
      size <- params$size
      prob <- params$prob
      mu <- params$mu
      if (!is.null(prob)) {
        samples <- stats::rnbinom(n = n, size = size, prob = prob)
      } else if (!is.null(mu)) {
        samples <- stats::rnbinom(n = n, size = size, mu = mu)
      }
    },
  )
  return(samples)
}

# Sample from a multivariate Gaussian distribution with specified mean and cov. matrix
.MVN_sample <- function(n_samples, mu, Sigma) {
  n <- length(mu)
  if (any(dim(Sigma) != c(n, n))) {
    stop("Dimensions of mean and covariance matrix are not compatible!")
  }
  .check_cov(Sigma, "Covariance matrix", pd_check = FALSE, symm_check = FALSE)

  Z <- matrix(stats::rnorm(n * n_samples), ncol = n)

  Ch <- tryCatch(base::chol(Sigma),
    error = function(e) stop(paste0(e, "check the covariance in .MVN_sample, the Cholesky fails."))
  )

  samples <- Z %*% Ch + matrix(mu, nrow = n_samples, ncol = n, byrow = TRUE)
  return(samples)
}

# Compute the MVN density
.MVN_density <- function(x, mu, Sigma, max_size_x = 5e3, suppress_warnings = TRUE) {
  # save dimension of mu
  n <- length(mu)

  # Check Sigma
  if (any(dim(Sigma) != c(n, n))) {
    stop("Dimension of mu and Sigma are not compatible!")
  }
  .check_cov(Sigma, "Sigma", pd_check = FALSE, symm_check = FALSE)

  # x must be a matrix with ncol = n (nrow is the number of points to evaluate)
  # or a vector with length n (in which case it is transformed into a matrix)
  if (is.vector(x)) {
    if (length(x) != n) stop("Length of x must be the same of mu")
    x <- matrix(x, ncol = length(x))
  } else if (is.matrix(x)) {
    if (ncol(x) != n) stop("The number of columns of x must be equal to the length of mu")
  } else {
    stop("x must be either a vector or a matrix")
  }

  # Compute Cholesky of Sigma
  chol_S <- tryCatch(base::chol(Sigma),
    error = function(e) stop(paste0(e, "check the covariance in .MVN_density, the Cholesky fails."))
  )

  # Constant of the loglikelihood (computed here because it is always the same)
  const <- -sum(log(diag(chol_S))) - 0.5 * n * log(2 * pi)

  # This part breaks down the evaluation of the density eval into batches, for memory
  rows_x <- nrow(x)

  if (rows_x > max_size_x) {
    logval <- rep(0, rows_x)

    # Compute how many batches we need
    num_backsolves <- rows_x %/% max_size_x

    if (!suppress_warnings) {
      warning_msg <- paste0("x has ", rows_x, " rows, the density evaluation is broken down into ", num_backsolves, " pieces for memory preservation.")
      warning(warning_msg)
    }

    for (j in seq(num_backsolves)) {
      idx_to_select <- (1 + (j - 1) * max_size_x):((j) * max_size_x)
      # Do one backsolve for each batch
      tmp <- backsolve(chol_S, t(x[idx_to_select, ]) - mu, transpose = TRUE)
      rss <- colSums(tmp^2)

      # Update the logval for those indices
      logval[idx_to_select] <- const - 0.5 * rss
    }

    # Last indices: if the number of rows of x is not exactly divided by the size of the batches
    remainder <- rows_x %% max_size_x
    if (remainder != 0) {
      idx_to_select <- (1 + (num_backsolves) * max_size_x):(remainder + (num_backsolves) * max_size_x)
      # Do backsolve on the remaining indices
      tmp <- backsolve(chol_S, t(x[idx_to_select, ]) - mu, transpose = TRUE)
      rss <- colSums(tmp^2)

      logval[idx_to_select] <- const - 0.5 * rss
    }
  } else {
    tmp <- backsolve(chol_S, t(x) - mu, transpose = TRUE)
    rss <- colSums(tmp^2)

    logval <- const - 0.5 * rss
  }

  return(exp(logval))
}

# Resample from weighted sample
.resample <- function(S_, weights, num_samples = NA) {
  if (is.na(num_samples)) {
    num_samples <- length(weights)
  }

  if (nrow(S_) != length(weights)) {
    stop("Error in .resample: nrow(S_) != length(weights)")
  }

  tmp_idx <- sample(x = 1:nrow(S_), num_samples, replace = TRUE, prob = weights)
  return(S_[tmp_idx, ])
}

################################################################################
# Miscellaneous

# Compute the pmf of the distribution specified by distr and params at the points x
.distr_pmf <- function(x, params, distr) {
  .check_distr_params(distr, params)
  switch(distr,
    "gaussian" = {
      mean <- params$mean
      sd <- params$sd
      pmf <- stats::dnorm(x = x, mean = mean, sd = sd)
    },
    "poisson" = {
      lambda <- params$lambda
      pmf <- stats::dpois(x = x, lambda = lambda)
    },
    "nbinom" = {
      size <- params$size
      prob <- params$prob
      mu <- params$mu
      if (!is.null(prob)) {
        pmf <- stats::dnbinom(x = x, size = size, prob = prob)
      } else if (!is.null(mu)) {
        pmf <- stats::dnbinom(x = x, size = size, mu = mu)
      }
    },
  )
  return(pmf)
}

.shape <- function(m) {
  print(paste0("(", nrow(m), ",", ncol(m), ")"))
}

################################################################################
# Functions for tests

.gen_gaussian <- function(params_file, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  params <- utils::read.csv(file = params_file, header = FALSE)
  out <- list()
  for (i in 1:nrow(params)) {
    out[[i]] <- stats::rnorm(n = 1e6, mean = params[[1]][i], sd = params[[2]][i])
  }
  return(out)
}

.gen_poisson <- function(params_file, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  params <- utils::read.csv(file = params_file, header = FALSE)
  out <- list()
  for (i in 1:nrow(params)) {
    out[[i]] <- stats::rpois(n = 1e6, lambda = params[[1]][i])
  }
  return(out)
}
