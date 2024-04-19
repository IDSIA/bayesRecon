################################################################################

.DISTR_TYPES = c("continuous", "discrete")
.DISCR_DISTR = c("poisson", "nbinom")
.CONT_DISTR = c("gaussian")

################################################################################
# CHECK INPUT

# Function to check values allowed in S.
.check_S <- function(S) {
  if(!identical(sort(unique(as.vector(S))), c(0,1)) ){
    stop("Input error: S must be a matrix containing only 0s and 1s.")
  }
}

# Checks if a matrix is a covariance matrix (i.e. symmetric p.d.)
.check_cov <- function(cov_matrix, Sigma_str) {
  # Check if the matrix is square
  if (!is.matrix(cov_matrix) || nrow(cov_matrix) != ncol(cov_matrix)) {
    stop(paste0(Sigma_str, " is not square"))
  }
  # Check if the matrix is positive semi-definite
  eigen_values <- eigen(cov_matrix, symmetric = TRUE)$values
  if (any(eigen_values <= 0)) {
    stop(paste0(Sigma_str, " is not positive semi-definite"))
  }
  # Check if the matrix is symmetric
  if (!isSymmetric(cov_matrix)) {
    stop("base_forecasts.Sigma not symmetric")
  }
  # Check if the diagonal elements are non-negative
  if (any(diag(cov_matrix) <= 0)) {
    stop(paste0(Sigma_str, ": some elements on the diagonal are non-positive"))
  }
  # If all checks pass, return TRUE
  return(TRUE)
}

# Checks if the input is a real number
.check_real_number <- function(x) {
  return(length(x)==1 & is.numeric(x))
}

# Checks if the input is a positive number
.check_positive_number <- function(x) {
  return(length(x)==1 && is.numeric(x) && x > 0)
}

# Check that the distr is implemented
.check_implemented_distr <- function(distr) {
  if (!(distr %in% c(.DISCR_DISTR, .CONT_DISTR))) {
    stop(paste(
      "Input error: the distribution must be one of {",
      paste(c(.DISCR_DISTR, .CONT_DISTR), collapse = ', '), "}"))
  }
}

# Check the parameters of distr
.check_distr_params <- function(distr, params) {
  .check_implemented_distr(distr)
  switch(
    distr,
    "gaussian" = {
      mean = params$mean
      sd = params$sd
      if (!.check_real_number(mean)) {
        stop("Input error: mean of Gaussian must be a real number")
      }
      if (!.check_positive_number(sd)) {
        stop("Input error: sd of Gaussian must be a positive number")
      }
    },
    "poisson"  = {
      lambda = params$lambda
      if (!.check_positive_number(lambda)) {
        stop("Input error: lambda of Poisson must be a positive number")
      }
    },
    "nbinom"   = {
      size = params$size,
      prob = params$prob,
      mu = params$mu
      # Check that size is specified, and that is a positive number
      if (is.null(size)) {
        stop("Input error: size parameter for the nbinom distribution must be specified")
      } 
      if (!.check_positive_number(size)) {
        stop("Input error: size of nbinom must be a positive number")
      }
      # Check that exactly one of prob, mu is specified
      n_prob_mu = !is.null(prob) + !is.null(mu)
      if (n_prob_mu == 2) {
        stop("Input error: prob and mu for the nbinom distribution are both specified ")
      } else if (n_prob_mu == 0) {
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
  if (!all.equal(x, as.integer(x))) {
    stop("Input error: samples are not all discrete")
  }
}

# Check input for BUIS (and for MH)
# base_forecasts, in_type, and distr must be list
.check_input_BUIS <- function(S, base_forecasts, in_type, distr) {
  
  .check_S(S)
  
  # Check in_type
  if (!is.list(in_type)) {
    stop("Input error: in_type must be a list")
  }
  if (!(nrow(S) == length(in_type))) {
    stop("Input error: nrow(S) != length(in_type)")
  }
  for(i in 1:nrow(S)){
    if (!(in_type[[i]] %in% c("params", "samples"))) {
      stop("Input error: in_type[[",i,"]] must be either 'samples' or 'params'")
    }
  }
  
  # Check distr and base forecasts
  if (!is.list(distr)) {
    stop("Input error: distr must be a list")
  }
  if (!(nrow(S) == length(distr))) {
    stop("Input error: nrow(S) != length(distr)")
  }
  if (!is.list(base_forecasts)) {
    stop("Input error: base_forecasts must be a list")
  }
  if (!(nrow(S) == length(base_forecasts))) {
    stop("Input error: nrow(S) != length(base_forecasts)")
  }
  for(i in 1:nrow(S)){
    if (in_type[[i]] == "params") {
      .check_distr_params(distr[[i]], base_forecasts[[i]])
    } else if (in_type[[i]] == "samples") {
      if (!(distr[[i]] %in% .DISTR_TYPES)) {
        stop(paste(
          "Input error: the distribution must be one of {",
          paste(.DISTR_TYPES, collapse = ', '), "}"))
      }
      if (distr[[i]] == "discrete")
      # TODO: check sample size?
    } else {
      stop("Input error: in_type[[",i,"]] must be either 'samples' or 'params'")
    }
  }
}

# Check input for TDcond
.check_input_TD <- function(S, fc_bottom, fc_upper, 
                           bottom_in_type, distr,
                           return_pmf, return_samples) {
  
  .check_S(S)
  
  n_b = ncol(S)        # number of bottom TS
  n_u = nrow(S) - n_b  # number of upper TS
  
  if (length(fc_bottom) != n_b) {
    stop("Input error: length of fc_bottom does not match with S")
  }
  if (!(bottom_in_type %in% c("pmf", "samples", "params"))) {
    stop("Input error: bottom_in_type must be either 'pmf', 'samples', or 'params'")
  }
  if (!(return_pmf | return_samples)) {
    stop("Input error: at least one of 'return_pmf' and 'return_samples' must be TRUE")
  } 
  # Check the dimensions of mu and Sigma
  if (length(fc_upper$mu) != n_u | any(dim(fc_upper$Sigma) != c(n_u, n_u))) {
    stop("Input error: the dimensions of the upper parameters do not match with S")
  }
  # Check that Sigma is a covariance matrix (symmetric positive semi-definite)
  .check_cov(fc_upper$Sigma, "Upper covariance matrix")
  
  # If bottom_in_type is not "params" but distr is specified, throw a warning
  if (bottom_in_type %in% c("pmf", "samples") & !is.null(distr)) {
    warning(paste0("Since bottom_in_type = '", bottom_in_type, "', the input distr is ignored"))
  }
  # If bottom_in_type is params, distr must be one of the implemented discrete distr.
  # Also, check the parameters
  if (bottom_in_type == "params") {
    if (!(distr %in% .DISCR_DISTR)) {
      stop(paste0("Input error: distr must be one of {",
                  paste(DISCR_DISTR, collapse = ', '), "}"))
    }
    for (i in 1:n_b) {
      .check_distr_params(distr, fc_bottom[[i]])
    }
  }
}

################################################################################
# SAMPLE

# Sample from one of the implemented distributions
.distr_sample <- function(params, distr, n) {
  .check_distr_params(distr, params)
  switch(
    distr,
    "gaussian" = {
      mean = params$mean
      sd   = params$sd
      samples = stats::rnorm(n = n, mean = mean, sd = sd) },
    "poisson"  = {
      lambda = params$lambda
      samples = stats::rpois(n = n, lambda = lambda) },
    "nbinom"   = {
      size = params$size
      prob = params$prob
      mu   = params$mu
      if (!is.null(prob)) {
        samples = rnbinom(n = n, size = size, prob = prob)
      } else if (!is.null(mu)) {
        samples = rnbinom(n = n, size = size, mu = mu)
      } 
      },
  )
  return(samples)
}

# Sample from a multivariate Gaussian distribution with specified mean and cov. matrix
.MVN_sample <- function(n_samples, mu, Sigma) {
  n = length(mu)
  if (any(dim(Sigma) != c(n,n))) {
    stop("Dimension of mu and Sigma are not compatible!")
  } 
  .check_cov <- function(Sigma, "Sigma")
  
  Z = matrix(rnorm(n*n_samples), ncol = n)
  Ch = chol(Sigma)
  samples = Z %*% Ch + matrix(mu, nrow = n_samples, ncol = n, byrow = TRUE)
  return(samples)
}

################################################################################
# Miscellaneous

# Compute the pmf of the distribution specified by distr and params at the points x
.distr_pmf <- function(x, params, distr) {
  .check_distr_params(distr, params)
  switch(
    distr,
    "gaussian" = {
      mean = params$mean
      sd   = params$sd
      pmf = stats::dnorm(x = x, mean = mean, sd = sd) },
    "poisson"  = {
      lambda = params$lambda
      pmf = stats::dpois(x = x, lambda = lambda) },
    "nbinom"   = {
      size = params$size
      prob = params$prob
      mu   = params$mu
      if (!is.null(prob)) {
        pmf = dnbinom(x = x, size = size, prob = prob)
      } else if (!is.null(mu)) {
        pmf = dnbinom(x = x, size = size, mu = mu)
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

.gen_gaussian <- function(params_file, seed=NULL) {
  set.seed(seed)
  params = utils::read.csv(file = params_file, header = FALSE)
  out = list()
  for (i in 1:nrow(params)) {
    out[[i]] = stats::rnorm(n=1e6, mean=params[[1]][i], sd=params[[2]][i])
  }
  return(out)
}

.gen_poisson <- function(params_file, seed=NULL) {
  set.seed(seed)
  params = utils::read.csv(file = params_file, header = FALSE)
  out = list()
  for (i in 1:nrow(params)) {
    out[[i]] = stats::rpois(n=1e6, lambda=params[[1]][i])
  }
  return(out)
}
