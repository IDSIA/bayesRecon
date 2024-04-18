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

# Check that distr is one of the allowed values, depending on in_type 
# (see .DISTR_TYPES, .DISCR_DISTR, .CONT_DISTR)
.check_distr <- function(in_type, distr, i=NULL) {
  
  add_string = ""
  if(!is.null(i)){
    add_string = paste("[[",i,"]]")
  }
  
  if (in_type == "params" & !(distr %in% c(.DISCR_DISTR, .CONT_DISTR))) {
    stop(paste(
      "Input error: if in_type='params', distr", add_string, " must be one of {",
      paste(c(.DISCR_DISTR, .CONT_DISTR), collapse = ', '),
      "}"
    ))
  }
  if (in_type == "samples" & !(distr %in% .DISTR_TYPES)) {
    stop(paste(
      "Input error: if in_type='samples', distr", add_string, " must be one of {",
      paste(.DISTR_TYPES, collapse = ', '),
      "}"
    ))
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

# Check the parameters of one of the implemented distribution
.check_distr_params(distr, params) {
  if (!(distr %in% c(.DISCR_DISTR, .CONT_DISTR))) {
    # Check that the distr is implemented
    stop(paste(
      "Input error: the distribution must be one of {",
      paste(c(.DISCR_DISTR, .CONT_DISTR), collapse = ', '), "}"))
  }
  switch(
    distr,
    "gaussian" = {
      mean = params[[1]]
      sd = params[[2]]
      # mean must be a real number
      if (length(mean)!=1 | !is.numeric(mean)) {
        stop("Input error: the mean of a Gaussian must be a real number")
      }
      # std must be a positive number
      if (length(sd)!=1 || !is.numeric(sd) || sd <= 0) {
        stop("Input error: the mean of a Gaussian must be a real number")
      }
    },
    "poisson"  = {
      lambda = params[[1]]
      # lambda must be a positive number
      if (length(lambda)!=1 || !is.numeric(lambda) || lambda <= 0) {
        stop("Input error: the mean of a Gaussian must be a real number")
      }
    },
    "nbinom"   = {
      # TODO
    },
  )
}

# Check input for BUIS (and for MH)
.check_input_BUIS <- function(S, base_forecasts, in_type, distr) {
  
  .check_S(S)
  
  if (!(nrow(S) == length(base_forecasts))) {
    stop("Input error: nrow(S) != length(base_forecasts)")
  }
  
  if (is.character(in_type)) {
    if (!(in_type %in% c("params", "samples"))) {
      stop("Input error: in_type must be either 'samples' or 'params'")
    }
    flag_list_intype = FALSE
  }else if (is.list(in_type)) {
    if (!(nrow(S) == length(in_type))) {
      stop("Input error: nrow(S) != length(in_type)")
    }
    for(i in 1:nrow(S)){
      if (!(in_type[[i]] %in% c("params", "samples"))) {
        stop("Input error: in_type[[",i,"]] must be either 'samples' or 'params'")
      }
    }
    flag_list_intype = TRUE
  }else{
    stop("Input error: in_type must be either a string or a list, check documentation")
  }
  
  if (is.character(distr) &
      length(distr) == 1) {
    # if distr is a string...
    if (flag_list_intype){
      for(i in 1:nrow(S)){
        .check_distr(in_type[[i]], distr)
      }
    }else{
      .check_distr(in_type, distr)
    }
  }else  if (is.list(distr)) {
    if (!(nrow(S) == length(distr))) {
      stop("Input error: nrow(S) != length(distr)")
    }
    for(i in 1:nrow(S)){
      if (flag_list_intype){
        .check_distr(in_type[[i]], distr[[i]],i)
      }else{
        .check_distr(in_type, distr[[i]],i)
      }
      
    }
  }else{
    stop("Input error: distr must be either a string or a list, check documentation")
  }

  # Eventually:
  # TODO if in_type=='samples' check sample sizes
  # TODO if in_type=='params' check that:
  # - gaussian: 2 parameters, (mu, sd)
  # - poisson:  1 parameter,  (lambda)
  # - nbinom:   2 parameters, (n, p)
}

# Check input for TDcond
.check_input_TD = function(S, fc_bottom, fc_upper, 
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
  # If bottom_in_type is params, distr must be one of the implemented discrete distr
  if (bottom_in_type == "params" & !(distr %in% .DISCR_DISTR)) {
    stop(paste0("Input error: distr must be one of {",
                paste(DISCR_DISTR, collapse = ', '), "}"))
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

# Sample from a multivariate Gaussian distribution with specified mean and cov. matrix
.MVN_sample = function(n_samples, mu, Sigma) {
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
