# Input checks
.DISTR_SET = c("gaussian", "poisson", "nbinom")
.DISTR_SET2 = c("continuous", "discrete")
.check_input <- function(S, base_forecasts, in_type, distr) {
  
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
  # TODO if distr is a list, check that entries are coherent
}


# Checks if a matrix is a covariance matrix (i.e. symmetric p.d.)
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


# Function to check values allowed in S.
.check_S <- function(S) {
  if(!identical(sort(unique(as.vector(S))), c(0,1)) ){
    stop("Input error: S must be a matrix containing only 0s and 1s.")
  }
}

# Individual check on the parameter distr
.check_distr <- function(in_type, distr, i=NULL) {
  
  add_string = ""
  if(!is.null(i)){
    add_string = paste("[[",i,"]]")
  }
  
  if (in_type == "params" & !(distr %in% .DISTR_SET)) {
    stop(paste(
      "Input error: if in_type='params', distr", add_string, " must be {",
      paste(.DISTR_SET, collapse = ', '),
      "}"
    ))
  }
  if (in_type == "samples" & !(distr %in% .DISTR_SET2)) {
    stop(paste(
      "Input error: if in_type='samples', distr", add_string, " must be {",
      paste(.DISTR_SET2, collapse = ', '),
      "}"
    ))
  }
}

# Returns TRUE if A is a hierarchy matrix 
# (according to Definition 1 in "Find Maximal Hierarchy")
# If this returns TRUE we avoid solving the integer linear programming problem
.check_hierarchical <- function(A) {
  
  k <- nrow(A)
  m <- ncol(A)
  
  for (i in 1:k) {
    for (j in 1:k) {
      if (i < j) {
        cond1 = A[i,] %*% A[j,] != 0     # Upper i and j have some common descendants
        cond2 = sum(A[i,] - A[j,] >= 0) < m  # Upper j is not a descendant of upper i
        cond3 = sum(A[i,] - A[j,] <= 0) < m  # Upper i is not a descendant of upper j
        if (cond1 & cond2 & cond3) {
          return(FALSE)
        }
      }
    }
  }
  
  return(TRUE)
  
}

# Checks that there is no bottom continuous variable child of a 
# discrete upper variable
.check_hierfamily_rel <- function(sh.res, distr, debug=FALSE) {
  for (bi in seq_along(distr[sh.res$bottom_idxs])) {
    distr_bottom = distr[sh.res$bottom_idxs][[bi]]
    rel_upper_i = sh.res$A[,bi]
    rel_distr_upper = unlist(distr[sh.res$upper_idxs])[rel_upper_i == 1]
    err_message = "A continuous bottom distribution is child of a discrete one."
    if (distr_bottom == .DISTR_SET2[1]) {              
      if (sum(rel_distr_upper == .DISTR_SET2[2]) | 
          sum(rel_distr_upper == .DISTR_SET[2]) | sum(rel_distr_upper == .DISTR_SET[3])) {
        if (debug) { return(-1) } else { stop(err_message) }
      }
    }
    if (distr_bottom == .DISTR_SET[1]) {                
      if (sum(rel_distr_upper == .DISTR_SET2[2]) |
          sum(rel_distr_upper == .DISTR_SET[2]) | sum(rel_distr_upper == .DISTR_SET[3])) {
        if (debug) { return(-1) } else { stop(err_message) }
      }
    }
  }
  if (debug) { return(0) }
}


# Misc
.shape <- function(m) {
  print(paste0("(", nrow(m), ",", ncol(m), ")"))
}

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
