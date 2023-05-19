#-------------------------------------------------------------------------------
# # Minimal example
# S <- rbind(matrix(1, 1, 2), diag(1,2))
# distr <- rep(list("poisson"), 3)
# params <- list(c(8), c(3), c(4))
# out <- reconc_MCMC(S, distr, params)
#-------------------------------------------------------------------------------
#' @title Probabilistic reconciliation with MCMC
#' @param S aggregating matrix
#' @param distr list of strings specifying the distribution of each variable
#' @param params list of the parameters of the distributions
#' @param N_samples number of samples to draw using MCMC
#' @param tuning_int number of iterations between scale updates of the proposal
#' @param init_scale initial scale of the proposal
#' @param burn_in number of initial samples to be discarded
#' @param seed Seed for randomness reproducibility
#' @return samples from the reconciled distribution
#' @export
reconc_MCMC <- function(S,
                        distr,
                        params,
                        N_samples = 10000,
                        tuning_int = 100,
                        init_scale = 1,
                        burn_in = 1000,
                        seed = NULL) {

  set.seed(seed)

  # Check input
  # TODO: use function to check input
  # .check_input of the package?
  # ok only if Poisson or NegBin?

  m <- ncol(S)
  n <- nrow(S)

  # Set the covariance matrix of the proposal (identity matrix)
  cov_mat_prop <- diag(rep(1,m))

  #Set the counter for tuning
  c_tuning <- tuning_int

  # Initialize acceptance counter
  accept_count <- 0

  # Create empty matrix for the samples from MCMC
  b <- matrix(nrow = N_samples, ncol = m)

  # Initialize first sample (draw from base distribution)
  #TODO: use function to get distr_b, params_b
  distr_b <- distr[(n-m+1):n]
  params_b <- params[(n-m+1):n]
  b[1,] <- .initialize_b(distr_b, params_b)

  # Initialize prop list
  old_prop <- list(
    "b" = b[1,],
    "scale" = init_scale
  )

  # Run the chain
  for (i in 2:N_samples) {

    if (c_tuning == 0) {
      old_prop$acc_rate <- accept_count / tuning_int  #set acc_rate
      accept_count <- 0                               #reset acceptance counter
      c_tuning <- tuning_int                          #reset tuning counter
    }

    prop <- .proposal(old_prop, cov_mat_prop)
    b_prop <- prop$b
    alpha <- .accept_prob(b_prop, b[i-1,], S, distr, params)

    if (stats::runif(1) < alpha) {
      b[i,] <- b_prop
      accept_count <- accept_count + 1
    } else {
      b[i,] <- b[i-1,]
    }

    old_prop <- list(
      "b" = b[i,],
      "scale" = prop$scale
    )

    c_tuning <- c_tuning - 1

  }

  burn_in <- min(burn_in, N_samples/2)  # keep at least half of the samples

  b_samples <- b[(burn_in+1) : N_samples, ]
  # u<-samples <- b_samples %*% t(A)
  y_samples <- b_samples %*% t(S)

  out = list(
    bottom_reconciled_samples = b_samples,
    # upper_reconciled_samples = u_samples,
    reconciled_samples = y_samples
  )

  return(out)

}


#-------------------------------------------------------------------------------
##################################
.initialize_b <- function(distr_b, params_b) {

  m <- length(distr_b)

  b <- c()
  for (i in 1:m) {
    b[i] <- .distr_sample(params_b[[i]], distr_b[[i]], 1)
  }

  return(b)

}



#-------------------------------------------------------------------------------
##################################
#' @title Compute acceptance probability
#' @param b proposal state
#' @param b0 current state
#' @param S aggregating matrix
#' @param distr list of strings specifying the distribution of each variable
#' @param params list of the parameters of the distributions
#' @return the acceptance probability alpha
.accept_prob <- function(b, b0, S, distr, params) {

  alpha <- .target_pmf(b, S, distr, params) / .target_pmf(b0, S, distr, params)

  return(min(1,alpha))

}


##################################
.target_pmf <- function(b, S, distr, params) {

  m <- ncol(S)
  n <- nrow(S)
  k <- n-m

  if (length(distr)!=n | length(params)!= n){
    print("Dimensions are wrong")
  }
  if (length(b)!=m) {
    print("Matrix S does not match with b")
  }

  y <- S %*% b

  pmf <- 1
  for (j in 1:n) {
    pmf <- pmf * .distr_pmf(y[[j]], distr[[j]], params[[j]])
  }

  return(pmf)

}


##################################
.distr_pmf <- function(x, distr_, par) {
  switch(
    distr_,
    "poisson"  = {
      pmf = stats::dpois(x=x, lambda = par[[1]]) },
    "negbin"   = {
      pmf = stats::dnbinom(x=x, size = par[[1]], prob = par[[2]]) },
  )
  return(pmf)
}


#-------------------------------------------------------------------------------
##################################
#' @title Generate a proposal for MCMC step
#' @param prev_prop a list containing b, scale, acc_rate
#' @param cov_mat_prop the covariance matrix (ncol(cov_mat_prop)=length(b)) for the normal proposal
#' @return a list containing the new b, scale, acc_rate
.proposal <- function(prev_prop, cov_mat_prop){

  # extract previous proposal
  b0 <- prev_prop$b
  old_scale <- prev_prop$scale
  acc_rate <- prev_prop$acc_rate

  if(!is.null(acc_rate)){
    # Compute the scaling parameter
    scale <- .tune(old_scale,acc_rate)
  }else{
    scale <- old_scale
  }


  n_x <- length(b0)

  if(n_x != ncol(cov_mat_prop)){
    stop(sprintf("Error in .proposal: previous state dim (%s) and
                 covariance dim (%s) do not match", n_x, ncol(cov_mat_prop)))
  }

  # We always do an independent Gaussian proposal
  dd <- stats::rnorm(n_x)*diag(cov_mat_prop)*scale

  # CHANGE HERE for mixed case
  dd <- round(dd,0)
  b = b0+dd

  prop <- list(b=b, acc_rate=acc_rate, scale=scale)
  return(prop)
}


# scaling parameter
# we use the same variance adaptation used in pymc3
# Rate    Variance adaptation
# ----    -------------------
# <0.001        x 0.1
# <0.05         x 0.5
# <0.2          x 0.9
# >0.5          x 1.1
# >0.75         x 2
# >0.95         x 10
.tune <- function(scale, acc_rate){
  if(acc_rate<0.001){
    return(scale*0.1)
  }else if(acc_rate<0.05){
    return(scale*0.5)
  }else if(acc_rate<0.2){
    return(scale*0.9)
  }else if(acc_rate>0.5){
    return(scale*1.1)
  }else if(acc_rate>0.75){
    return(scale*2)
  }else if(acc_rate>0.95){
    return(scale*10)
  }

  return(scale)

}

