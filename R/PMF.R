# A pmf is represented as normalized numeric vector v: 
# for each j = 0, ..., M, the probability of j is the value v[[j+1]]

###

# Compute the empirical pmf from a vector of samples
PMF.from_samples = function(v) {
  .check_discrete_samples(v)
  pmf = tabulate(v+1) / length(v)  # the support starts from 0 
  # Tabulate only counts values above 1: if sum(tabulate(v+1)) > length(v),
  # it means that there were negative samples
  if (!isTRUE(all.equal(sum(pmf), 1))) {
    stop("Input error: same samples are negative")
  }
  return(pmf)
}

# Compute the pmf from a parametric distribution
PMF.from_params = function(params, distr, Rtoll = .RTOLL) {
  # Check that the distribution is implemented, and that the params are ok
  if (!(distr %in% .DISCR_DISTR)) {
    stop(paste0("Input error: distr must be one of {",
                paste(.DISCR_DISTR, collapse = ', '), "}"))
  }
  .check_distr_params(distr, params)
  # Compute the pmf
  switch(
    distr,
    "poisson"  = {
      lambda = params$lambda
      M = stats::qpois(1-Rtoll, lambda)
      pmf = stats::dpois(0:M, lambda)},
    "nbinom"   = {
      size = params$size
      prob = params$prob
      mu   = params$mu
      if (!is.null(prob)) {
        M   = stats::qnbinom(1-Rtoll, size=size, prob=prob)
        pmf = stats::dnbinom(0:M, size=size, prob=prob)
      } else if (!is.null(mu)) {
        M   = stats::qnbinom(1-Rtoll, size=size, mu=mu)
        pmf = stats::dnbinom(0:M, size=size, mu=mu)
        }
      },
  )
  pmf = pmf / sum(pmf)
  return(pmf)
}

#' @title Sample from the distribution given as a PMF object
#'
#' @description
#' 
#' Samples (with replacement) from the probability distribution specified by `pmf`.
#' 
#' @param pmf the PMF object.
#' @param N_samples number of samples. 
#' 
#' @return Samples drawn from the distribution specified by `pmf`.  
#' @seealso [PMF.get_mean()], [PMF.get_var()], [PMF.get_quantile()], [PMF.summary()]
#' 
#' @examples 
#' library(bayesRecon)
#' 
#' # Let's build the pmf of a Binomial distribution with parameters n and p
#' n <- 10
#' p <- 0.6 
#' pmf_binomial <- apply(matrix(seq(0,n)),MARGIN=1,FUN=function(x) dbinom(x,size=n,prob=p))
#' 
#' # Draw samples from the PMF object
#' set.seed(1)
#' samples <- PMF.sample(pmf=pmf_binomial,N_samples = 1e4)
#' 
#' # Plot the histogram computed with the samples and the true value of the PMF
#' hist(samples,breaks=seq(0,n),freq=FALSE)
#' points(seq(0,n)-0.5,pmf_binomial,pch=16)
#' 
#' @export
PMF.sample = function(pmf, N_samples) {
  s = sample(0:(length(pmf)-1), prob = pmf, replace = TRUE, size = N_samples)
  return(s)
}

#' @title Get the mean of the distribution from a PMF object
#'
#' @description
#' 
#' Returns the mean from the PMF specified by `pmf`.
#' 
#' @param pmf the PMF object.
#' 
#' @examples
#' library(bayesRecon)
#' 
#' # Let's build the pmf of a Binomial distribution with parameters n and p
#' n <- 10
#' p <- 0.6 
#' pmf_binomial <- apply(matrix(seq(0,10)),MARGIN=1,FUN=function(x) dbinom(x,size=n,prob=p))
#' 
#' 
#' # The true mean corresponds to n*p
#' true_mean <- n*p
#' mean_from_PMF <- PMF.get_mean(pmf=pmf_binomial)
#' cat("True mean:", true_mean, "\nMean from PMF:", mean_from_PMF)
#' 
#' @return A numerical value for mean of the distribution.  
#' @seealso [PMF.get_var()], [PMF.get_quantile()], [PMF.sample()], [PMF.summary()]
#' @export
PMF.get_mean = function(pmf) {
  supp = 0:(length(pmf)-1)
  m = pmf %*% supp
  return(m)
}

#' @title Get the variance of the distribution from a PMF object
#'
#' @description
#' 
#' Returns the variance from the PMF specified by `pmf`.
#' 
#' @param pmf the PMF object.
#' 
#' @examples 
#' library(bayesRecon)
#' 
#' # Let's build the pmf of a Binomial distribution with parameters n and p
#' n <- 10
#' p <- 0.6 
#' pmf_binomial <- apply(matrix(seq(0,10)),MARGIN=1,FUN=function(x) dbinom(x,size=n,prob=p))
#' 
#' # The true variance corresponds to n*p*(1-p)
#' true_var <- n*p*(1-p)
#' var_from_PMF <- PMF.get_var(pmf=pmf_binomial)
#' cat("True variance:", true_var, "\nVariance from PMF:", var_from_PMF)
#' 
#' @return A numerical value for variance.  
#' @seealso [PMF.get_mean()], [PMF.get_quantile()], [PMF.sample()], [PMF.summary()]
#' @export
PMF.get_var = function(pmf) {
  supp = 0:(length(pmf)-1)
  v = pmf %*% supp^2 - PMF.get_mean(pmf)^2
  return(v)
}


#' @title Get quantile from a PMF object
#'
#' @description
#' 
#' Returns the `p` quantile from the PMF specified by `pmf`.
#' 
#' @param pmf the PMF object.
#' @param p the probability of the required quantile.
#' 
#' @examples 
#' library(bayesRecon)
#' 
#' # Let's build the pmf of a Binomial distribution with parameters n and p
#' n <- 10
#' p <- 0.6 
#' pmf_binomial <- apply(matrix(seq(0,10)),MARGIN=1,FUN=function(x) dbinom(x,size=n,prob=p))
#' 
#' # The true median is ceiling(n*p)
#' quant_50 <- PMF.get_quantile(pmf=pmf_binomial,p=0.5)
#' cat("True median:", ceiling(n*p), "\nMedian from PMF:", quant_50)
#'
#' 
#' @return A numeric value for the quantile.  
#' @seealso [PMF.get_mean()], [PMF.get_var()], [PMF.sample()], [PMF.summary()]
#' @export
PMF.get_quantile = function(pmf, p) {
  if (p <= 0 | p >= 1) {
    stop("Input error: probability p must be between 0 and 1")
  }
  cdf = cumsum(pmf)
  x = min(which(cdf >= p))
  return(x-1)
}


#' @title Returns summary of a PMF object
#'
#' @description
#' 
#' Returns the summary (min, max, IQR, median, mean) of the PMF specified by `pmf`.
#' 
#' @param pmf the PMF object.
#' @param Ltoll used for computing the min of the PMF: the min is the smallest value
#'        with probability greater than Ltoll (default: 1e-15) 
#' @param Rtoll used for computing the max of the PMF: the max is the largest value
#'        with probability greater than Rtoll (default: 1e-9)
#' 
#' @examples 
#' library(bayesRecon)
#' 
#' # Let's build the pmf of a Binomial distribution with parameters n and p
#' n <- 10
#' p <- 0.6 
#' pmf_binomial <- apply(matrix(seq(0,10)),MARGIN=1,FUN=function(x) dbinom(x,size=n,prob=p))
#' 
#' # Print the summary of this distribution
#' PMF.summary(pmf=pmf_binomial)
#'
#' 
#' @return A summary data.frame
#' @seealso [PMF.get_mean()], [PMF.get_var()], [PMF.get_quantile()], [PMF.sample()]
#' @export
PMF.summary = function(pmf, Ltoll=.TOLL, Rtoll=.RTOLL) {
  min_pmf = min(which(pmf > Ltoll)) - 1
  max_pmf = max(which(pmf > Rtoll)) - 1
  all_summaries <- data.frame("Min."    = min_pmf,
                              `1st Qu.` = PMF.get_quantile(pmf,0.25),
                              "Median"  = PMF.get_quantile(pmf,0.5),
                              "Mean"    = PMF.get_mean(pmf), 
                              `3rd Qu.` = PMF.get_quantile(pmf,0.75),
                              "Max"     = max_pmf,
                              check.names = FALSE)
  return(all_summaries)
}

# Apply smoothing to a pmf to "cover the holes" in the support.
# If there is no hole, it doesn't do anything.
# If the smoothing parameter alpha is not specified, it is set to the min of pmf. 
# If laplace is set to TRUE, add alpha to all the points. 
# Otherwise, add alpha only to points with zero mass.
PMF.smoothing = function(pmf, alpha = 1e-9, laplace=FALSE) {
  
  if (is.null(alpha)) alpha = min(pmf[pmf!=0])
  
  # apply smoothing only if there are holes
  if (sum(pmf==0)) {
    if (laplace) { pmf = pmf + rep(alpha, length(pmf))
    } else pmf[pmf==0] = alpha
  }
  
  return(pmf / sum(pmf))
}



# Compute convolution between 2 pmfs. Then, for numerical reasons: 
# -removes small values at the end of the vector (< Rtoll)
# -set to zero all the values to the left of the support
# -set to zero small values (< toll)
PMF.conv = function(pmf1, pmf2, toll=.TOLL, Rtoll=.RTOLL) {
  pmf = stats::convolve(pmf1, rev(pmf2), type="open")
  # Look for last value > Rtoll and remove all the elements after it:
  last_pos = max(which(pmf > Rtoll))  
  pmf = pmf[1:last_pos]
  # Set to zero values smaller than toll:
  pmf[pmf<toll] = 0  
  # Set to zero elements at the left of m1 + m2, which are the minima of the supports
  # of pmf1 and pmf2: guarantees that supp(v) is not "larger" than supp(v1) + supp(v2) 
  m1 = min(which(pmf1>0))
  m2 = min(which(pmf2>0))
  m = m1 + m2 -1
  if (m>1) pmf[1:(m-1)] = 0
  pmf = pmf / sum(pmf)  # re-normalize
  return(pmf)
}

# Computes the pmf of the bottom-up distribution analytically
# l_pmf: list of bottom pmfs
# toll and Rtoll: used during convolution (see PMF.conv)
# smoothing: whether to apply smoothing to the bottom pmfs to "cover the holes"  
# al_smooth, lap_smooth: smoothing parameters (see PMF.smoothing)
# Returns:
# -the bottom-up pmf, if return_all=FALSE
# -otherwise, a list of lists of pmfs for all the steps of the algorithm; 
#  they correspond to the variables of the "auxiliary binary tree"
PMF.bottom_up = function(l_pmf, toll=.TOLL, Rtoll=.RTOLL, return_all=FALSE,
                         smoothing=TRUE, al_smooth=.ALPHA_SMOOTHING, lap_smooth=.LAP_SMOOTHING) {
  
  # Smoothing to "cover the holes" in the supports of the bottom pmfs
  if (smoothing) l_pmf = lapply(l_pmf, PMF.smoothing, 
                                alpha=al_smooth, laplace=lap_smooth)
  
  # In case we have an upper which is a duplicate of a bottom,
  # the bottom up is simply that bottom.
  if(length(l_pmf)==1){
    if (return_all) {
      return(list(l_pmf))
    } else {
      return(l_pmf[[1]])
    }
  }
  
  # Doesn't do convolutions sequentially 
  # Instead, for each iteration (while) it creates a new list of vectors 
  # by doing convolution between 1 and 2, 3 and 4, ...
  # Then, the new list has length halved (if L is odd, just copy the last element)
  # Ends when the list has length 1: contains just 1 vector that is the convolution 
  # of all the vectors of the list 
  old_v = l_pmf
  l_l_v = list(old_v)   # list with all the step-by-step lists of pmf
  L = length(old_v)
  while (L > 1) {
    new_v = c()
    for (j in 1:(L%/%2)) {
      new_v = c(new_v, list(PMF.conv(old_v[[2*j-1]], old_v[[2*j]], 
                                     toll=toll, Rtoll=Rtoll)))
    }
    if (L%%2 == 1) new_v = c(new_v, list(old_v[[L]]))
    old_v = new_v
    l_l_v = c(l_l_v, list(old_v))
    L = length(old_v)
  }
  
  if (return_all) {
    return(l_l_v)
  } else {
    return(new_v[[1]])
  }
}

# Given a vector v_u and a list of bottom pmf l_pmf,
# checks if the elements of v_u are contained in the support of the bottom-up distr 
# Returns a vector with the same length of v_u with TRUE if it is contained and FALSE otherwise 
PMF.check_support = function(v_u, l_pmf, toll=.TOLL, Rtoll=.RTOLL,
                             smoothing=TRUE, al_smooth=.ALPHA_SMOOTHING, lap_smooth=.LAP_SMOOTHING) {
  pmf_u = PMF.bottom_up(l_pmf, toll=toll, Rtoll=Rtoll, return_all=FALSE,
                        smoothing=smoothing, al_smooth=al_smooth, lap_smooth=lap_smooth)
  # The elements of v_u must be in the support of pmf_u
  supp_u = which(pmf_u > 0) - 1
  mask = v_u %in% supp_u
  return(mask)
}

# Compute the tempered pmf
# The pmf is raised to the power of 1/temp, and then normalized
# temp must be a positive number
PMF.tempering = function(pmf, temp) {
  
  if (temp <= 0) stop("temp must be positive")
  if (temp == 1) return(pmf)
  
  temp_pmf = pmf**(1/temp)
  return(temp_pmf / sum(temp_pmf))
}
