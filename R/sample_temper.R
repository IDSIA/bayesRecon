################################################################################
.RTOLL = 1e-9

# Compute the approximated PMF for a Gaussian distribution.
# The mass of point n is computed as P(n-0.5 < Gauss < n+0.5)
PMF.from_gauss = function(mu, sigma, Rtoll = .RTOLL) {
  M = stats::qnorm(1-Rtoll, mu, sigma)
  cdf = pnorm(0:M + 0.5, mu, sigma)
  pmf = c(cdf[1], diff(cdf))
  pmf = pmf / sum(pmf)
  return(pmf)
}

# Compute the product of 2 PMFs 
PMF.prod = function(pmf1, pmf2) {
  m = min(length(pmf1), length(pmf2))
  pmf = pmf1[1:m] * pmf2[1:m]
  pmf = pmf / sum(pmf)
  return(pmf)
}

# Sample from PMF, given a vector of uniform samples u
# Calling:
# u = runif(n_samp); s = PMF.sample_from_u(pmf, u)
# is equivalent to 
# s = PMF.sample(pmf, n_samp)
PMF.sample_from_u = function(pmf, u) {
  cdf = cumsum(pmf)
  samples = rowSums(outer(u, cdf, ">"))
  return(samples)
}


################################################################################
# Sample from distribution MVN(mu, Sigma) * pmf_1 * ... * pmf_n
# This is a distribution with 
# *i-th marginal = N(mu_i, Sigma_ii) * pmf_i 
# *Gaussian copula (specified by Sigma)
#
# First, sample from Gaussian copula; then, transform the marginals
# Output dim: n_samples x n
temper_sample = function(MVN_params, pmf_list, n_samples) {
  
  mu    = MVN_params$mu
  Sigma = MVN_params$Sigma
  
  n = length(pmf_list)
  # TODO: check that len(mu)=n and dim(Sigma)=(n,n)
  
  # Compute list of marginal PMFs
  gauss_marginal_pmfs = mapply(PMF.from_gauss, mu, diag(Sigma)**0.5)
  marginal_pmfs = mapply(PMF.prod, pmf_list, gauss_marginal_pmfs)
  
  # Sample from Gaussian copula
  R = stats::cov2cor(Sigma)         
  U = stats::pnorm(.MVN_sample(n_samples, rep(0,n), R))  # samples from Gaussian copula (dim: n_samples x n)
  
  # Get multivariate samples
  samples = mapply(PMF.sample_from_u, marginal_pmfs, asplit(U,2))
  samples = matrix(unlist(samples), nrow = n_samples)
  
  return(samples)
}


################################################################################
# Example

n = 10
mu = 1:n
sigma = seq(1, 2, length.out = n)
R = 0.3 * matrix(1, nrow=n, ncol=n)
diag(R) = 1
Sigma = R * outer(sigma, sigma, "*")
MVN_params = list(mu = mu, Sigma = Sigma)

lambdas = (1:n) / 2
pmf_list = lapply(lambdas, function(l) PMF.from_params(list(lambda=l), "poisson"))

n_samples = 1e5

samples = temper_sample(MVN_params, pmf_list, n_samples)

# Check correlation
cor(samples)
R
# They are different! Is it ok?

# Check marginals
# TODO






