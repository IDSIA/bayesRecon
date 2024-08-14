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
temper_sample = function(MVN_params, pmf_list, num_samples) {
  
  mu    = MVN_params$mu
  Sigma = MVN_params$Sigma
  
  n = length(pmf_list)
  
  # CHECK ON MU, SIGMA
  # If Sigma is a number, transform into a matrix 
  if (length(Sigma) == 1) Sigma = as.matrix(Sigma)
  # Check the dimensions of mu and Sigma
  if (length(mu) != n | any(dim(Sigma) != c(n, n))) {
    stop("Input error: the dimensions of the upper parameters do not match with S")
  }
  # Check that Sigma is a covariance matrix (symmetric positive semi-definite)
  .check_cov(Sigma, "Covariance matrix", symm_check=TRUE)
  
  # Compute list of marginal PMFs
  gauss_marginal_pmfs = mapply(PMF.from_gauss, mu, diag(Sigma)**0.5)  # discretize gaussian marginal
  marginal_pmfs = mapply(PMF.prod, pmf_list, gauss_marginal_pmfs)     # multiply by pmf in pmf_list
  
  # Sample from Gaussian copula
  R = stats::cov2cor(Sigma)         
  U = stats::pnorm(.MVN_sample(num_samples, rep(0,n), R))  # samples from Gaussian copula (dim: n_samples x n)
  
  # Get multivariate samples
  samples = mapply(PMF.sample_from_u, marginal_pmfs, asplit(U,2))
  samples = matrix(unlist(samples), nrow = num_samples)
  
  return(samples)
}


################################################################################
# # Example
# # 
# n = 10
# mu = 1:n
# sigma = seq(1, 2, length.out = n)
# R = 0.3 * matrix(1, nrow=n, ncol=n)
# diag(R) = 1
# Sigma = R * outer(sigma, sigma, "*")
# MVN_params = list(mu = mu, Sigma = Sigma)
# 
# lambdas = (1:n) / 2
# pmf_list = lapply(lambdas, function(l) PMF.from_params(list(lambda=l), "poisson"))
# 
# n_samples = 1e5
# 
# samples = temper_sample(MVN_params, pmf_list, n_samples)
# 
# # Check correlation
# cor(samples)
# R
# # They are different! Is it ok?
# Sp_rho = matrix(nrow = n, ncol = n)
# for (i in 1:n) {
#   for (j in 1:n) {
#     Sp_rho[i,j] = cor.test(samples[,i], samples[,j], method = "spearman")$estimate
#   }
# }

# # Check marginals
# for (i in 1:n) {
#   print(i)
#   
#   pmf_g = PMF.from_gauss(mu[i], sigma[i])
#   pmf = PMF.prod(pmf_g, pmf_list[[i]])
#   pmf_emp = PMF.from_samples(samples[,i])
#   
#   print(PMF.summary(pmf))
#   print(PMF.summary(pmf_emp))
# }









