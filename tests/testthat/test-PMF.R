test_that("PMF.conv", {
  
  # Generate a size and a prob each binomial
  s1 <- 20; p1 <- 0.6
  s2 <- 30; p2 <- 0.7
  
  # Compute the pmf for the two binomials
  pmf1 <- dbinom(seq(0,s1),size=s1,prob=p1)
  pmf2 <- dbinom(seq(0,s2),size=s2,prob=p2)
  
  # True mean of the convolution
  expected_conv_mean = s1*p1+s2*p2
  expected_conv_var = s1*p1*(1-p1)+s2*p2*(1-p2)
  
  # Compute the PMF of the convolution
  pmf_conv <- PMF.conv(pmf1 = pmf1,pmf2 = pmf2)
  
  abs(PMF.get_mean(pmf_conv)-expected_conv_mean)
  abs(PMF.get_var(pmf_conv)-expected_conv_var)
  
  # Check if the convolution mean is equal to the truth
  expect_lt(abs(PMF.get_mean(pmf_conv)-expected_conv_mean),1e-6)
  
  # Check if the convolution var is equal to the truth
  expect_lt(abs(PMF.get_var(pmf_conv)-expected_conv_var),8e-5)
})

test_that("PMF.bottom_up", {
  
  # Test with 10 bottom
  n_bottom <- 10
  
  # Create sizes and probs for nbinom bottom distributions
  sizes <- c(seq(11,15),seq(19,15))
  probs <- c(rep(0.4,5),rep(0.7,5))
  distr <- "nbinom"
  
  # Compute true BU params (mean/var) and bottom pmfs
  true_bu_mean <- 0
  true_bu_var <- 0
  bott_pmfs <- list()
  for(i in seq(n_bottom)){
    true_bu_mean <- true_bu_mean + sizes[i]*(1-probs[i])/probs[i]
    true_bu_var <- true_bu_var + sizes[i]*(1-probs[i])/probs[i]^2
    params <- list(size=sizes[i],prob=probs[i])
    bott_pmfs[[i]] <- PMF.from_params(params=params,distr = distr)
  }
  # Run PMF.bottom_up to get the bottom up pmf
  bottom_up_pmf <- PMF.bottom_up(l_pmf=bott_pmfs)
  
  # Check if true mean is close enough to bottom up pmf mean
  expect_lt(abs(PMF.get_mean(bottom_up_pmf)-true_bu_mean)/true_bu_mean,2e-6)
  
  # Check if true var is close enough to bottom up pmf var
  expect_lt(abs(PMF.get_var(bottom_up_pmf)-true_bu_var)/true_bu_var,6e-5)
})

test_that("PMF.quantile",{
  n_samples = 1e5
  size = 10
  prob = 0.6
  p = 0.01
  
  x = rnbinom(n_samples, size = size, prob = prob)
  pmf = PMF.from_samples(x)
  q = PMF.get_quantile(pmf, p)
  qq = qnbinom(p, size = size, prob = prob)
  
  expect_equal(q,qq)
})
