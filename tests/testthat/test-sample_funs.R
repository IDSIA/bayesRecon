test_that("sampling from univariate normal", {
  
  #for(i in seq(5e4)){
  # Generate 1e4 samples from univariate Gaussian
  params <- list(mean=42, sd=1)
  distr <- "gaussian"
  n <- 1e4
  samples <- .distr_sample(params, distr, n)
  
  # Compute empirical mean and sd 
  sam_mean <- mean(samples)
  sam_sd <- sd(samples)
  
  # Check how close empirical values are to the truth
  m <- abs(sam_mean-42)/42
  s <- abs(sam_sd-1)
  #}
  #mean(s<4e-2)
  #mean(m<2e-3)
  
  expect_equal(m < 2e-3, TRUE)
  expect_equal(s < 4e-2, TRUE)
})

test_that("sampling from univariate nbinom", {
  
  #m <- seq(5e4)
  #for(i in seq(5e4)){
  # Generate 1e4 samples from negative binomial (size, prob)
  params <- list(size=12,prob=0.8)
  distr <- "nbinom"
  n <- 1e4
  samples <- .distr_sample(params, distr, n)
  
  # Compute empirical mean
  sam_mean <- mean(samples)
  true_mean <- params$size*(1-params$prob)/params$prob
  
  # Check how close empirical values are to the truth
  m <- abs(sam_mean-true_mean)/true_mean
  #}
  #mean(m<3e-2)

  expect_equal(m < 3e-2, TRUE)
  
  #m <- seq(5e4)
  #for(i in seq(5e4)){
  # Generate 1e4 samples from negative binomial (size, mu)
  params <- list(size=12,mu=true_mean)
  distr <- "nbinom"
  n <- 1e4
  samples <- .distr_sample(params, distr, n)
  
  # Compute empirical mean
  sam_mean <- mean(samples)
  
  # Check how close empirical values are to the truth
  m <- abs(sam_mean-params$mu)/params$mu
  #}
  #mean(m<3e-2)
  
  expect_equal(m < 3e-2, TRUE)
})

test_that("sampling from univariate poisson", {
  
  #m <- seq(5e4)
  #for(i in seq(5e4)){
  # Generate 1e4 samples from poisson
  params <- list(lambda=10)
  distr <- "poisson"
  n <- 1e4
  samples <- .distr_sample(params, distr, n)
  
  # Compute empirical mean
  sam_mean <- mean(samples)
  
  # Check how close empirical values are to the truth
  m <- abs(sam_mean-10)/10
  #}
  
  #mean(m<3e-2)
  
  expect_equal(m < 3e-2, TRUE)
})

test_that("sampling from multivariate normal", {
  #m <- matrix(nrow=5e4,ncol=2)
  #for(i in seq(5e4)){
  # Generate 1e4 samples from bivariate normal
  mu=c(10,10)
  Sigma= matrix(c(1,0.7,0.7,1),nrow = 2)
  n <- 1e4
  samples <- .MVN_sample(n, mu, Sigma)
  
  # Compute empirical mean
  sam_mean <- colMeans(samples)
  
  # Check how close empirical values are to the truth
  m <- abs(sam_mean-10)/10
  #}
  
  #all(colMeans(m<8e-3)==1)
  
  expect_equal(all(m < 8e-3), TRUE)
})
