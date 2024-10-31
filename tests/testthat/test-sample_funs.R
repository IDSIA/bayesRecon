test_that("sampling from univariate normal", {

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

  expect_equal(m < 2e-3, TRUE)
  expect_equal(s < 4e-2, TRUE)
})

test_that("sampling from univariate nbinom", {

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

  expect_equal(m < 3e-2, TRUE)

  # Generate 1e4 samples from negative binomial (size, mu)
  params <- list(size=12,mu=true_mean)
  distr <- "nbinom"
  n <- 1e4
  samples <- .distr_sample(params, distr, n)

  # Compute empirical mean
  sam_mean <- mean(samples)

  # Check how close empirical values are to the truth
  m <- abs(sam_mean-params$mu)/params$mu

  expect_equal(m < 3e-2, TRUE)

  # Check if it returns an error when all 3 parameters are specified
  params <- list(size=12,mu=true_mean,prob=0.8)
  distr <- "nbinom"
  n <- 1e4
  expect_error(.distr_sample(params, distr, n))

  # Check if it returns an error when size is not specified
  params <- list(mu=true_mean,prob=0.8)
  distr <- "nbinom"
  n <- 1e4
  expect_error(.distr_sample(params, distr, n))



})

test_that("sampling from univariate poisson", {

  # Generate 1e4 samples from poisson
  params <- list(lambda=10)
  distr <- "poisson"
  n <- 1e4
  samples <- .distr_sample(params, distr, n)

  # Compute empirical mean
  sam_mean <- mean(samples)

  # Check how close empirical values are to the truth
  m <- abs(sam_mean-10)/10

  expect_equal(m < 3e-2, TRUE)
})

test_that("sampling from multivariate normal", {

  # Generate 1e4 samples from bivariate normal
  mu=c(10,10)
  Sigma= matrix(c(1,0.7,0.7,1),nrow = 2)
  n <- 1e4
  samples <- .MVN_sample(n, mu, Sigma)

  # Compute empirical mean
  sam_mean <- colMeans(samples)

  # Check how close empirical values are to the truth
  m <- abs(sam_mean-10)/10

  expect_equal(all(m < 8e-3), TRUE)
})

test_that("MVN density works", {

  # Create 3x3 covariance matrix
  L <- matrix(0,nrow=3,ncol=3)
  L[lower.tri(L,diag=TRUE)] <- c(0.9,0.8,0.5,0.9,0.2,0.6)
  Sigma <- L%*%t(L)

  # create mean vector
  mu <- c(0,1,-1)

  # matrix where to evaluate the MVN
  xx <- matrix(c(0,2,1,
                 2,3,4,
                 0.5,0.5,0.5,
                 0,1,-1), ncol=3,byrow=TRUE)
  
  res <- .MVN_density(x=xx,mu=mu,Sigma=Sigma)

  true_val <- c(8.742644e-04, 1.375497e-11, 3.739985e-03, 1.306453e-01)
  expect_equal(res,true_val, tolerance = 1e-6)

  # Check if block-evaluation works
  xx <- matrix(runif(3*1e4),ncol=3,byrow=TRUE)

  res_chuncks <- .MVN_density(x=xx,mu=mu,Sigma=Sigma)
  res_all <- .MVN_density(x=xx,mu=mu,Sigma=Sigma,max_size_x = 1e4)

  expect_equal(res_chuncks,res_all)

})
