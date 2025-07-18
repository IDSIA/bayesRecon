test_that("Test shrinkage estimator", {
  # Parameters
  nSamples <- 500
  pTrue <- 2
  
  # True moments
  trueMean <- c(0,0)
  trueSigma <- matrix(c(3,2,2,2), nrow=2)
  chol_trueSigma <- chol(trueSigma)
  
  # we run 100 shrinkage estimators
  lambdas <- rep(0, 100)
  for(i in seq(100)){
    rr <- replicate(nSamples, trueMean) +  
      t(chol_trueSigma)%*%matrix(stats::rnorm(pTrue*nSamples), nrow=pTrue,ncol=nSamples)
    # Estimate mean and covariance from samples
    mean_est <- rowMeans(rr)
    Sigma_est <- cov(t(rr))
    
    mean_est
    trueMean
    
    Sigma_est
    
    rr_centered <- rr- replicate(nSamples, mean_est)
    rowMeans(rr_centered)
    
    x <- t(rr)
    #x <- t(rr_centered)
    
    lambdas[i] <- schaferStrimmer_cov(x)$lambda_star
    
  }
  # The average over 100 runs must be within a certain range
  mean_lambdas <- mean(lambdas)
  
  # We expect this to be lower than the 99.99% quantile over 200 simulations of this scenario
  expect_lte(mean_lambdas, 0.005191634)
  
  # We expect this to be greater than the 0.01% quantile over 200 simulations of this scenario
  expect_gte(mean_lambdas, 0.004843328)
  
})

test_that("Test behavior for p=1", {
  # Parameters
  nSamples <- 500
  pTrue <- 1
  
  samples <- matrix(stats::rnorm(n=nSamples,sd=2),ncol=pTrue)
  
  exp_sample_squared <- matrix(mean(samples^2),1,1)
  
  # The function should return E[x^2] and a warning
  expect_warning(expect_equal(schaferStrimmer_cov(samples)$shrink_cov, exp_sample_squared))
})
