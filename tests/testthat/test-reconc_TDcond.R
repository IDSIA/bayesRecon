test_that("reconc_TDcond simple example", {

  # Simple example with 
  # - 12 bottom
  # - 10 upper: year, 6 bi-monthly, 3 quarterly
  A <- matrix(data=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                     1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
                     1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1),
              nrow=10,byrow = TRUE)
  
  
  S <- rbind(A,diag(nrow=12))
  
  # Define means and vars for the forecasts
  means <- c(90,31,32,31,33,31,32,62,63,64,rep(15,12))
  vars  <- c(20,4,4,4,4,4,4,8,8,8,rep(2,12))^2
  
  # create the lists for reconciliation
  ## upper 
  fc_upper <- list(mu = means[1:10],
                   Sigma = diag(vars[1:10]))
  
  ## bottom
  fc_bottom <- list()
  for(i in seq(ncol(S))){
    fc_bottom[[i]] <-as.integer(.distr_sample(list(mean=means[i+10],sd = vars[i+10]), "gaussian", 2e4))
    fc_bottom[[i]][which(fc_bottom[[i]]<0)] <- 0 # set-negative-to-zero
  }
  

  # reconciliation with TDcond
  res.TDcond <- reconc_TDcond(S, fc_bottom, fc_upper,
                bottom_in_type = "samples",
                num_samples = 2e4, return_type = "pmf",
                seed = 42)
  
  res.TDcond2 <- reconc_TDcond(S, fc_bottom, fc_upper,
                              bottom_in_type = "samples",
                              num_samples = 2e4, return_type = "samples",
                              seed = 42)
  
  res.TDcond3 <- reconc_TDcond(S, fc_bottom, fc_upper,
                              bottom_in_type = "samples",
                              num_samples = 2e4, return_type = "all",
                              seed = 42)
  
  # Check if all return_type return identical results
  expect_identical(res.TDcond$bottom_reconciled$pmf,res.TDcond3$bottom_reconciled$pmf)
  expect_identical(res.TDcond2$bottom_reconciled$samples,res.TDcond3$bottom_reconciled$samples)
  
  # Compute the reconciliation analytically (everything Gaussian)
  ## Save bottom forecast parameters
  fc_bott_gauss <- list(mu = means[11:22],
                        Sigma = diag(vars[11:22]))
  
  # Compute the reconciled precision  
  inv_B <- diag(1/diag(fc_bott_gauss$Sigma))
  inv_U <- diag(1/diag(fc_upper$Sigma))
  At_inv_U_A <- crossprod(A,inv_U)%*%A
  # Here we use the reduced A with only the lowest level
  Au <- A[.lowest_lev(A),]
  inv_A_B_At <- solve(Au%*%tcrossprod(fc_bott_gauss$Sigma,Au))
  
  # formulas for the reconciled precision, covariance and mean 
  bott_reconc_Prec <- inv_B+At_inv_U_A-crossprod(Au,inv_A_B_At)%*%Au
  bott_reconc_cov  <- solve(bott_reconc_Prec)
  bott_reconc_mean <- fc_bott_gauss$mu + tcrossprod(bott_reconc_cov,A)%*%inv_U%*%(fc_upper$mu-A%*%fc_bott_gauss$mu)
  
  # compute the difference between empirical and analytical
  m_diff <- unlist(lapply(res.TDcond$bottom_reconciled$pmf,PMF.get_mean)) - bott_reconc_mean
  
  expect_true(all(abs(m_diff/bott_reconc_mean)<8e-3))
  
  
  # The variances are different
  #unlist(lapply(res.TDcond$pmf$bottom,PMF.get_var))
  #diag(bott_reconc_cov)

})
