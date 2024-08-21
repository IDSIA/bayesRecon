test_that("reconc_MixCond simple example", {
  
  # Simple example with 
  # - 12 bottom
  # - 10 upper: year, 6 bi-monthly, 3 quarterly
  A <- matrix(data=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                     1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
                     1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1),
              nrow=10,byrow = TRUE)
  
  # Define means and vars for the forecasts
  means <- c(90,62,63,64,31,32,31,33,31,32,rep(15,12))
  vars  <- c(20,8,8,8,4,4,4,4,4,4,rep(2,12))^2
  
  # create the lists for reconciliation
  ## upper 
  fc_upper <- list(mu = means[1:10],
                   Sigma = diag(vars[1:10]))
  
  ## bottom
  fc_bottom <- list()
  for(i in seq(ncol(A))){
    fc_bottom[[i]] <-as.integer(.distr_sample(list(mean=means[i+10],sd = vars[i+10]), "gaussian", 2e4))
    fc_bottom[[i]][which(fc_bottom[[i]]<0)] <- 0 # set-negative-to-zero
  }
  
  
  res.MixCond <- reconc_MixCond(A,fc_bottom,fc_upper,bottom_in_type = "samples",seed=42)
  
  bott_rec_means <- unlist(lapply(res.MixCond$bottom_reconciled$pmf,PMF.get_mean))
  bott_rec_vars <- unlist(lapply(res.MixCond$bottom_reconciled$pmf,PMF.get_var))
  
  
  # Create PMF from samples
  fc_bottom_pmf <- list()
  for(i in seq(ncol(A))){
    fc_bottom_pmf[[i]] <-PMF.from_samples(fc_bottom[[i]])
  }
  
  # Reconcile from bottom PMF
  res.MixCond_pmf <- reconc_MixCond(A,fc_bottom_pmf,fc_upper,seed=42)
  
  bott_rec_means_pmf <- unlist(lapply(res.MixCond_pmf$bottom_reconciled$pmf,PMF.get_mean))
  bott_rec_vars_pmf <- unlist(lapply(res.MixCond_pmf$bottom_reconciled$pmf,PMF.get_var))

  expect_equal(bott_rec_means,bott_rec_means_pmf,tolerance = "3e")
  expect_equal(bott_rec_vars,bott_rec_vars_pmf,tolerance = "3e")

})
