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

  expect_equal(bott_rec_means,bott_rec_means_pmf, tolerance = 0.01)
  expect_equal(bott_rec_vars,bott_rec_vars_pmf, tolerance = 0.1)
  
})

test_that("reconc_MixCond and reconc_TDcond with temporal hier and params", {
  

  # Read samples from dataForTests (reproducibility)
  vals <- read.csv(file = "dataForTests/Monthly-Count_ts.csv", header = FALSE)
  
  # Create a count time series with monthly observations for 10 years
  y <- ts(data=vals,frequency = 12)
  
  # Create the aggregated yearly time series
  y_agg <- temporal_aggregation(y,agg_levels = c(1,12))
  
  # We use a marginal forecast that computes for each month 
  # the empirical mean and forecasts a Poisson with that value
  fc_bottom <- list()
  for(i in seq(12)){
    fc_bottom[[i]] <- list(lambda=mean(y_agg$`f=12`[seq(i,120,12)]))
  }
  
  # We compute the empirical mean and variance of the yearly ts
  # we forecast with a Gaussian with those parameters
  fc_upper <- list(mu=mean(y_agg$`f=1`), Sigma=matrix(var(y_agg$`f=1`)))
  
  # Obtain the aggregation matrix for this hierarchy
  rec_mat <- get_reconc_matrices(c(1,12),12)
  
  # Do a couple of checks on S and A
  expect_no_error(.check_S(rec_mat$S))
  expect_error(.check_S(rec_mat$A))
  expect_true(.check_BU_matr(rec_mat$A))
  expect_false(.check_BU_matr(rec_mat$S))
    
  # We can reconcile with reconc_MixCond
  res.mixCond <- reconc_MixCond(rec_mat$A, fc_bottom, fc_upper, bottom_in_type = "params", distr = 'poisson')
  
  # We can reconcile with reconc_TDcond
  res.TDcond <- reconc_TDcond(rec_mat$A, fc_bottom, fc_upper, bottom_in_type = "params", distr = 'poisson')
  
  # Summary of the upper reconciled with TDcond 
  pmfSum <- PMF.summary(res.TDcond$upper_reconciled$pmf[[1]])
  # We expect that the reconciled mean is very similar to the initial mean (should be equal)
  expect_equal(pmfSum$Mean,fc_upper$mu,tolerance = 0.01)
  
  # Check that all bottom and upper reconciled PMF sum to 1
  check_pmf_bott_mixCond <- sum(unlist(lapply(res.mixCond$bottom_reconciled$pmf, function(x){sum(x)}))) 
  check_pmf_upp_mixCond <- sum(unlist(lapply(res.mixCond$upper_reconciled$pmf, function(x){sum(x)}))) 
  expect_equal(check_pmf_bott_mixCond,12)
  expect_equal(check_pmf_upp_mixCond,1)
  
  # Check that all bottom and upper reconciled PMF sum to 1
  check_pmf_bott_TDcond <- sum(unlist(lapply(res.TDcond$bottom_reconciled$pmf, function(x){sum(x)}))) 
  check_pmf_upp_TDcond <- sum(unlist(lapply(res.TDcond$upper_reconciled$pmf, function(x){sum(x)}))) 
  expect_equal(check_pmf_bott_TDcond,12)
  expect_equal(check_pmf_upp_TDcond,1)

  
})
