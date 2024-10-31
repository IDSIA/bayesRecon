test_that("Monthly, in_type=='params', distr='gaussian'",{
  # Run IS Reconc
  A = read.csv(file = "dataForTests/Monthly-Gaussian_A.csv", header = FALSE)
  A = as.matrix(A)
  base_forecasts_in = read.csv(file = "dataForTests/Monthly-Gaussian_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(mean = base_forecasts_in[i,1],
                               sd = base_forecasts_in[i,2])
  }
  res.buis = reconc_BUIS(A, base_forecasts,
               in_type = "params", distr = "gaussian", num_samples = 100000, seed=42)
  # Run Gauss Reconc
  res.gauss = reconc_gaussian(A, base_forecasts_in[[1]], diag(base_forecasts_in[[2]]^2))
  # Test on bottom
  n_upper=nrow(A)
  n_bottom=ncol(A)
  m = mean(rowMeans(res.buis$reconciled_samples)[(n_upper+1):(n_upper+n_bottom)] - as.numeric(res.gauss$bottom_reconciled_mean))
  expect_equal(abs(m) < 8e-3, TRUE)
})

test_that("Weekly, in_type=='params', distr='gaussian'",{
  # Run IS Reconc
  A = read.csv(file = "dataForTests/Weekly-Gaussian_A.csv", header = FALSE)
  A = as.matrix(A)
  mode(A) <- "numeric"
  base_forecasts_in <- read.csv(file = "dataForTests/Weekly-Gaussian_basef.csv", header = FALSE)
  base_forecasts <- list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] <- list(mean = base_forecasts_in[i,1],
                               sd = base_forecasts_in[i,2])
  }
  res.buis <- reconc_BUIS(A, base_forecasts,
                  in_type = "params", distr = "gaussian", num_samples = 100000, seed=42)
  # Run Gauss Reconc
  res.gauss <- reconc_gaussian(A, base_forecasts_in[[1]], diag(base_forecasts_in[[2]]^2))
  # Test on bottom
  n_upper=nrow(A)
  n_bottom=ncol(A)
  m <- mean(rowMeans(res.buis$reconciled_samples)[(n_upper+1):(n_upper+n_bottom)] - as.numeric(res.gauss$bottom_reconciled_mean))
  expect_equal(abs(m) < 2e-2, TRUE)
})

test_that("Monthly, in_type=='params', distr='poisson'",{
  A = read.csv(file = "dataForTests/Monthly-Poisson_A.csv", header = FALSE)
  A = as.matrix(A)
  base_forecasts_in = read.csv(file = "dataForTests/Monthly-Poisson_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(lambda = base_forecasts_in[i,1])
  }
  res.buis = reconc_BUIS(A, base_forecasts,
                  in_type = "params", distr = "poisson", num_samples = 100000, seed=42)
  expect_no_error(res.buis)
})


test_that("Monthly, in_type=='params', distr='nbinom'",{
  A = read.csv(file = "dataForTests/Monthly-NegBin_A.csv", header = FALSE)
  A = as.matrix(A)
  base_forecasts_in = read.csv(file = "dataForTests/Monthly-NegBin_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(size = base_forecasts_in[i,2],
                               mu = base_forecasts_in[i,1])
  }
  res.buis = reconc_BUIS(A, base_forecasts,
                         in_type = "params", distr = "nbinom", num_samples = 10000, seed=42)
  expect_no_error(res.buis)
})



test_that("Monthly, in_type=='samples', distr='continuous'",{
  # Run IS Reconc from samples
  A = read.csv(file = "dataForTests/Monthly-Gaussian_A.csv", header = FALSE)
  A = as.matrix(A)
  base_forecasts = .gen_gaussian("dataForTests/Monthly-Gaussian_basef.csv", seed=42)
  res.buis_samples = reconc_BUIS(A, base_forecasts, in_type = "samples", distr = "continuous", seed=42)
  # Run IS Reconc
  base_forecasts_in = read.csv(file = "dataForTests/Monthly-Gaussian_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(mean = base_forecasts_in[i,1],
                               sd = base_forecasts_in[i,2])
  }
  res.buis = reconc_BUIS(A, base_forecasts, in_type = "params", distr = "gaussian", num_samples = 100000, seed=42)
  m = mean(rowMeans(res.buis$reconciled_samples) - rowMeans(res.buis_samples$reconciled_samples))
  expect_equal(abs(m) < 1e-2, TRUE)
})

test_that("Monthly, in_type=='samples', distr='discrete'",{
  # Run IS Reconc from samples
  A = read.csv(file = "dataForTests/Monthly-Poisson_A.csv", header = FALSE)
  A = as.matrix(A)
  base_forecasts = .gen_poisson("dataForTests/Monthly-Poisson_basef.csv", seed=42)
  res.buis_samples = reconc_BUIS(A, base_forecasts, in_type = "samples", distr = "discrete", seed=42)
  # Run IS Reconc
  base_forecasts_in = read.csv(file = "dataForTests/Monthly-Poisson_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(lambda = base_forecasts_in[i,1])
  }
  res.buis = reconc_BUIS(A, base_forecasts, in_type = "params", distr = "poisson", num_samples = 100000, seed=42)
  m = mean(rowMeans(res.buis$reconciled_samples) - rowMeans(res.buis_samples$reconciled_samples))
  expect_equal(abs(m) < 1.5e-2, TRUE)
})

test_that("Monthly simple, in_type=='params', distr='nbinom'",{

  # Read samples from dataForTests (reproducibility)
  vals <- read.csv(file = "dataForTests/Monthly-Count_ts.csv", header = FALSE)
  
  # Create a count time series with monthly observations for 10 years
  y <- ts(data=vals,frequency = 12)
  
  # Create the aggregated yearly time series
  y_agg <- temporal_aggregation(y,agg_levels = c(1,12))
  
  # We use a marginal forecast that computes for each month 
  # the empirical mean and variance
  # the forecast is a negative binomial with those params
  fc_bottom <- list()
  for(i in seq(12)){
    mm <- mean(y_agg$`f=12`[seq(i,120,12)])
    vv <- max(var(y_agg$`f=12`[seq(i,120,12)]), mm+0.5)
    #cat("i: ",i, "mean: ",mm, "var: ",vv, "size: ",mm^2/(vv-mm), "prob: ",mm/vv, "\n")
    
    fc_bottom[[i]] <- list(size=mm^2/(vv-mm),mu=mm)
  }
  
  # We compute the empirical mean and variance of the yearly ts
  # we forecast with a negative binomial with those parameters
  mm <- mean(y_agg$`f=1`)
  vv <- var(y_agg$`f=1`)
  fc_upper <- list(size=mm^2/(vv-mm), prob= mm/vv)
  
  # Obtain the aggregation matrix for this hierarchy
  rec_mat <- get_reconc_matrices(c(1,12),12)
  
  base_forecasts = append(list(fc_upper),fc_bottom)
  res.buis_params = reconc_BUIS(rec_mat$A, base_forecasts, in_type = "params", distr = "nbinom", seed=42)
  
  
  fc_upper_gauss <- list(mu=mm, Sigma = matrix(vv))
  res.mixCond <- reconc_MixCond(rec_mat$A, fc_bottom, fc_upper_gauss, bottom_in_type = "params", distr = 'nbinom')
  upp_pmf <- PMF.from_samples(as.integer(res.buis_params$upper_reconciled_samples))
  
  expect_equal(res.mixCond$upper_reconciled$pmf[[1]],upp_pmf,tolerance = 0.1)
  
})

##############################################################################
