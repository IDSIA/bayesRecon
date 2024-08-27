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

##############################################################################
