test_that("Monthly, in_type=='params', distr='gaussian'",{
  # Run IS Reconc
  S = read.csv(file = "dataForTests/Monthly-Gaussian_S.csv", header = FALSE)
  S = as.matrix(S)
  base_forecasts_in = read.csv(file = "dataForTests/Monthly-Gaussian_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(as.numeric(base_forecasts_in[i,]))[[1]]
  }
  res.buis = reconc_BUIS(S, base_forecasts,
               in_type = "params", distr = "gaussian", num_samples = 100000, seed=42)
  # Run Gauss Reconc
  res.gauss = reconc_gaussian(S, base_forecasts_in[[1]], diag(base_forecasts_in[[2]]^2))
  # Test
  m = mean(rowMeans(res.buis$reconciled_samples) - as.numeric(rbind(res.gauss$upper_reconciled_mean, res.gauss$bottom_reconciled_mean)))
  expect_equal(abs(m) < 8e-3, TRUE)
})

test_that("Weekly, in_type=='params', distr='gaussian'",{
  # Run IS Reconc
  S = read.csv(file = "dataForTests/Weekly-Gaussian_S.csv", header = FALSE)
  S = as.matrix(S)
  base_forecasts_in = read.csv(file = "dataForTests/Weekly-Gaussian_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(as.numeric(base_forecasts_in[i,]))[[1]]
  }
  res.buis = reconc_BUIS(S, base_forecasts,
                  in_type = "params", distr = "gaussian", num_samples = 100000, seed=42)
  # Run Gauss Reconc
  res.gauss = reconc_gaussian(S, base_forecasts_in[[1]], diag(base_forecasts_in[[2]]^2))
  # Test
  m = mean(rowMeans(res.buis$reconciled_samples) - as.numeric(rbind(res.gauss$upper_reconciled_mean, res.gauss$bottom_reconciled_mean)))
  expect_equal(abs(m) < 2e-2, TRUE)
})

test_that("Monthly, in_type=='params', distr='poisson'",{
  S = read.csv(file = "dataForTests/Monthly-Poisson_S.csv", header = FALSE)
  S = as.matrix(S)
  base_forecasts_in = read.csv(file = "dataForTests/Monthly-Poisson_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(as.numeric(base_forecasts_in[i,]))[[1]]
  }
  res.buis = reconc_BUIS(S, base_forecasts,
                  in_type = "params", distr = "poisson", num_samples = 100000, seed=42)
  expect_no_error(res.buis)
})

test_that("Monthly, in_type=='samples', distr='continuous'",{
  # Run IS Reconc from samples
  S = read.csv(file = "dataForTests/Monthly-Gaussian_S.csv", header = FALSE)
  S = as.matrix(S)
  base_forecasts = .gen_gaussian("dataForTests/Monthly-Gaussian_basef.csv", seed=42)
  res.buis_samples = reconc_BUIS(S, base_forecasts, in_type = "samples", distr = "continuous", seed=42)
  # Run IS Reconc
  base_forecasts_in = read.csv(file = "dataForTests/Monthly-Gaussian_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(as.numeric(base_forecasts_in[i,]))[[1]]
  }
  res.buis = reconc_BUIS(S, base_forecasts, in_type = "params", distr = "gaussian", num_samples = 100000, seed=42)
  m = mean(rowMeans(res.buis$reconciled_samples) - rowMeans(res.buis_samples$reconciled_samples))
  expect_equal(abs(m) < 1e-2, TRUE)
})

test_that("Monthly, in_type=='samples', distr='discrete'",{
  # Run IS Reconc from samples
  S = read.csv(file = "dataForTests/Monthly-Poisson_S.csv", header = FALSE)
  S = as.matrix(S)
  base_forecasts = .gen_poisson("dataForTests/Monthly-Poisson_basef.csv", seed=42)
  res.buis_samples = reconc_BUIS(S, base_forecasts, in_type = "samples", distr = "discrete", seed=42)
  # Run IS Reconc
  base_forecasts_in = read.csv(file = "dataForTests/Monthly-Poisson_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(as.numeric(base_forecasts_in[i,]))[[1]]
  }
  res.buis = reconc_BUIS(S, base_forecasts, in_type = "params", distr = "poisson", num_samples = 100000, seed=42)
  m = mean(rowMeans(res.buis$reconciled_samples) - rowMeans(res.buis_samples$reconciled_samples))
  expect_equal(abs(m) < 1.5e-2, TRUE)
})

##############################################################################
