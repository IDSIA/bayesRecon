test_that("Monthly, in_type=='params', distr='gaussian'",{
  # Run IS Reconc
  S = data.table::fread(file = "dataForTests/Monthly-Gaussian_S.csv", header = FALSE)
  S = as.matrix(S)
  base_forecasts_in = data.table::fread(file = "dataForTests/Monthly-Gaussian_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(as.numeric(base_forecasts_in[i,]))[[1]]
  }
  res.is = reconc(S, base_forecasts,
               in_type = "params", distr = "gaussian", num_samples = 100000)
  # Run Gauss Reconc
  res.gauss = reconc_gaussian(base_forecasts_in[,1], diag(base_forecasts_in[,2]^2), S)
  # Test
  m = mean(colMeans(res.is$reconciled_samples) - as.numeric(rbind(res.gauss$upper_reconciled_mean, res.gauss$bottom_reconciled_mean)))
  expect_equal(abs(m) < 1e-3, TRUE)
})

test_that("Weekly, in_type=='params', distr='gaussian'",{
  # Run IS Reconc
  S = data.table::fread(file = "dataForTests/Weekly-Gaussian_S.csv", header = FALSE)
  S = as.matrix(S)
  base_forecasts_in = data.table::fread(file = "dataForTests/Weekly-Gaussian_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(as.numeric(base_forecasts_in[i,]))[[1]]
  }
  res.is = reconc(S, base_forecasts,
                  in_type = "params", distr = "gaussian", num_samples = 100000)
  # Run Gauss Reconc
  res.gauss = reconc_gaussian(base_forecasts_in[,1], diag(base_forecasts_in[,2]^2), S)
  # Test
  m = mean(colMeans(res.is$reconciled_samples) - as.numeric(rbind(res.gauss$upper_reconciled_mean, res.gauss$bottom_reconciled_mean)))
  expect_equal(abs(m) < 1e-3, TRUE)
})

test_that("Monthly, in_type=='params', distr='poisson'",{
  S = data.table::fread(file = "dataForTests/Monthly-Poisson_S.csv", header = FALSE)
  S = as.matrix(S)
  base_forecasts_in = data.table::fread(file = "dataForTests/Monthly-Poisson_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(as.numeric(base_forecasts_in[i,]))[[1]]
  }
  res.is = reconc(S, base_forecasts,
                  in_type = "params", distr = "poisson", num_samples = 100000)
  expect_no_error(res.is)
})

test_that("Monthly, in_type=='samples', distr='continuous'",{
  # Run IS Reconc from samples
  S = data.table::fread(file = "dataForTests/Monthly-Gaussian_S.csv", header = FALSE)
  S = as.matrix(S)
  base_forecasts_in = data.table::fread(file = "dataForTests/Monthly-Gaussian_basef-samples.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(as.numeric(base_forecasts_in[i,]))[[1]]
  }
  res.is_samples = reconc(S, base_forecasts, in_type = "samples", distr = "continuous")
  # Run IS Reconc
  base_forecasts_in = data.table::fread(file = "dataForTests/Monthly-Gaussian_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(as.numeric(base_forecasts_in[i,]))[[1]]
  }
  res.is = reconc(S, base_forecasts, in_type = "params", distr = "gaussian", num_samples = 100000)
  m = mean(colMeans(res.is$reconciled_samples) - colMeans(res.is_samples$reconciled_samples))
  expect_equal(abs(m) < 0.5, TRUE) # Suspicious...
})

test_that("Monthly, in_type=='samples', distr='discrete'",{
  # Run IS Reconc from samples
  S = data.table::fread(file = "dataForTests/Monthly-Poisson_S.csv", header = FALSE)
  S = as.matrix(S)
  base_forecasts_in = data.table::fread(file = "dataForTests/Monthly-Poisson_basef-samples.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(as.numeric(base_forecasts_in[i,]))[[1]]
  }
  res.is_samples = reconc(S, base_forecasts, in_type = "samples", distr = "discrete")
  # Run IS Reconc
  base_forecasts_in = data.table::fread(file = "dataForTests/Monthly-Poisson_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(as.numeric(base_forecasts_in[i,]))[[1]]
  }
  res.is = reconc(S, base_forecasts, in_type = "params", distr = "poisson", num_samples = 100000)
  m = mean(colMeans(res.is$reconciled_samples) - colMeans(res.is_samples$reconciled_samples))
  expect_equal(abs(m) < 1e-2, TRUE)
})

##############################################################################
