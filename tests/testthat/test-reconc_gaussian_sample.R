test_that("Monthly, in_type=='samples', distr='continuous'", {

  # Run IS Reconc
  S = data.table::fread(file = "dataForTests/Monthly-Gaussian_S.csv", header = FALSE)
  S = as.matrix(S)
  base_forecasts_in = data.table::fread(file = "dataForTests/Monthly-Gaussian_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(as.numeric(base_forecasts_in[i,]))[[1]]
  }

  # Run Gauss Reconc
  res.gauss = bayesReco::reconc_gaussian(base_forecasts_in[[1]], diag(base_forecasts_in[[2]]^2), S)

  # Test


  res.gauss_sample = reconc_gaussian_sample(base_forecasts,
                                            S,
                                            cores = 1,
                                            chains = 12,
                                            iter = 2000,
                                            seed = 42)
  m = abs(rowMeans(res.gauss_sample$reconciled_samples) - as.numeric(rbind(res.gauss$upper_reconciled_mean, res.gauss$bottom_reconciled_mean)))
  m = max(m)


  expect_equal(m<5e-2, TRUE)
})
