test_that("MCMC Monthly, in_type=='params', distr='poisson'", {
  A = read.csv(file = "dataForTests/Monthly-Poisson_A.csv", header = FALSE)
  A = as.matrix(A)
  base_forecasts_in = read.csv(file = "dataForTests/Monthly-Poisson_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(lambda = base_forecasts_in[i,1])
  }
  res.buis = reconc_BUIS(A, base_forecasts,
                     in_type = "params", distr = "poisson", num_samples = 100000,seed=42)

  res.mcmc = reconc_MCMC(A, base_forecasts = base_forecasts, distr = "poisson", num_samples = 100000, seed=42)

  m = (rowMeans(res.buis$reconciled_samples) - rowMeans(res.mcmc$reconciled_samples))/rowMeans(res.buis$reconciled_samples)
  expect_equal(max(abs(m)) < 0.02, TRUE)
})
