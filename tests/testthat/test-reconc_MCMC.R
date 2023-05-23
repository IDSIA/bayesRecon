test_that("MCMC Monthly, in_type=='params', distr='poisson'", {
  S = read.csv(file = "dataForTests/Monthly-Poisson_S.csv", header = FALSE)
  S = as.matrix(S)
  base_forecasts_in = read.csv(file = "dataForTests/Monthly-Poisson_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(as.numeric(base_forecasts_in[i,]))[[1]]
  }
  res.buis = reconc_BUIS(S, base_forecasts,
                     in_type = "params", distr = "poisson", num_samples = 100000,seed=42)

  res.mcmc = reconc_MCMC(S, base_forecasts = base_forecasts, distr = "poisson", num_samples = 100000, seed=42)

  m = (rowMeans(res.buis$reconciled_samples) - rowMeans(res.mcmc$reconciled_samples))/rowMeans(res.buis$reconciled_samples)
  expect_equal(max(abs(m)) < 0.02, TRUE)
})
