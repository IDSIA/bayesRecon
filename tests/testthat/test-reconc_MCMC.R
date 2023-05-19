test_that("MCMC Monthly, in_type=='params', distr='poisson'", {
  S = data.table::fread(file = "dataForTests/Monthly-Poisson_S.csv", header = FALSE)
  S = as.matrix(S)
  base_forecasts_in = data.table::fread(file = "dataForTests/Monthly-Poisson_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(as.numeric(base_forecasts_in[i,]))[[1]]
  }
  res.is = reconc_IS(S, base_forecasts,
                     in_type = "params", distr = "poisson", num_samples = 100000,seed=42)

  res.mcmc = reconc_MCMC(S,params = base_forecasts,distr = rep("poisson",28),N_samples = 100000,seed=42)

  m = (colMeans(res.is$reconciled_samples) - colMeans(res.mcmc$reconciled_samples))/colMeans(res.is$reconciled_samples)
  expect_equal(max(abs(m)) < 0.02, TRUE)
})
