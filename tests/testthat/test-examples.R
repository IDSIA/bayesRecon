test_that("Monthly, in_type=='params', distr='gaussian'",{
  S = read.csv(file = "dataForTests/Monthly-Gaussian_S.csv", header = FALSE)
  S = as.matrix(S)
  base_forecasts_in = read.csv(file = "dataForTests/Monthly-Gaussian_basef.csv", header = FALSE)
  base_forecasts = list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] = list(as.numeric(base_forecasts_in[i,]))[[1]]
  }
  res = reconc(S, base_forecasts,
               in_type = "params", distr = "gaussian", num_samples = 100000)
  expect_no_error(res)
})



