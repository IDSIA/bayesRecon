test_that("Test effective sample size", {
  S = matrix(data = c(1,1,1,0,0,1), nrow=3, byrow = TRUE)
  
  # -----------
  n = 200
  b1 = stats::rpois(n=n, lambda = 3)
  b2 = stats::rpois(n=n, lambda = 4)
  u = stats::rnorm(n=n, mean = 30, sd = 1)
  B = cbind(b1,b2)
  c = matrix(S[1,])
  b = (B %*% c)
  
  w = .compute_weights(b, u, "samples", "continuous")
  
  check_w = .check_weigths(w, n_eff_min=200)
  expect_equal(check_w$warning, TRUE)
  expect_equal(check_w$warning_code, 1)
  expect_equal(check_w$n_eff, n)
  
  # Try the warning message
  # base_forecast = list(u,b1,b2)
  # a = reconc_BUIS(S, base_forecast, in_type = "samples", distr = list("continuous","discrete","discrete"), seed=42)
  
  # -----------
  n = 199
  b1 = stats::rpois(n=n, lambda = 3)
  b2 = stats::rpois(n=n, lambda = 4)
  u = stats::rnorm(n=n, mean = 30, sd = 1)
  B = cbind(b1,b2)
  c = matrix(S[1,])
  b = (B %*% c)
  
  w = .compute_weights(b, u, "samples", "continuous")
  
  check_w = .check_weigths(w, n_eff_min=200)
  expect_equal(check_w$warning, TRUE)
  expect_equal(check_w$warning_code, 1)
  expect_equal(check_w$n_eff, n)
  
  # Try the warning message
  # base_forecast = list(u,b1,b2)
  # a = bayesRecon::reconc_BUIS(S, base_forecast, in_type = "samples", distr = list("continuous","discrete","discrete"), seed=42)

  # -----------
  n = 2000
  b1 = stats::rpois(n=n, lambda = 3)
  b2 = stats::rpois(n=n, lambda = 4)
  u = stats::rnorm(n=n, mean = 18, sd = 1)
  B = cbind(b1,b2)
  c = matrix(S[1,])
  b = (B %*% c)
  
  w = .compute_weights(b, u, "samples", "continuous")
  
  check_w = .check_weigths(w, n_eff_min=200, p_n_eff=0.01)
  expect_equal(check_w$warning, TRUE)
  expect_equal(check_w$warning_code, c(2,3))
  expect_equal(check_w$n_eff < 200, TRUE)
  
  # Try the warning message
  # base_forecast = list(u,b1,b2)
  # a = bayesRecon::reconc_BUIS(S, base_forecast, in_type = "samples", distr = list("continuous","discrete","discrete"), seed=42)
})
