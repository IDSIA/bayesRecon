test_that("Monthly, in_type=='params', distr='gaussian'", {
  # Run IS Reconc
  A <- read.csv(file = "dataForTests/Monthly-Gaussian_A.csv", header = FALSE)
  A <- as.matrix(A)
  base_forecasts_in <- read.csv(file = "dataForTests/Monthly-Gaussian_basef.csv", header = FALSE)
  base_forecasts <- list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] <- list(
      mean = base_forecasts_in[i, 1],
      sd = base_forecasts_in[i, 2]
    )
  }
  res.buis <- reconc_BUIS(A, base_forecasts,
    in_type = "params", distr = "gaussian", num_samples = 100000, seed = 42
  )
  # Run Gauss Reconc
  res.gauss <- reconc_gaussian(A, base_forecasts_in[[1]], diag(base_forecasts_in[[2]]^2))
  # Test on bottom
  n_upper <- nrow(A)
  n_bottom <- ncol(A)
  m <- mean(rowMeans(res.buis$bottom_rec_samples) - as.numeric(res.gauss$bottom_rec_mean))
  expect_equal(abs(m) < 8e-3, TRUE)
})

test_that("Weekly, in_type=='params', distr='gaussian'", {
  # Run IS Reconc
  A <- read.csv(file = "dataForTests/Weekly-Gaussian_A.csv", header = FALSE)
  A <- as.matrix(A)
  mode(A) <- "numeric"
  base_forecasts_in <- read.csv(file = "dataForTests/Weekly-Gaussian_basef.csv", header = FALSE)
  base_forecasts <- list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] <- list(
      mean = base_forecasts_in[i, 1],
      sd = base_forecasts_in[i, 2]
    )
  }
  res.buis <- reconc_BUIS(A, base_forecasts,
    in_type = "params", distr = "gaussian", num_samples = 100000, seed = 42
  )
  # Run Gauss Reconc
  res.gauss <- reconc_gaussian(A, base_forecasts_in[[1]], diag(base_forecasts_in[[2]]^2))
  # Test on bottom
  n_upper <- nrow(A)
  n_bottom <- ncol(A)
  m <- mean(rowMeans(res.buis$bottom_rec_samples) - as.numeric(res.gauss$bottom_rec_mean))
  expect_equal(abs(m) < 2e-2, TRUE)
})

test_that("Monthly, in_type=='params', distr='poisson'", {
  A <- read.csv(file = "dataForTests/Monthly-Poisson_A.csv", header = FALSE)
  A <- as.matrix(A)
  base_forecasts_in <- read.csv(file = "dataForTests/Monthly-Poisson_basef.csv", header = FALSE)
  base_forecasts <- list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] <- list(lambda = base_forecasts_in[i, 1])
  }
  res.buis <- reconc_BUIS(A, base_forecasts,
    in_type = "params", distr = "poisson", num_samples = 100000, seed = 42
  )
  expect_no_error(res.buis)
})


test_that("Monthly, in_type=='params', distr='nbinom'", {
  A <- read.csv(file = "dataForTests/Monthly-NegBin_A.csv", header = FALSE)
  A <- as.matrix(A)
  base_forecasts_in <- read.csv(file = "dataForTests/Monthly-NegBin_basef.csv", header = FALSE)
  base_forecasts <- list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] <- list(
      size = base_forecasts_in[i, 2],
      mu = base_forecasts_in[i, 1]
    )
  }
  res.buis <- reconc_BUIS(A, base_forecasts,
    in_type = "params", distr = "nbinom", num_samples = 10000, seed = 42
  )
  expect_no_error(res.buis)
})


test_that("Monthly, in_type=='samples', distr='continuous'", {
  # Run IS Reconc from samples
  A <- read.csv(file = "dataForTests/Monthly-Gaussian_A.csv", header = FALSE)
  A <- as.matrix(A)
  base_forecasts <- .gen_gaussian("dataForTests/Monthly-Gaussian_basef.csv", seed = 42)
  res.buis_samples <- reconc_BUIS(A, base_forecasts, in_type = "samples", distr = "continuous", seed = 42)
  # Run IS Reconc
  base_forecasts_in <- read.csv(file = "dataForTests/Monthly-Gaussian_basef.csv", header = FALSE)
  base_forecasts <- list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] <- list(
      mean = base_forecasts_in[i, 1],
      sd = base_forecasts_in[i, 2]
    )
  }
  res.buis <- reconc_BUIS(A, base_forecasts, in_type = "params", distr = "gaussian", num_samples = 100000, seed = 42)
  m <- mean(rowMeans(res.buis$bottom_rec_samples) - rowMeans(res.buis_samples$bottom_rec_samples))
  expect_equal(abs(m) < 1e-2, TRUE)
})

test_that("Monthly, in_type=='samples', distr='discrete'", {
  # Run IS Reconc from samples
  A <- read.csv(file = "dataForTests/Monthly-Poisson_A.csv", header = FALSE)
  A <- as.matrix(A)
  base_forecasts <- .gen_poisson("dataForTests/Monthly-Poisson_basef.csv", seed = 42)
  res.buis_samples <- reconc_BUIS(A, base_forecasts, in_type = "samples", distr = "discrete", seed = 42)
  # Run IS Reconc
  base_forecasts_in <- read.csv(file = "dataForTests/Monthly-Poisson_basef.csv", header = FALSE)
  base_forecasts <- list()
  for (i in 1:nrow(base_forecasts_in)) {
    base_forecasts[[i]] <- list(lambda = base_forecasts_in[i, 1])
  }
  res.buis <- reconc_BUIS(A, base_forecasts, in_type = "params", distr = "poisson", num_samples = 100000, seed = 42)
  m <- mean(rowMeans(res.buis$bottom_rec_samples) - rowMeans(res.buis_samples$bottom_rec_samples))
  expect_equal(abs(m) < 1.5e-2, TRUE)
})

test_that("Monthly simple, in_type=='params', distr='nbinom'", {
  # Read samples from dataForTests (reproducibility)
  vals <- read.csv(file = "dataForTests/Monthly-Count_ts.csv", header = FALSE)

  # Create a count time series with monthly observations for 10 years
  y <- ts(data = vals, frequency = 12)

  # Create the aggregated yearly time series
  y_agg <- temporal_aggregation(y, agg_levels = c(1, 12))

  # We use a marginal forecast that computes for each month
  # the empirical mean and variance
  # the forecast is a negative binomial with those params
  fc_bottom <- list()
  for (i in seq(12)) {
    mm <- mean(y_agg$`f=12`[seq(i, 120, 12)])
    vv <- max(var(y_agg$`f=12`[seq(i, 120, 12)]), mm + 0.5)
    # cat("i: ",i, "mean: ",mm, "var: ",vv, "size: ",mm^2/(vv-mm), "prob: ",mm/vv, "\n")

    fc_bottom[[i]] <- list(size = mm^2 / (vv - mm), mu = mm)
  }

  # We compute the empirical mean and variance of the yearly ts
  # we forecast with a negative binomial with those parameters
  mm <- mean(y_agg$`f=1`)
  vv <- var(y_agg$`f=1`)
  fc_upper <- list(size = mm^2 / (vv - mm), prob = mm / vv)

  # Obtain the aggregation matrix for this hierarchy
  rec_mat <- get_reconc_matrices(c(1, 12), 12)

  base_forecasts <- append(list(fc_upper), fc_bottom)
  res.buis_params <- reconc_BUIS(rec_mat$A, base_forecasts, in_type = "params", distr = "nbinom", seed = 42)


  fc_upper_gauss <- list(mean = mm, cov = matrix(vv))
  res.mixCond <- reconc_MixCond(rec_mat$A, fc_bottom, fc_upper_gauss, bottom_in_type = "params", distr = "nbinom")
  upp_pmf <- PMF_from_samples(as.integer(res.buis_params$upper_rec_samples))

  expect_equal(res.mixCond$upper_rec$pmf[[1]], upp_pmf, tolerance = 0.1)
})


# reconc_gaussian with residuals -----------------------------------------------

# Shared fixture
set.seed(42)
.n_gauss   <- 3L
.L_gauss   <- 60L
.A_gauss   <- matrix(c(1, 1), nrow = 1)
.res_gauss <- matrix(rnorm(.L_gauss * .n_gauss), nrow = .L_gauss, ncol = .n_gauss)
.mu_gauss  <- c(11, 4, 6)  # incoherent: 11 != 4 + 6 = 10

test_that("reconc_gaussian runs without error when residuals are provided", {
  expect_no_error(
    reconc_gaussian(.A_gauss, base_fc_mean = .mu_gauss, residuals = .res_gauss)
  )
})

test_that("reconc_gaussian output structure is correct with residuals", {
  result <- reconc_gaussian(.A_gauss, base_fc_mean = .mu_gauss, residuals = .res_gauss)

  expect_named(result, c("bottom_rec_mean", "bottom_rec_covariance"), ignore.order = TRUE)
  expect_length(result$bottom_rec_mean, ncol(.A_gauss))
  expect_equal(dim(result$bottom_rec_covariance), c(ncol(.A_gauss), ncol(.A_gauss)))
})

test_that("reconc_gaussian reconciled covariance is symmetric and positive definite with residuals", {
  result <- reconc_gaussian(.A_gauss, base_fc_mean = .mu_gauss, residuals = .res_gauss)

  S <- result$bottom_rec_covariance
  expect_equal(S, t(S), tolerance = 1e-12)
  eigs <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that("reconc_gaussian satisfies hierarchical constraint with residuals", {
  result <- reconc_gaussian(.A_gauss, base_fc_mean = .mu_gauss, residuals = .res_gauss,
    return_upper = TRUE
  )

  expect_equal(
    as.vector(.A_gauss %*% result$bottom_rec_mean),
    as.vector(result$upper_rec_mean),
    tolerance = 1e-10
  )
  expect_equal(
    .A_gauss %*% result$bottom_rec_covariance %*% t(.A_gauss),
    result$upper_rec_covariance,
    tolerance = 1e-10
  )
})

test_that("return_upper=TRUE returns upper fields with residuals", {
  result <- reconc_gaussian(.A_gauss, base_fc_mean = .mu_gauss, residuals = .res_gauss,
    return_upper = TRUE
  )
  expect_true(all(c("upper_rec_mean", "upper_rec_covariance") %in% names(result)))
})

# Regression test: hard-coded reference values
# Reference computed with:
#   set.seed(42)
#   n_series <- 3L; L_train <- 60L
#   A <- matrix(c(1,1), nrow=1)
#   res_mat <- matrix(rnorm(L_train * n_series), nrow=L_train, ncol=n_series)
#   mu <- c(11, 4, 6)
#   reconc_gaussian(A, base_fc_mean=mu, residuals=res_mat, return_upper=TRUE)

test_that("reconc_gaussian with residuals matches hard-coded reference values", {
  set.seed(42)
  n_series <- 3L
  L_train  <- 60L
  A        <- matrix(c(1, 1), nrow = 1)
  res_mat  <- matrix(rnorm(L_train * n_series), nrow = L_train, ncol = n_series)
  mu       <- c(11, 4, 6)

  result <- reconc_gaussian(A, base_fc_mean = mu, residuals = res_mat, return_upper = TRUE)

  # bottom_rec_mean
  expect_equal(result$bottom_rec_mean, c(4.284368, 6.268177), tolerance = 1e-5)

  # bottom_rec_covariance
  expect_equal(
    result$bottom_rec_covariance,
    matrix(c(0.5930640, -0.2222453, -0.2222453, 0.5719506), nrow = 2),
    tolerance = 1e-5
  )

  # upper_rec_mean
  expect_equal(result$upper_rec_mean, 10.55254, tolerance = 1e-4)

  # upper_rec_covariance
  expect_equal(result$upper_rec_covariance, matrix(0.7205241), tolerance = 1e-5)
})

##############################################################################
