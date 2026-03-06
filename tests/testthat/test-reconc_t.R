# Shared fixtures --------------------------------------------------------------

# Simple 1-upper, 2-bottom hierarchy: u = b1 + b2
.A_simple <- matrix(c(1, 1), nrow = 1)

# Coherent base forecasts: upper = sum of bottoms
.mu_coh <- c(10, 4, 6) # u=10, b1=4, b2=6  =>  delta = A%*%b - u = 0

# Incoherent base forecasts: upper != sum of bottoms
.mu_inc <- c(11, 4, 6) # u=11, b1=4, b2=6  =>  delta = -1

# Posterior parameters
.nu_post <- 50
.Psi_post <- diag(c(3, 1, 1.5)) * .nu_post  # scale ~ nu * sigma^2 per series


# 1. Output structure ----------------------------------------------------------

test_that("reconc_t output contains expected fields (bottom only)", {
  res <- reconc_t(.A_simple, .mu_coh,
    posterior = list(nu = .nu_post, Psi = .Psi_post)
  )

  expect_named(res, c("bottom_rec_mean", "bottom_rec_scale_matrix", "bottom_rec_df"),
    ignore.order = TRUE
  )
  expect_length(res$bottom_rec_mean, ncol(.A_simple))
  expect_equal(dim(res$bottom_rec_scale_matrix), c(ncol(.A_simple), ncol(.A_simple)))
  expect_length(res$bottom_rec_df, 1)
})


# 2. Degrees of freedom formula ------------------------------------------------

test_that("bottom_rec_df equals nu_post - n_bottom + 1", {
  res <- reconc_t(.A_simple, .mu_coh,
    posterior = list(nu = .nu_post, Psi = .Psi_post)
  )
  expect_equal(res$bottom_rec_df, .nu_post - ncol(.A_simple) + 1)
})


# 3. Coherent forecasts: reconciled mean unchanged -----------------------------

test_that("reconciled mean is unchanged when base forecasts are already coherent", {
  res <- reconc_t(.A_simple, .mu_coh,
    posterior = list(nu = .nu_post, Psi = .Psi_post)
  )
  expect_equal(res$bottom_rec_mean, .mu_coh[(nrow(.A_simple) + 1):length(.mu_coh)],
    tolerance = 1e-10
  )
})


# 4. return_upper: upper_rec_mean == A %*% bottom_rec_mean ----------------------

test_that("return_upper=TRUE returns correct upper fields consistent with bottom", {
  res <- reconc_t(.A_simple, .mu_inc,
    posterior = list(nu = .nu_post, Psi = .Psi_post),
    return_upper = TRUE
  )

  expect_true(all(c("upper_rec_mean", "upper_rec_scale_matrix", "upper_rec_df") %in% names(res)))

  # upper_rec_mean must equal A %*% bottom_rec_mean (closure property)
  expect_equal(
    as.vector(res$upper_rec_mean),
    as.vector(.A_simple %*% res$bottom_rec_mean),
    tolerance = 1e-10
  )

  # upper_rec_scale_matrix must equal A %*% bottom_rec_scale_matrix %*% t(A)
  expect_equal(
    res$upper_rec_scale_matrix,
    .A_simple %*% res$bottom_rec_scale_matrix %*% t(.A_simple),
    tolerance = 1e-10
  )

  # degrees of freedom are the same for upper and bottom
  expect_equal(res$upper_rec_df, res$bottom_rec_df)
})

test_that("return_upper=FALSE (default) does not return upper fields", {
  res <- reconc_t(.A_simple, .mu_inc,
    posterior = list(nu = .nu_post, Psi = .Psi_post)
  )
  expect_false("upper_rec_mean" %in% names(res))
  expect_false("upper_rec_scale_matrix" %in% names(res))
})


# 5. Scale matrix is symmetric and positive-definite ---------------------------

test_that("bottom_rec_scale_matrix is symmetric and positive definite", {
  res <- reconc_t(.A_simple, .mu_inc,
    posterior = list(nu = .nu_post, Psi = .Psi_post)
  )
  S <- res$bottom_rec_scale_matrix

  # Symmetry
  expect_equal(S, t(S), tolerance = 1e-12)

  # Positive definiteness: all eigenvalues > 0
  eigs <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})


# 6. Incoherence: reconciled mean is coherent ----------------------------------
# After reconciliation, A %*% bottom_rec_mean always equals upper_rec_mean

test_that("reconciled forecasts satisfy the hierarchical constraint", {
  res <- reconc_t(.A_simple, .mu_inc,
    posterior = list(nu = .nu_post, Psi = .Psi_post),
    return_upper = TRUE
  )
  expect_equal(
    as.vector(.A_simple %*% res$bottom_rec_mean),
    as.vector(res$upper_rec_mean),
    tolerance = 1e-10
  )
})


# 7. Input via residuals (no error) --------------------------------------------

test_that("reconc_t runs without error when residuals and prior are provided", {
  set.seed(42)
  n <- 3
  L <- 30  # training length

  A <- .A_simple
  mu <- .mu_inc

  # Simulate residuals
  res_mat <- matrix(rnorm(L * n), nrow = L, ncol = n)

  # Simple diagonal prior
  nu_prior <- 10
  Psi_prior <- diag(n) * (nu_prior - n - 1)

  expect_no_error(
    reconc_t(A, mu, prior = list(nu = nu_prior, Psi = Psi_prior), residuals = res_mat)
  )
})


# 8. Posterior input reproduces prior+residuals path -------------------------

test_that("direct posterior input gives same result as prior+residuals path", {
  set.seed(7)
  n <- 3
  L <- 20

  A <- .A_simple
  mu <- .mu_inc
  res_mat <- matrix(rnorm(L * n), nrow = L, ncol = n)

  nu_prior <- 15
  Psi_prior <- diag(n) * (nu_prior - n - 1)

  # Path A: via prior + residuals
  res_via_prior <- reconc_t(A, mu,
    prior = list(nu = nu_prior, Psi = Psi_prior),
    residuals = res_mat
  )

  # Manual posterior computation (mirrors internal logic)
  Samp_cov <- crossprod(res_mat) / L
  lam <- .L_SHRINK_RECONC_T
  Samp_cov_shr <- (1 - lam) * Samp_cov + lam * diag(diag(Samp_cov))
  nu_post_man <- nu_prior + L
  Psi_post_man <- Psi_prior + L * Samp_cov_shr

  # Path B: directly provide posterior
  res_via_post <- reconc_t(A, mu,
    posterior = list(nu = nu_post_man, Psi = Psi_post_man)
  )

  expect_equal(res_via_prior$bottom_rec_mean, res_via_post$bottom_rec_mean, tolerance = 1e-10)
  expect_equal(res_via_prior$bottom_rec_scale_matrix, res_via_post$bottom_rec_scale_matrix,
    tolerance = 1e-10
  )
  expect_equal(res_via_prior$bottom_rec_df, res_via_post$bottom_rec_df)
})


# 9. Unknown argument triggers a warning -------------------------------------

test_that("unknown arguments in ... trigger a warning", {
  expect_warning(
    reconc_t(.A_simple, .mu_coh,
      posterior = list(nu = .nu_post, Psi = .Psi_post),
      unknown_arg = 99
    ),
    regexp = "unknown_arg"
  )
})


# 10. Standard usage: y_train + residuals (CASE 2b) ---------------------------
# This path runs .compute_naive_cov + multi_log_score_optimization to estimate
# nu_prior and Psi_prior, then updates to the posterior and reconciles.

# Shared fixture for CASE 2b tests (mts y_train with frequency=1, no seasonality)
set.seed(123)
.n_series <- 3L   # 1 upper + 2 bottom
.L_train  <- 50L
.A_2b     <- matrix(c(1, 1), nrow = 1)
.y_mat    <- ts(matrix(rnorm(.L_train * .n_series, mean = 10), nrow = .L_train, ncol = .n_series),
                frequency = 4)  # mts object with frequency=4, but no actual seasonality
.res_mat  <- matrix(rnorm(.L_train * .n_series), nrow = .L_train, ncol = .n_series)
.mu_2b    <- c(11, 4, 6)  # incoherent: 11 != 4 + 6 = 10

test_that("reconc_t runs without error with y_train + residuals (standard usage)", {
  expect_no_error(
    reconc_t(.A_2b, .mu_2b, y_train = .y_mat, residuals = .res_mat)
  )
})

test_that("reconc_t output structure is correct with y_train + residuals", {
  result <- reconc_t(.A_2b, .mu_2b, y_train = .y_mat, residuals = .res_mat)

  expect_named(result, c("bottom_rec_mean", "bottom_rec_scale_matrix", "bottom_rec_df"),
    ignore.order = TRUE
  )
  expect_length(result$bottom_rec_mean, ncol(.A_2b))
  expect_equal(dim(result$bottom_rec_scale_matrix), c(ncol(.A_2b), ncol(.A_2b)))
  expect_length(result$bottom_rec_df, 1)
  expect_true(result$bottom_rec_df > 0)
})

test_that("reconciled forecasts satisfy hierarchical constraint with y_train + residuals", {
  result <- reconc_t(.A_2b, .mu_2b, y_train = .y_mat, residuals = .res_mat,
    return_upper = TRUE
  )

  expect_equal(
    as.vector(.A_2b %*% result$bottom_rec_mean),
    as.vector(result$upper_rec_mean),
    tolerance = 1e-10
  )
})

test_that("bottom_rec_df equals posterior_nu - n_bottom + 1 with y_train + residuals", {
  result <- reconc_t(.A_2b, .mu_2b, y_train = .y_mat, residuals = .res_mat,
    return_parameters = TRUE
  )

  expect_equal(result$bottom_rec_df, result$posterior_nu - ncol(.A_2b) + 1)
})

test_that("posterior_nu equals optimized nu_prior + L with y_train + residuals", {
  result <- reconc_t(.A_2b, .mu_2b, y_train = .y_mat, residuals = .res_mat,
    return_parameters = TRUE
  )

  # posterior_nu = nu_prior + L; since nu_prior >= n + 2, posterior_nu > L
  expect_true(result$posterior_nu > .L_train)
  # Also must be valid for the t-distribution (> n_series - 1)
  expect_true(result$posterior_nu > .n_series - 1)
})

test_that("prior and posterior parameters are returned with y_train + residuals (return_parameters=TRUE)", {
  result <- reconc_t(.A_2b, .mu_2b, y_train = .y_mat, residuals = .res_mat,
    return_parameters = TRUE
  )

  expect_true(all(c("prior_nu", "prior_Psi", "posterior_nu", "posterior_Psi") %in% names(result)))
  expect_true(is.numeric(result$prior_nu))
  expect_true(is.matrix(result$prior_Psi))
  expect_equal(dim(result$prior_Psi), c(.n_series, .n_series))
  expect_true(is.numeric(result$posterior_nu))
  expect_true(is.matrix(result$posterior_Psi))
  expect_equal(dim(result$posterior_Psi), c(.n_series, .n_series))
})

test_that("bottom_rec_scale_matrix is symmetric and positive definite with y_train + residuals", {
  result <- reconc_t(.A_2b, .mu_2b, y_train = .y_mat, residuals = .res_mat)

  S <- result$bottom_rec_scale_matrix
  expect_equal(S, t(S), tolerance = 1e-12)
  eigs <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that("return_upper=TRUE works correctly with y_train + residuals", {
  result <- reconc_t(.A_2b, .mu_2b, y_train = .y_mat, residuals = .res_mat,
    return_upper = TRUE
  )

  expect_true(all(c("upper_rec_mean", "upper_rec_scale_matrix", "upper_rec_df") %in% names(result)))
  expect_equal(
    result$upper_rec_scale_matrix,
    .A_2b %*% result$bottom_rec_scale_matrix %*% t(.A_2b),
    tolerance = 1e-10
  )
  expect_equal(result$upper_rec_df, result$bottom_rec_df)
})

test_that("reconc_t standard usage works with explicit freq argument", {
  set.seed(42)
  freq <- 4L
  L <- 40L
  y_seas <- matrix(rnorm(L * .n_series, mean = 10), nrow = L, ncol = .n_series)
  res_seas <- matrix(rnorm(L * .n_series), nrow = L, ncol = .n_series)

  expect_no_error(
    reconc_t(.A_2b, .mu_2b, y_train = y_seas, residuals = res_seas, freq = freq)
  )
})

test_that("reconc_t standard usage works when y_train is an mts object", {
  set.seed(42)
  freq <- 4L
  L <- 40L
  y_mts <- ts(matrix(rnorm(L * .n_series, mean = 10), nrow = L, ncol = .n_series),
    frequency = freq
  )
  res_seas <- matrix(rnorm(L * .n_series), nrow = L, ncol = .n_series)

  expect_no_error(
    reconc_t(.A_2b, .mu_2b, y_train = y_mts, residuals = res_seas)
  )
})

test_that("reconc_t standard usage gives deterministic results for fixed inputs", {
  result1 <- reconc_t(.A_2b, .mu_2b, y_train = .y_mat, residuals = .res_mat)
  result2 <- reconc_t(.A_2b, .mu_2b, y_train = .y_mat, residuals = .res_mat)

  expect_equal(result1$bottom_rec_mean, result2$bottom_rec_mean)
  expect_equal(result1$bottom_rec_scale_matrix, result2$bottom_rec_scale_matrix)
  expect_equal(result1$bottom_rec_df, result2$bottom_rec_df)
})


# 12. Regression test: hard-coded reference values ----------------------------
# Reference computed with:
#   set.seed(123)
#   n_series <- 3L; L_train <- 50L
#   A_2b <- matrix(c(1,1), nrow=1)
#   y_mat <- matrix(rnorm(L_train * n_series, mean=10), nrow=L_train, ncol=n_series)
#   res_mat <- matrix(rnorm(L_train * n_series), nrow=L_train, ncol=n_series)
#   mu_2b <- c(11, 4, 6)
#   reconc_t(A_2b, mu_2b, y_train=y_mat, residuals=res_mat, return_parameters=TRUE)

test_that("reconc_t standard usage matches hard-coded reference values", {
  set.seed(123)
  n_series <- 3L
  L_train  <- 50L
  A        <- matrix(c(1, 1), nrow = 1)
  y_mat    <- matrix(rnorm(L_train * n_series, mean = 10), nrow = L_train, ncol = n_series)
  res_mat  <- matrix(rnorm(L_train * n_series),            nrow = L_train, ncol = n_series)
  mu       <- c(11, 4, 6)

  result <- suppressWarnings(
    reconc_t(A, mu, y_train = y_mat, residuals = res_mat, return_parameters = TRUE)
  )

  # bottom_rec_mean
  expect_equal(result$bottom_rec_mean,
    c(4.314540, 6.391355),
    tolerance = 1e-5
  )

  # bottom_rec_scale_matrix
  expect_equal(result$bottom_rec_scale_matrix,
    matrix(c(0.6503992, -0.3002090, -0.3002090, 0.5975068), nrow = 2),
    tolerance = 1e-5
  )

  # bottom_rec_df
  expect_equal(result$bottom_rec_df, 54)

  # posterior_nu
  expect_equal(result$posterior_nu, 55)
})
