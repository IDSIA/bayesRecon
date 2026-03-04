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


# 5. return_parameters ---------------------------------------------------------

test_that("return_parameters=TRUE returns posterior_nu, posterior_Psi, C", {
  res <- reconc_t(.A_simple, .mu_inc,
    posterior = list(nu = .nu_post, Psi = .Psi_post),
    return_parameters = TRUE
  )

  expect_true(all(c("posterior_nu", "posterior_Psi", "C") %in% names(res)))
  expect_equal(res$posterior_nu, .nu_post)
  expect_equal(res$posterior_Psi, .Psi_post)
  expect_true(res$C > 0)
})

test_that("C equals 1/nu_tilde when base forecasts are coherent (delta=0)", {
  res <- reconc_t(.A_simple, .mu_coh,
    posterior = list(nu = .nu_post, Psi = .Psi_post),
    return_parameters = TRUE
  )
  nu_tilde <- .nu_post - ncol(.A_simple) + 1
  expect_equal(res$C, 1 / nu_tilde, tolerance = 1e-10)
})


# 6. Scale matrix is symmetric and positive-definite ---------------------------

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


# 7. Incoherence: reconciled mean is coherent ----------------------------------
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


# 8. Input via residuals (no error) --------------------------------------------

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


# 9. Posterior input reproduces prior+residuals path -------------------------

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


# 10. Unknown argument triggers a warning -------------------------------------

test_that("unknown arguments in ... trigger a warning", {
  expect_warning(
    reconc_t(.A_simple, .mu_coh,
      posterior = list(nu = .nu_post, Psi = .Psi_post),
      unknown_arg = 99
    ),
    regexp = "unknown_arg"
  )
})
