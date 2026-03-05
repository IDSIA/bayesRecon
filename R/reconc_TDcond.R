# Sample from the distribution p(B_1, B_2 | B_1 + B_2 = u),
# where B_1 and B_2 are distributed as pmf1 and pmf2.
# u is a vector
.cond_biv_sampling <- function(u, pmf1, pmf2) {
  # In this way then we iterate over the one with shorter support:
  sw <- FALSE
  if (length(pmf1) > length(pmf2)) {
    pmf_ <- pmf1
    pmf1 <- pmf2
    pmf2 <- pmf_
    sw <- TRUE
  }

  b1 <- rep(NA, length(u)) # initialize empty vector

  for (u_uniq in unique(u)) { # loop over different values of u

    len_supp1 <- length(pmf1)
    supp1 <- 0:(len_supp1 - 1)
    p1 <- pmf1

    supp2 <- u_uniq - supp1
    supp2[supp2 < 0] <- Inf # trick to get NA when we access pmf2 outside the support
    p2 <- pmf2[supp2 + 1] # +1 because support starts from 0, but vector indexing from 1
    p2[is.na(p2)] <- 0 # set NA to zero

    p <- p1 * p2
    p <- p / sum(p)

    u_posit <- (u == u_uniq)
    b1[u_posit] <- sample(supp1, size = sum(u_posit), replace = TRUE, prob = p)
  }

  if (sw) b1 <- u - b1 # if we have switched, switch back

  return(list(b1, u - b1))
}

# Given a vector u of the upper values and a list of the bottom distr pmfs,
# returns samples (dim: n_bottom x length(u)) from the conditional distr
# of the bottom given the upper values
.TD_sampling <- function(u, bott_pmf,
                         toll = .TOLL, Rtoll = .RTOLL, smoothing = TRUE,
                         al_smooth = .ALPHA_SMOOTHING, lap_smooth = .LAP_SMOOTHING) {
  # If the bottom pmf list contains only 1 element,
  # then the TD samples are simply a copy of the upper samples.
  if (length(bott_pmf) == 1) {
    return(matrix(u, nrow = 1))
  }

  l_l_pmf <- rev(PMF_bottom_up(bott_pmf,
    toll = toll, Rtoll = Rtoll, return_all = TRUE,
    smoothing = smoothing, al_smooth = al_smooth, lap_smooth = lap_smooth
  ))

  b_old <- matrix(u, nrow = 1)
  for (l_pmf in l_l_pmf[2:length(l_l_pmf)]) {
    L <- length(l_pmf)
    b_new <- matrix(ncol = length(u), nrow = L)
    for (j in 1:(L %/% 2)) {
      b <- .cond_biv_sampling(b_old[j, ], l_pmf[[2 * j - 1]], l_pmf[[2 * j]])
      b_new[2 * j - 1, ] <- b[[1]]
      b_new[2 * j, ] <- b[[2]]
    }
    if (L %% 2 == 1) b_new[L, ] <- b_old[L %/% 2 + 1, ]
    b_old <- b_new
  }

  return(b_new)
}

#' @rdname reconc_mixed
#' @examples
#'
#' library(bayesRecon)
#'
#' # Consider a simple hierarchy with two bottom and one upper
#' A <- matrix(c(1, 1), nrow = 1)
#' # The bottom forecasts are Poisson with lambda=15
#' lambda <- 15
#' n_tot <- 60
#' base_fc_bottom <- list()
#' base_fc_bottom[[1]] <- apply(matrix(seq(0, n_tot)), MARGIN = 1,
#'                              FUN = \(x) dpois(x, lambda = lambda))
#' base_fc_bottom[[2]] <- apply(matrix(seq(0, n_tot)), MARGIN = 1,
#'                              FUN = \(x) dpois(x, lambda = lambda))
#'
#' # The upper forecast is a Normal with mean 40 and std 5
#' base_fc_upper <- list(mean = 40, cov = matrix(c(5^2)))
#'
#' # Reconcile with reconc_TDcond
#' res.TDcond <- reconc_TDcond(A, base_fc_bottom, base_fc_upper)
#'
#' # Note that the bottom distributions are shifted to the right
#' PMF_summary(res.TDcond$bottom_rec$pmf[[1]])
#' PMF_summary(base_fc_bottom[[1]])
#'
#' PMF_summary(res.TDcond$bottom_rec$pmf[[2]])
#' PMF_summary(base_fc_bottom[[2]])
#'
#' # The upper distribution remains similar
#' PMF_summary(res.TDcond$upper_rec$pmf[[1]])
#' PMF_get_var(res.TDcond$upper_rec$pmf[[1]])
#'
#' ## Example 2: reconciliation with unbalanced hierarchy
#' # We consider the example in Fig. 9 of Zambon et al. (2024).
#'
#' # The hierarchy has 5 bottoms and 3 uppers
#' A <- matrix(c(
#'   1, 1, 1, 1, 1,
#'   1, 1, 0, 0, 0,
#'   0, 0, 1, 1, 0
#' ), nrow = 3, byrow = TRUE)
#' # Note that the 5th bottom only appears in the highest level, this is an unbalanced hierarchy.
#' n_upper <- nrow(A)
#' n_bottom <- ncol(A)
#'
#' # The bottom forecasts are Poisson with lambda=15
#' lambda <- 15
#' n_tot <- 60
#' base_fc_bottom <- list()
#' for (i in seq(n_bottom)) {
#'   base_fc_bottom[[i]] <- apply(matrix(seq(0, n_tot)), MARGIN = 1,
#'                                FUN = \(x) dpois(x, lambda = lambda))
#' }
#'
#' # The upper forecasts are a multivariate Gaussian
#' mean <- c(75, 30, 30)
#' cov <- matrix(c(
#'   5^2, 5, 5,
#'   5, 10, 0,
#'   5, 0, 10
#' ), nrow = 3, byrow = TRUE)
#'
#' base_fc_upper <- list(mean = mean, cov = cov)
#' \dontrun{
#' # If we reconcile with reconc_TDcond it won't work (unbalanced hierarchy)
#' res.TDcond <- reconc_TDcond(A, base_fc_bottom, base_fc_upper)
#' }
#'
#' # We can balance the hierarchy by duplicating the node b5:
#' # i) consider the time series observations for b5 as the upper u4,
#' # ii) fit the multivariate ts model for u1, u2, u3, u4.
#'
#' # In this example we simply assume that the forecast for u1-u4 is
#' # Gaussian with the mean and variance of u4 given by the parameters in b5.
#' mean_b5 <- lambda
#' var_b5 <- lambda
#' mean <- c(75, 30, 30, mean_b5)
#' cov <- matrix(c(
#'   5^2, 5, 5, 5,
#'   5, 10, 0, 0,
#'   5, 0, 10, 0,
#'   5, 0, 0, var_b5
#' ), nrow = 4, byrow = TRUE)
#' base_fc_upper <- list(mean = mean, cov = cov)
#'
#' # We also need to update the aggregation matrix
#' A <- matrix(c(
#'   1, 1, 1, 1, 1,
#'   1, 1, 0, 0, 0,
#'   0, 0, 1, 1, 0,
#'   0, 0, 0, 0, 1
#' ), nrow = 4, byrow = TRUE)
#'
#' # We can now reconcile with TDcond
#' res.TDcond <- reconc_TDcond(A, base_fc_bottom, base_fc_upper)
#'
#' # Note that the reconciled distribution of b5 and u4 are identical,
#' # keep this in mind when using the results of your reconciliation!
#' max(abs(res.TDcond$bottom_rec$pmf[[5]] - res.TDcond$upper_rec$pmf[[4]]))
#'
#' @export
reconc_TDcond <- function(A, base_fc_bottom, base_fc_upper,
                          bottom_in_type = "pmf", distr = NULL,
                          num_samples = 2e4, return_type = "pmf",
                          return_upper = TRUE,
                          suppress_warnings = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Check inputs
  .check_input_TD(
    A, base_fc_bottom, base_fc_upper,
    bottom_in_type, distr,
    return_type
  )

  # Get mean and covariance matrix of the MVN upper base forecasts
  mean_upper <- base_fc_upper$mean
  cov_upper <- as.matrix(base_fc_upper$cov)

  # Prepare list of bottom pmf
  if (bottom_in_type == "pmf") {
    L_pmf <- base_fc_bottom
  } else if (bottom_in_type == "samples") {
    L_pmf <- lapply(base_fc_bottom, PMF_from_samples)
  } else if (bottom_in_type == "params") {
    L_pmf <- lapply(base_fc_bottom, PMF_from_params, distr = distr)
  }

  out <- .core_reconc_TDcond(
    A, mean_upper, cov_upper, L_pmf, num_samples,
    return_type = return_type, 
    return_upper = return_upper,
    suppress_warnings = suppress_warnings
  )

  return(out)
}

#' Core Reconciliation via Top-Down Conditioning for Mixed Hierarchies
#'
#' Internal function that performs the core reconciliation logic using top-down conditioning
#' (TD) for mixed hierarchies. First, upper forecasts are reconciled analytically via conditioning
#' (if necessary), then bottom distributions are updated through probabilistic top-down procedure
#' by conditioning on the reconciled upper values.
#'
#' @param A Matrix (n_upper x n_bottom) defining the hierarchy where upper = A %*% bottom.
#' @param mean_upper Vector of mean forecasts for upper level.
#' @param cov_upper Covariance matrix of upper level forecasts.
#' @param L_pmf List of PMF objects representing the bottom level base forecasts.
#' @param num_samples Number of samples to draw from the reconciled distribution.
#' @param return_type Character string specifying return format: 'pmf', 'samples', or 'all'.
#' @param return_upper Logical, whether to return reconciled upper forecasts (default TRUE).
#' @param suppress_warnings Logical, whether to suppress warnings about samples outside support (default TRUE).
#' @param min_fraction_samples_ok Numeric between 0 and 1, minimum fraction of reconciled 
#' upper samples that must lie in the support of the bottom-up distribution (default 0.5). 
#' If the fraction is below this threshold, the function returns an error.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `bottom_rec`: List with reconciled bottom forecasts (pmf and/or samples).
#'     \item `upper_rec`:  (only if `return_upper = TRUE`) List with reconciled upper forecasts (pmf and/or samples).
#'   }
#'
#' @details
#' The function internally:
#' \enumerate{
#'   \item Identifies the "lowest upper" nodes in the hierarchy.
#'   \item If all uppers are lowest-uppers, samples directly from the upper MVN. Otherwise,
#'         analytically reconciles the upper hierarchy and samples from the lowest level.
#'   \item Reconciles bottom distributions by conditioning on the sampled/reconciled upper values
#'         using the probabilistic top-down algorithm.
#'   \item Discards samples that fall outside the support of the bottom-up distribution.
#' }
#'
#' @keywords internal
#' @export
.core_reconc_TDcond <- function(A, mean_upper, cov_upper, L_pmf, num_samples,
                                return_type, return_upper = TRUE,
                                suppress_warnings = TRUE,
                                min_fraction_samples_ok = .MIN_FRACTION_SAMPLES_OK) {
  # Find the "lowest upper"
  n_u <- nrow(A)
  n_b <- ncol(A)
  lowest_rows <- .lowest_lev(A)
  n_u_low <- length(lowest_rows) # number of lowest upper
  n_u_upp <- n_u - n_u_low # number of "upper upper"

  ### Get upper samples
  if (n_u == n_u_low) {
    # If all the upper are lowest-upper, just sample from the base distribution
    U <- .MVN_sample(num_samples, mean_upper, cov_upper) # (dim: num_samples x n_u_low)
    U <- round(U) # round to integer
    mode(U) <- "integer" # convert to integer
    U_js <- asplit(U, MARGIN = 2) # split into list of column vectors
  } else {
    # Else, analytically reconcile the upper and then sample from the lowest-uppers

    # Get the aggregation matrix A_u for the upper sub-hierarchy
    A_u <- .get_Au(A, lowest_rows)

    # Analytically reconcile the upper
    # The entries of mean_upper must be in the correct order, i.e. rows of A_u (upper), columns of A_u (bottom)
    mean_upper_ord <- c(mean_upper[-lowest_rows], mean_upper[lowest_rows])
    # Same for cov_upper
    cov_upper_ord <- matrix(nrow = n_u, ncol = n_u)
    cov_upper_ord[1:n_u_upp, 1:n_u_upp] <- cov_upper[-lowest_rows, -lowest_rows]
    cov_upper_ord[1:n_u_upp, (n_u_upp + 1):n_u] <- cov_upper[-lowest_rows, lowest_rows]
    cov_upper_ord[(n_u_upp + 1):n_u, 1:n_u_upp] <- cov_upper[lowest_rows, -lowest_rows]
    cov_upper_ord[(n_u_upp + 1):n_u, (n_u_upp + 1):n_u] <- cov_upper[lowest_rows, lowest_rows]
    rec_gauss_u <- reconc_gaussian(A_u, mean_upper_ord, cov_upper_ord)

    # Sample from reconciled MVN on the lowest level of the upper (dim: num_samples x n_u_low)
    U <- .MVN_sample(
      n_samples = num_samples,
      mu = rec_gauss_u$bottom_rec_mean,
      Sigma = rec_gauss_u$bottom_rec_covariance
    )
    U <- round(U) # round
    mode(U) <- "integer" # convert to integer
    U_js <- asplit(U, MARGIN = 2) # split into list of column vectors
  }

  # Prepare list of lists of bottom pmf relative to each lowest upper
  L_pmf_js <- list()
  for (j in lowest_rows) {
    Aj <- A[j, ]
    L_pmf_js <- c(L_pmf_js, list(L_pmf[as.logical(Aj)]))
  }

  # Check that each multiv. sample of U is contained in the supp of the bottom-up distr
  samp_ok <- mapply(PMF_check_support, U_js, L_pmf_js)
  samp_ok <- rowSums(samp_ok) == n_u_low

  # Check if the fraction of samples that are in the support of the bottom-up distribution 
  # is above the threshold; if not, stop 
  if (mean(samp_ok) < min_fraction_samples_ok) {
    stop(paste0("The fraction of reconciled upper samples that lie in the support 
                 of the bottom-up distribution is",
                 round(mean(samp_ok * 100), 1), 
                "which is below the minimum threshold (",
                 round(min_fraction_samples_ok * 100, 1),  "). ",
                "Consider increasing the variance of the base forecasts."))
  }

  # Resample the upper samples that are in the support of the bottom-up distribution, if necessary.
  # The number of output samples is equal to the number of input samples, but some of the output samples 
  # are duplicates of the "good" input samples.
  if (sum(samp_ok) < num_samples) {
    res_idxs <- sample(which(samp_ok), size = num_samples, replace = TRUE)
    U_js <- lapply(U_js, "[", res_idxs) 
    if (!suppress_warnings) {
      warning(paste0(
        "Only ", floor(sum(samp_ok) / num_samples * 1000) / 10, "% of the upper samples ",
        "are in the support of the bottom-up distribution; ",
        "the others are discarded and the remaining ones are resampled with replacement."
      ))
    }
  }

  # Get bottom samples via the prob top-down
  B <- matrix(nrow = n_b, ncol = num_samples)
  for (j in 1:n_u_low) {
    mask_j <- as.logical(A[lowest_rows[j], ]) # mask for the position of the bottom referring to lowest upper j
    B[mask_j, ] <- .TD_sampling(U_js[[j]], L_pmf_js[[j]])
  }
  U <- A %*% B # dim: n_upper x num_samples

  # Prepare output: include the marginal pmfs and/or the samples (depending on "return" inputs)
  out <- list(bottom_rec = list(), upper_rec = list())
  if (return_type %in% c("pmf", "all")) {
    bottom_pmf <- lapply(1:n_b, function(i) PMF_from_samples(B[i, ]))
    out$bottom_rec$pmf <- bottom_pmf
    if (return_upper) {
      upper_pmf <- lapply(1:n_u, function(i) PMF_from_samples(U[i, ]))
      out$upper_rec$pmf <- upper_pmf
    }
  }
  if (return_type %in% c("samples", "all")) {
    out$bottom_rec$samples <- B
    if (return_upper) {
      out$upper_rec$samples <- U
    }
  }
  return(out)
}
