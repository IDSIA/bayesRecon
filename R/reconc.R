source("utils.R")
source("hierarchy.R")

distr_sample <- function(params, distr_, n) {
  switch(distr_,
         "gaussian" = {samples = stats::rnorm(n=n, mean=params[[1]], sd=params[[2]])},
         "poisson"  = {samples = stats::rpois(n=n, lambda=params[[1]])},
         "negbin"   = {samples = stats::rnbinom(n=n, size=params[[1]], prob=params[[2]])},
  )
  return(samples)
}

emp_pmf <- function(l, density_samples) {
  empirical_pmf = sapply(0:max(density_samples), function(i) sum(density_samples == i)/length(density_samples))
  w = sapply(l, function(i) empirical_pmf[i+1])
  return( w )
}

compute_weights <- function(b, u, in_type_, distr_) {
  if (in_type_ == "samples") {
    if (distr_ %in% c("poisson","negbin")) {
      # Discrete samples
      w = emp_pmf(b, u)
    } else if (distr_ == "gaussian") {
      # TODO w = kde(...)
    }
  } else if (in_type_ == "params") {
    switch(distr_,
           "gaussian" = {w = stats::dpois(x=b, mean=u[[1]], sd=u[[2]])},
           "poisson"  = {w = stats::dpois(x=b, lambda=u[[1]])},
           "negbin"   = {w = stats::dnbinom(x=b, size=u[[1]], prob=u[[2]], mu=u[[3]])}
           )
  }
  return(w)
}

resample <- function(S, weights, num_samples=NA) {
  if (is.na(num_samples)) {
    num_samples = length(weights)
  }
  return( S[sample(x=1:num_samples, num_samples, replace=TRUE, prob=weights),] )
}

reconc <- function(
    S,
    base_forecasts,
    in_type,
    distr,
    num_samples = 2e4) {

  # Ensure that data inputs are valid #########################################
  check_input(S, base_forecasts, in_type, distr)
  if (!is.list(distr)) {
    distr = rep(list(distr), nrow(S))
  }

  # Split bottoms, uppers #####################################################
  split_hierarchy.res = split_hierarchy(S, base_forecasts)
  A = split_hierarchy.res$A
  upper_base_forecasts = split_hierarchy.res$upper
  bottom_base_forecasts = split_hierarchy.res$bottom

  # H, G ######################################################################
  get_H.res = get_H(A, upper_base_forecasts, distr[split_hierarchy.res$upper_idxs])
  H = get_H.res$H
  upper_base_forecasts_H = get_H.res$Hv
  G = get_H.res$G
  upper_base_forecasts_G = get_H.res$Gv

  # Reconciliation using BUIS #################################################
  n_upper = nrow(A)
  n_bottom = ncol(A)
  # 1. Bottom samples
  B = list()
  for (bi in 1:n_bottom) {
    if (in_type == "samples") {
      B[[bi]] = unlist(bottom_base_forecasts[[bi]])
    } else if (in_type == "params") {
      B[[bi]] = distr_sample(
        bottom_base_forecasts[[bi]],
        distr[split_hierarchy.res$bottom_idxs][[bi]],
        num_samples)
    }
  }
  B = do.call("cbind", B)

  # Bottom-Up IS on the hierarchical part #####################################
  for (hi in 1:nrow(H)) {
    c = H[hi,]
    b_mask = (c != 0)
    weights = compute_weights(
      b = (B %*% c),
      u = unlist(upper_base_forecasts_H[[hi]]),
      in_type_ = in_type,
      distr_ = get_H.res$Hdistr[[hi]])
    B[,b_mask] = resample(B[,b_mask], weights)
  }

  if (!is.null(G)) {
    # Plain IS on the additional constraints ##################################
    # TODO
  }

  U = B %*% t(A)
  Y_reconc = cbind(U, B)

  out = list(bottom_reconciled_samples=B,
             upper_reconciled_samples=U,
             reconciled_forecasts=Y_reconc)
  return(out)
}


