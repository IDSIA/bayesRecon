# Input checks ################################################################
.DISTR_SET = c("gaussian", "poisson", "nbinom")
.DISTR_SET2 = c("continuous", "discrete")
.check_input <- function(S, base_forecasts, in_type, distr) {

  if (!(nrow(S) == length(base_forecasts))) {
    stop("Input error: nrow(S) != length(base_forecasts)")
  }

  if (!(in_type %in% c("params", "samples"))) {
    stop("Input error: in_type must be a 'samples' or 'params'")
  }

  if (is.character(distr) & length(distr)==1) { # if distr is a string...
    if (in_type=="params" & !(distr %in% .DISTR_SET)) {
      stop(paste("Input error: if in_type='params', distr must be {", paste(.DISTR_SET, collapse = ', '), "}"))
    }
    if (in_type=="samples" & !(distr %in% .DISTR_SET2)) {
      stop(paste("Input error: if in_type='samples', distr must be {", paste(.DISTR_SET2, collapse = ', '), "}"))
    }
  }

  if (is.list(distr)) {
    if (!(nrow(S)==length(distr))) {
      stop("Input error: nrow(S) != length(distr)")
    }
  }

  # TODO if in_type=='samples' check sample sizes
  # TODO if in_type=='params' check that:
  # - gaussian: 2 parameters, (mu, sd)
  # - poisson:  1 parameter,  (lambda)
  # - nbinom:   2 parameters, (n, p)
  # TODO if distr is a list, check that entries are coherent
}

# Split bottoms, uppers #######################################################
.split_hierarchy <- function(S, Y) {
  bottom_idxs = which(rowSums(S)==1)
  upper_idxs = setdiff(1:nrow(S), bottom_idxs)
  A = S[upper_idxs,]
  upper = Y[upper_idxs]
  bottom = Y[bottom_idxs]
  out = list(A=A, upper=upper, bottom=bottom, upper_idxs=upper_idxs, bottom_idxs=bottom_idxs)
}

# Reconc utils ################################################################
.distr_sample <- function(params, distr_, n) {
  switch(distr_,
         "gaussian" = {samples = stats::rnorm(n=n, mean=params[[1]], sd=params[[2]])},
         "poisson"  = {samples = stats::rpois(n=n, lambda=params[[1]])},
         "negbin"   = {samples = stats::rnbinom(n=n, size=params[[1]], prob=params[[2]])},
  )
  return(samples)
}

.emp_pmf <- function(l, density_samples) {
  empirical_pmf = sapply(0:max(density_samples), function(i) sum(density_samples == i)/length(density_samples))
  w = sapply(l, function(i) empirical_pmf[i+1])
  return( w )
}

.fix_weights <- function(w) {
  w[is.na(w)] = 0
  if (sum(w)==0) {
    w = w + 1
  }
  return( w )
}

.compute_weights <- function(b, u, in_type_, distr_) {
  if (in_type_ == "samples") {
    if (distr_ == "discrete") {
      # Discrete samples
      w = .emp_pmf(b, u)
    } else if (distr_ == "continuous") {
      # KDE
      d = stats::density(u, bw="SJ", n=2**16)
      df = stats::approxfun(d)
      w = df(b)
    }
  } else if (in_type_ == "params") {
    switch(distr_,
           "gaussian" = {w = stats::dpois(x=b, mean=u[[1]], sd=u[[2]])},
           "poisson"  = {w = stats::dpois(x=b, lambda=u[[1]])},
           "negbin"   = {w = stats::dnbinom(x=b, size=u[[1]], prob=u[[2]])})
  }
  w = .fix_weights(w)
  return(w)
}

.resample <- function(S_, weights, num_samples=NA) {
  if (is.na(num_samples)) {
    num_samples = length(weights)
  }
  return( S_[sample(x=1:num_samples, num_samples, replace=TRUE, prob=weights),] )
}

# Misc ########################################################################
.shape <- function(m) {
  print(paste0("(",nrow(m),",",ncol(m),")"))
}
