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

  if (is.character(distr) &
      length(distr) == 1) {
    # if distr is a string...
    if (in_type == "params" & !(distr %in% .DISTR_SET)) {
      stop(paste(
        "Input error: if in_type='params', distr must be {",
        paste(.DISTR_SET, collapse = ', '),
        "}"
      ))
    }
    if (in_type == "samples" & !(distr %in% .DISTR_SET2)) {
      stop(paste(
        "Input error: if in_type='samples', distr must be {",
        paste(.DISTR_SET2, collapse = ', '),
        "}"
      ))
    }
  }

  if (is.list(distr)) {
    if (!(nrow(S) == length(distr))) {
      stop("Input error: nrow(S) != length(distr)")
    }
  }

  # Eventually:
  # TODO if in_type=='samples' check sample sizes
  # TODO if in_type=='params' check that:
  # - gaussian: 2 parameters, (mu, sd)
  # - poisson:  1 parameter,  (lambda)
  # - nbinom:   2 parameters, (n, p)
  # TODO if distr is a list, check that entries are coherent
}


# Non-overlapping temporal aggregation of a time series #######################
#' Time series temporal aggregation
#'
#' This function bla bla bla...
#'
#' @param y Univariate time series of class ts.
#' @param aggf User-selected list of aggregates to consider.
#'
#' @return A list of aggregates time series.
#' @export
temporal_aggregation <- function(y, aggf=NULL) {
  f = stats::frequency(y)
  L = length(y)
  s = stats::time(y)[1]
  if (is.null(aggf)) {
    aggf = c()
    for (i in 1:f) {
      if (f %% i == 0 && L >= i) {
        aggf = c(aggf, i)
      }
    }
  } else {
    aggf = aggf[aggf <= L]
  }
  out = list()
  for (i in 1:length(aggf)) {
    k = aggf[i]
    num_aggs = floor(L / k)
    y_trunc = y[(L - num_aggs*k + 1):L]
    y_matrix = matrix(y_trunc, nrow = k, ncol = num_aggs)
    y_start = s + (L - num_aggs * k) / f
    y_f = f / k
    y_agg = stats::ts(data = apply(y_matrix, 2, sum), frequency = y_f, start = y_start)
    out[[i]] = y_agg
  }
  names(out) <- paste0("f=", f / aggf)
  out = rev(out)
  return(out)
}
# Get A from S ################################################################
.getAfromS <- function(S) {
  bottom_idxs = which(rowSums(S) == 1)
  upper_idxs = setdiff(1:nrow(S), bottom_idxs)
  A = S[upper_idxs, ]
  out = list(A = A,
             upper_idxs = upper_idxs,
             bottom_idxs = bottom_idxs)
  return(out)
}


# Split bottoms, uppers #######################################################
.split_hierarchy <- function(S, Y) {
  getAfromS.res = .getAfromS(S)
  upper = Y[getAfromS.res$upper_idxs]
  bottom = Y[getAfromS.res$bottom_idxs]
  out = list(
    A = getAfromS.res$A,
    upper = upper,
    bottom = bottom,
    upper_idxs = getAfromS.res$upper_idxs,
    bottom_idxs = getAfromS.res$bottom_idxs
  )
  return(out)
}


# Build A / S matrices ########################################################
#' A/S matrices building
#'
#' This function bla bla bla...
#'
#' @param aggf User-selected list of aggregates to consider.
#' @param bottom.f Integer seasonal period of the bottom time series.
#' @param bottom.H Bottom time series forecasting steps.
#'
#' @return A, S matrices
#' @export
reconc_matrices <- function(aggf, bottom.f, bottom.H) {
  A = list()
  for (i in 1:length(aggf)) {
    k = aggf[i]
    if (k==1) {
      next
    }
    k.r = bottom.H / k
    k.A = matrix(data = 0, nrow = k.r, ncol = bottom.H)
    coli = 1
    for (r in 1:k.r) {
      k.A[r,coli:(coli+k-1)] = 1
      coli = coli + k
    }
    A[[i]] = k.A

  }
  A = do.call(rbind, rev(A))
  S = rbind(A, diag(bottom.H))
  out = list(A=A, S=S)
  return(out)
}



# IS utils ####################################################################
.distr_sample <- function(params, distr_, n) {
  switch(
    distr_,
    "gaussian" = {
      samples = stats::rnorm(n = n, mean = params[[1]], sd = params[[2]])
    },
    "poisson"  = {
      samples = stats::rpois(n = n, lambda = params[[1]])
    },
    "negbin"   = {
      samples = stats::rnbinom(n = n,
                               size = params[[1]],
                               prob = params[[2]])
    },
  )
  return(samples)
}

.emp_pmf <- function(l, density_samples) {
  empirical_pmf = sapply(0:max(density_samples), function(i)
    sum(density_samples == i) / length(density_samples))
  w = sapply(l, function(i)
    empirical_pmf[i + 1])
  return(w)
}

.fix_weights <- function(w) {
  w[is.na(w)] = 0
  if (sum(w) == 0) {
    w = w + 1
  }
  return(w)
}

.compute_weights <- function(b, u, in_type_, distr_) {
  if (in_type_ == "samples") {
    if (distr_ == "discrete") {
      # Discrete samples
      w = .emp_pmf(b, u)
    } else if (distr_ == "continuous") {
      # KDE
      d = stats::density(u, bw = "SJ", n = 2 ** 16)
      df = stats::approxfun(d)
      w = df(b)
    }
  } else if (in_type_ == "params") {
    switch(
      distr_,
      "gaussian" = {
        w = stats::dnorm(x = b,
                         mean = u[[1]],
                         sd = u[[2]])
      },
      "poisson"  = {
        w = stats::dpois(x = b, lambda = u[[1]])
      },
      "negbin"   = {
        w = stats::dnbinom(x = b,
                           size = u[[1]],
                           prob = u[[2]])
      }
    )
  }
  w = .fix_weights(w)
  return(w)
}

.resample <- function(S_, weights, num_samples = NA) {
  if (is.na(num_samples)) {
    num_samples = length(weights)
  }
  return(S_[sample(
    x = 1:num_samples,
    num_samples,
    replace = TRUE,
    prob = weights
  ), ])
}



# Misc ########################################################################
.shape <- function(m) {
  print(paste0("(", nrow(m), ",", ncol(m), ")"))
}
