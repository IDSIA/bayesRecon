# Input checks ################################################################
check_input <- function(S, base_forecasts, in_type, distr) {

  if (!(nrow(S) == length(base_forecasts))) {
    stop("Input error: nrow(S) != length(base_forecasts)")
  }

  if (!(in_type %in% c("params", "samples"))) {
    stop("Input error: in_type must be a 'samples' or 'params'")
  }

  if (is.character(distr) & length(distr)==1) {
    if (!(distr %in% c("gaussian", "poisson", "nbinom"))) {
      stop("Input error: distr must be {'gaussian', 'poisson', 'nbinom'}")
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
}

# Split bottoms, uppers #######################################################
split_hierarchy <- function(S, Y) {
  bottom_idxs = which(rowSums(S)==1)
  upper_idxs = setdiff(1:nrow(S), bottom_idxs)
  A = S[upper_idxs,]
  upper = Y[upper_idxs]
  bottom = Y[bottom_idxs]
  out = list(A=A, upper=upper, bottom=bottom, upper_idxs=upper_idxs, bottom_idxs=bottom_idxs)
}

# Misc ########################################################################
shape <- function(m) {
  print(paste0("(",nrow(m),",",ncol(m),")"))
}
