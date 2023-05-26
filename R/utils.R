# Input checks
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

# Misc
.shape <- function(m) {
  print(paste0("(", nrow(m), ",", ncol(m), ")"))
}

# Functions for tests
.gen_gaussian <- function(params_file, seed=NULL) {
  set.seed(seed)
  params = utils::read.csv(file = params_file, header = FALSE)
  out = list()
  for (i in 1:nrow(params)) {
    out[[i]] = stats::rnorm(n=1e6, mean=params[[1]][i], sd=params[[2]][i])
  }
  return(out)
}

.gen_poisson <- function(params_file, seed=NULL) {
  set.seed(seed)
  params = utils::read.csv(file = params_file, header = FALSE)
  out = list()
  for (i in 1:nrow(params)) {
    out[[i]] = stats::rpois(n=1e6, lambda=params[[1]][i])
  }
  return(out)
}
