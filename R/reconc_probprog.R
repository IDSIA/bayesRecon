# # Run IS Reconc
# S = data.table::fread(file = "../tests/testthat/dataForTests/Monthly-Gaussian_S.csv", header = FALSE)
# S = as.matrix(S)
# base_forecasts_in = data.table::fread(file = "../tests/testthat/dataForTests/Monthly-Gaussian_basef.csv", header = FALSE)
# base_forecasts = list()
# for (i in 1:nrow(base_forecasts_in)) {
#   base_forecasts[[i]] = list(as.numeric(base_forecasts_in[i,]))[[1]]
# }
# res.is = bayesReco::reconc_IS(S, base_forecasts,
#                    in_type = "params", distr = "gaussian", num_samples = 100000)
# # Run Gauss Reconc
# res.gauss = bayesReco::reconc_gaussian(base_forecasts_in[[1]], diag(base_forecasts_in[[2]]^2), S)
#
# # Test
# colMeans(res.is$reconciled_samples)
# as.numeric(rbind(res.gauss$upper_reconciled_mean, res.gauss$bottom_reconciled_mean))
#
# res.gauss_sample = reconc_gaussian_sample(base_forecasts,
#                                           S,
#                                           cores = 1,
#                                           chains = 1,
#                                           iter = 2000)
# rowMeans(res.gauss_sample$reconciled_samples)

##############################################################################

# data structure to know which bottom to sum in order to obtain each aggr
.build_ragged_array <- function(A){
  start <- vector(length = nrow(A))
  stop  <- vector(length = nrow(A))

  #sequence of indexes of the bottom to be summed to create all aggregates
  flat_array <- vector(length = sum(A))
  counter <- 1
  for (r_idx in 1:nrow(A)) {
    start[r_idx] <- counter
    for (c_idx in 1:ncol(A)) {
      if (A[r_idx, c_idx] ==1){
        flat_array[counter] <- c_idx
        counter <- counter +1
      }
      stop[r_idx] <- counter -1
    }
  }
  return(list(flat_array=flat_array, start=start, stop = stop, sum_a=sum(A)))
}

#' @title Probabilistic reconciliation via STAN
#' @param base_forecasts A list containing the base_forecasts, see details.
#' @param S Summing matrix (n x n_bottom)
#' @param cores Number of cores to use when executing the chains in parallel. If `NULL` use all available cores minus 1. See also \link[rstan]{sampling}.
#' @param chains A positive integer specifying the number of Markov chains. See also \link[rstan]{sampling}.
#' @param iter A positive integer specifying the number of samples for each chain. This number contains also the warm up iterations, by default `iter`/2. See also \link[rstan]{sampling}.
#' @param seed Seed for randomness reproducibility
#'
#' @return A list containing the reconciled forecasts. The list has the following named elements:
#'
#' * `reconciled_samples`: a matrix (n x (`iter` x `chains`/2 ) ) containing the reconciled samples for all time series
#'
#' @export
reconc_gaussian_sample <- function(base_forecasts,
                                   S,
                                   cores = NULL,
                                   chains = 4,
                                   iter = 2000,
                                   seed = NULL) {

  if (is.null(seed)) {
    seed = sample.int(.Machine$integer.max, 1)
  }

  set.seed(seed)

  if (is.null(cores) || cores > (parallel::detectCores()-1)) {
    options(mc.cores = parallel::detectCores()-1)
  } else {
    options(mc.cores = cores)
  }
  #rstan::rstan_options(auto_write = TRUE)

  split_hierarchy.res = .split_hierarchy(S, base_forecasts)
  A = split_hierarchy.res$A
  upper_base_forecasts = split_hierarchy.res$upper
  bottom_base_forecasts = split_hierarchy.res$bottom

  ragged_array = .build_ragged_array(A)

  n_b = ncol(A)
  n_a <- nrow(A)

  bottom = matrix(unlist(base_forecasts[split_hierarchy.res$bottom_idxs]), nrow=2, byrow = FALSE)
  aggr = matrix(unlist(base_forecasts[split_hierarchy.res$upper_idxs]), nrow=2, byrow = FALSE)

  normal_dat = list(n_b = n_b,
                     n_a = n_a,
                     bottom_means_base = bottom[1,],
                     bottom_sigma_base = bottom[2,],
                     aggr_means_base   = aggr[1,],
                     aggr_sigma_base   = aggr[2,],
                     sum_a             = ragged_array$sum_a,
                     flat_array        = ragged_array$flat_array,
                     start             = ragged_array$start,
                     stop              = ragged_array$stop
  )

  fit = rstan::sampling(stanmodels$gaussianReconc,
                    data = normal_dat,
                    algorithm = "NUTS",
                    chains = chains,
                    iter = iter,
                    seed = seed,
                    verbose = FALSE)

  matrix_of_draws  <- t(as.matrix(fit, pars=head(names(fit), -1)))
  matrix_of_draws_bottom = matrix_of_draws[1:n_b,]
  matrix_of_draws_aggr = matrix_of_draws[(n_b+1):(n_b+n_a),]

  order = rbind(data.frame(i=split_hierarchy.res$upper_idxs, t="aggr"),
                data.frame(i=split_hierarchy.res$bottom_idxs, t="bottom"))

  matrix_of_draws_ordered = matrix(data=NA, nrow = nrow(matrix_of_draws), ncol=ncol(matrix_of_draws))
  count_aggr = 1
  count_bottom = 1
  for (i in 1:nrow(order)) {
    pointer = order$i[i]
    if (order$t[i] == "aggr") {
      matrix_of_draws_ordered[i,] = matrix_of_draws_aggr[count_aggr,]
      count_aggr = count_aggr + 1
    }
    if (order$t[i] == "bottom") {
      matrix_of_draws_ordered[i,] = matrix_of_draws_bottom[count_bottom,]
      count_bottom = count_bottom + 1
    }
  }

  out = list(reconciled_samples = matrix_of_draws_ordered)
  return(out)
  }
