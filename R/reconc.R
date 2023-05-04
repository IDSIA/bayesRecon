#' Reconciliation
#'
#' This function bla bla bla...
#'
#' @param S Summing matrix
#' @param base_forecasts base_forecasts
#' @param in_type A string 'samples'/'params'
#' @param distr A string 'continuous'/'discrete' or 'gaussian'/'poisson'/'nbinom'
#' @param num_samples number of samples
#' @param seed seed for randomness reproducibility
#'
#' @return Reconciled forecasts
#' @export
reconc <- function(
    S,
    base_forecasts,
    in_type,
    distr,
    num_samples = 2e4,
    seed = 42) {

  set.seed(seed)

  # Ensure that data inputs are valid
  .check_input(S, base_forecasts, in_type, distr)
  if (!is.list(distr)) {
    distr = rep(list(distr), nrow(S))
  }

  # Split bottoms, uppers
  split_hierarchy.res = .split_hierarchy(S, base_forecasts)
  A = split_hierarchy.res$A
  upper_base_forecasts = split_hierarchy.res$upper
  bottom_base_forecasts = split_hierarchy.res$bottom

  # H, G
  get_HG.res = .get_HG(A, upper_base_forecasts, distr[split_hierarchy.res$upper_idxs])
  H = get_HG.res$H
  upper_base_forecasts_H = get_HG.res$Hv
  G = get_HG.res$G
  upper_base_forecasts_G = get_HG.res$Gv

  # Reconciliation using BUIS
  n_upper = nrow(A)
  n_bottom = ncol(A)
  # 1. Bottom samples
  B = list()
  for (bi in 1:n_bottom) {
    if (in_type == "samples") {
      B[[bi]] = unlist(bottom_base_forecasts[[bi]])
    } else if (in_type == "params") {
      B[[bi]] = .distr_sample(
        bottom_base_forecasts[[bi]],
        distr[split_hierarchy.res$bottom_idxs][[bi]],
        num_samples)
    }
  }
  B = do.call("cbind", B) # B is a matrix (num_samples x n_bottom)

  # Bottom-Up IS on the hierarchical part
  for (hi in 1:nrow(H)) {
    c = H[hi,]
    b_mask = (c != 0)
    weights = .compute_weights(
      b = (B %*% c), # (num_samples x 1)
      u = unlist(upper_base_forecasts_H[[hi]]),
      in_type_ = in_type,
      distr_ = get_HG.res$Hdistr[[hi]])
    B[,b_mask] = .resample(B[,b_mask], weights)
  }

  if (!is.null(G)) {
    # Plain IS on the additional constraints
    weights = matrix(1, nrow=nrow(B))
    for (gi in 1:nrow(G)) {
      c = G[gi,]
      weights = weights * .compute_weights(
        b = (B %*% c),
        u = unlist(upper_base_forecasts_G[[gi]]),
        in_type_ = in_type,
        distr_ = get_HG.res$Gdistr[[gi]])
    }
    B = .resample(B, weights)
  }

  U = B %*% t(A)
  Y_reconc = cbind(U, B)

  out = list(bottom_reconciled_samples=B,
             upper_reconciled_samples=U,
             reconciled_forecasts=Y_reconc)
  return(out)
}
