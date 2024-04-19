# # Sample from the distribution p(B_1, B_2 | B_1 + B_2 = u)
# # B_1 and B_2 are distributed as pmf1 and pmf2
# # u can be a vector!
# .cond_biv_sampling_old = function(u, pmf1, pmf2) {
#   
#   # In this way then we iterate over the one with shorter support:
#   sw = FALSE
#   if (length(pmf1) > length(pmf2)) {
#     pmf_ = pmf1
#     pmf1 = pmf2
#     pmf2 = pmf_
#     sw = TRUE
#   }
#   
#   len_u = length(u)
#   len_supp1 = length(pmf1)
#   supp1 = 0:(len_supp1-1)
#   
#   W1 = t(matrix(pmf1, ncol=len_u, nrow=len_supp1))
#   supp2 = outer(u, supp1, `-`)
#   supp2[supp2<0] = Inf  # trick to get NA when we access pmf2 outside the support
#   W2 = pmf2[supp2+1]    # add 1 because support starts from 0
#   W2[is.na(W2)] = 0     # set NA to zero
#   W2 = matrix(W2, nrow = len_u)  # back to matrix shape
#   W = W1 * W2
#   W = W / rowSums(W)  # normalize
#   cumW = W %*% upper.tri(diag(len_supp1), diag = TRUE)  # cumulative sum on each row
#   
#   unif = runif(len_u)
#   idxs = rowSums(unif > cumW) + 1
#   b1 = supp1[idxs]
#   if (sw) b1 = u - b1   # if we have switched, switch back
#   
#   return(list(b1, u-b1))
# }

# OPTIMIZED IMPLEMENTATION: LOOP ON DIFFERENT VALUES OF U
#
# Sample from the distribution p(B_1, B_2 | B_1 + B_2 = u),
# where B_1 and B_2 are distributed as pmf1 and pmf2.
# u is a vector
.cond_biv_sampling = function(u, pmf1, pmf2) {
  
  # In this way then we iterate over the one with shorter support:
  sw = FALSE
  if (length(pmf1) > length(pmf2)) {
    pmf_ = pmf1
    pmf1 = pmf2
    pmf2 = pmf_
    sw = TRUE
  }
  
  b1 = rep(NA, length(u))  # initialize empty vector
  
  for (u_uniq in unique(u)) {  # loop over different values of u
    
    len_supp1 = length(pmf1)
    supp1 = 0:(len_supp1-1)
    p1 = pmf1
    
    supp2 = u_uniq - supp1
    supp2[supp2<0] = Inf  # trick to get NA when we access pmf2 outside the support
    p2 = pmf2[supp2+1]    # +1 because support starts from 0, but vector indexing from 1
    p2[is.na(p2)] = 0     # set NA to zero
    
    p = p1 * p2
    p = p / sum(p)
    
    u_posit = (u == u_uniq)
    b1[u_posit] = sample(supp1, size = sum(u_posit), replace = TRUE, prob = p)
  }
  
  if (sw) b1 = u - b1   # if we have switched, switch back
  
  return(list(b1, u-b1))
}

# Given a vector u of the upper values and a list of the bottom distr pmfs,
# returns samples (dim: n_bottom x length(u)) from the conditional distr 
# of the bottom given the upper values
.TD_sampling = function(u, bott_pmf, toll=1e-16, Rtoll=1e-7, 
                       smoothing=TRUE, al_smooth=NULL, lap_smooth=FALSE) {
  
  l_l_pmf = rev(PMF.bottom_up(bott_pmf, toll = toll, Rtoll = Rtoll, return_all = TRUE, 
                              smoothing=smoothing, al_smooth=al_smooth, lap_smooth=lap_smooth))
  
  b_old = matrix(u, nrow = 1)
  for (l_pmf in l_l_pmf[2:length(l_l_pmf)]) {
    L = length(l_pmf)
    b_new = matrix(ncol = length(u), nrow = L)
    for (j in 1:(L%/%2)) {
      b = .cond_biv_sampling(b_old[j,], l_pmf[[2*j-1]], l_pmf[[2*j]])
      b_new[2*j-1,] = b[[1]]
      b_new[2*j,]   = b[[2]]
    }
    if (L%%2 == 1) b_new[L,] = b_old[L%/%2 + 1,]
    b_old = b_new
  }
  
  return(b_new)
}



# Reconciliation top-down using (...)
# 
reconc_TDcond = function(S, fc_bottom, fc_upper, 
                     bottom_in_type = "pmf", distr = NULL,
                     N_samples = 2e4, return_pmf = TRUE, return_samples = FALSE, 
                     ...,
                     suppress_warnings = FALSE, seed = NULL) {
  
  set.seed(seed)
  
  # Parameters for convolution
  # toll=1e-16
  # Rtoll=1e-7
  # smooth_bottom=TRUE
  # al_smooth=NULL
  # lap_smooth=FALSE 
  
  # After testing the convolution parameters:
  # remove dots, remove comment above, and set the "best parameters" as default in 
  # PMF.check_support and .TD_sampling
  
  # Check inputs
  .check_input_TD(S, fc_bottom, fc_upper, 
                  bottom_in_type, distr,
                  return_pmf, return_samples)
  
  # Get aggr. matrix A and find the "lowest upper" 
  A = .get_A_from_S(S)$A
  n_u = nrow(A)
  n_b = ncol(A)
  lowest_rows = .lowest_lev(A)
  n_u_low = length(lowest_rows)  # number of lowest upper
  
  # Get mean and covariance matrix of the MVN upper base forecasts
  mu_u    = fc_upper$mu
  Sigma_u = as.matrix(fc_upper$Sigma)
  
  ### Get upper samples
  if (n_u == n_u_low) {     
    # If all the upper are lowest-upper, just sample from the base distribution
    U = .MVN_sample(N_samples, mu_u, Sigma_u)   # (dim: N_samples x n_u_low)
    U = round(U)                 # round to integer
    U_js = asplit(U, MARGIN = 2) # split into list of column vectors
    
  } else {
    # Else, analytically reconcile the upper and then sample from the lowest-uppers
    
    # Get the aggregation matrix A_u and the summing matrix S_u for the upper sub-hierarchy
    A_u = .get_Au(A, lowest_rows)
    S_u = matrix(nrow = n_u, ncol = n_u_low)
    S_u[-lowest_rows,] = A_u
    S_u[lowest_rows,] = diag(n_u_low)
    
    # Analytically reconcile the upper
    rec_gauss_u = reconc_gaussian(S_u, mu_u, Sigma_u)
    
    # Sample from reconciled MVN on the lowest level of the upper (dim: N_samples x n_u_low)
    U = .MVN_sample(n_samples = N_samples,
                    mu    = rec_gauss_u$bottom_reconciled_mean, 
                    Sigma = rec_gauss_u$bottom_reconciled_covariance)  
    U = round(U)                 # round to integer
    U_js = asplit(U, MARGIN = 2) # split into list of column vectors
  }

  # Prepare list of bottom pmf
  if (bottom_in_type == "pmf") {
    L_pmf = fc_bottom
  } else if (bottom_in_type == "samples") {
    L_pmf = lapply(fc_bottom, PMF.from_samples)
  } else if (bottom_in_type == "params") {
    L_pmf = lapply(fc_bottom, PMF.from_params, distr = distr)
  }
  
  # Prepare list of lists of bottom pmf relative to each lowest upper
  L_pmf_js = list()   
  for (j in lowest_rows) {
    Aj = A[j,]
    L_pmf_js = c(L_pmf_js, list(L_pmf[as.logical(Aj)]))
  }
  
  # Check that each multiv. sample of U is contained in the supp of the bottom-up distr
  samp_ok = mapply(PMF.check_support, U_js, L_pmf_js, 
                   MoreArgs = list(...))
  samp_ok = rowSums(samp_ok) == n_u_low
  # Only keep the "good" upper samples, and throw a warning if some samples are discarded:
  U_js = lapply(U_js, "[", samp_ok) 
  if (sum(samp_ok) != N_samples & !suppress_warnings) {
    warning(paste0("Only ", round(sum(samp_ok)/N_samples, 3)*100, "% of the upper samples ",
                   "are in the support of the bottom-up distribution; ",
                   "the others are discarded."))
  }
  
  # Get bottom samples via the prob top-down
  B = list()
  for (j in 1:n_u_low) {
    B[[j]] = .TD_sampling(U_js[[j]], L_pmf_js[[j]], ...)
  }
  B = do.call("rbind", B)  # dim: n_bottom x N_samples
  U = A %*% B              # dim: n_upper x N_samples
  
  # Prepare output: include the marginal pmfs and/or the samples (depending on "return" inputs)
  out = list()
  if (return_pmf) {
    upper_pmf  = lapply(1:n_u, function(i) PMF.from_samples(U[i,]))
    bottom_pmf = lapply(1:n_b, function(i) PMF.from_samples(B[i,]))
    out$pmf = list(
      upper  = upper_pmf,
      bottom = bottom_pmf
    )
  }
  if (return_samples) {
    out$samples = list(
      upper = U,
      bottom = B
    )
  } 

  return(out)
}




























