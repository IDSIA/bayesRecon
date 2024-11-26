
# Compute list of lists of vectors of indices
# (...) TODO: add explanation
.BU_idxs = function(m) {
  
  l_l_v = list()
  
  old_l = as.list(1:m)
  L = length(old_l)
  while (L>1) {
    l_l_v = c(l_l_v, list(old_l))
    new_l = list()
    for (j in 1:(L%/%2)) {
      v = c(old_l[[2*j-1]], old_l[[2*j]])
      new_l = c(new_l, list(v))
    }
    if (L%%2 == 1) new_l = c(new_l, list(old_l[[L]]))
    old_l = new_l
    L = length(old_l)
  }
  
  return(l_l_v)
}

# Sample from the distribution p(B_1, B_2 | B_1 + B_2 = u),
# where B_1 and B_2 are jointly distributed as bpmf.
# u is a vector
# supp(B_1) = {0, ..., nrow(bpmf)-1}
# supp(B_2) = {0, ..., ncol(bpmf)-1}
.cond_joint_biv_sampling = function(u, bpmf) {
  
  s1 = nrow(bpmf) - 1
  s2 = ncol(bpmf) - 1
  
  if (any(u<0 | u>(s1+s2))) {
    stop("Samples are outside the support")
  }
  # TODO: manage these cases
  
  b1 = rep(NA, length(u))  # initialize empty vector
  
  Sums = row(bpmf)-1 + col(bpmf)-1  # matrix with the values of the support of the sum
  for (u_uniq in unique(u)) {  # loop over different values of u
    
    p = bpmf[Sums==u_uniq]  # probabilities
    p = p / sum(p)
    
    supp1 = row(bpmf)[Sums==u_uniq]
    
    u_posit = (u == u_uniq)
    b1[u_posit] = sample(supp1, size = sum(u_posit), replace = TRUE, prob = p)
  }
  return(list(b1, u-b1))
}

# (...) TODO: add explanation
# TODO: test the function
.TD_emp = function(u, b_train) {
  
  n_b = nrow(b_train)
  if (n_b == 1) {return(u)}
  
  l_l_idx = rev(.BU_idxs(n_b))
  
  b_old = matrix(u, nrow = 1)
  for (l_idx in l_l_idx) {
    L = length(l_idx)
    b_new = matrix(ncol = length(u), nrow = L)
    for (j in 1:(L%/%2)) {
      # TODO: bPMF.from_samples
      bpmf = bPMF.from_samples_kde(colSums(train[l_idx[[2*j-1]]]),
                                  colSums(train[l_idx[[2*j]]]))
      b = .cond_joint_biv_sampling(b_old[j,], bpmf)
      b_new[2*j-1,] = b[[1]]
      b_new[2*j,]   = b[[2]]
    }
    if (L%%2 == 1) b_new[L,] = b_old[L%/%2 + 1,]
    b_old = b_new
  }
  
  return(b_new)
}


# Given Gaussian upper base forecasts and historical observations of the bottom:
# -compute reconciled upper fc (Gaussian reconciliation)
# -compute bottom fc via probabilistic top-down, using empirical in-sample pmf
# fc_upper: list with mu and Sigma
# bottom_train: matrix n_b x T
#
hier_TD = function(A, fc_upper, bottom_train, 
                   num_samples = 2e4, return_type = "pmf", 
                   suppress_warnings = FALSE, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # TODO: check_input
  
  # Find the "lowest upper" 
  n_u = nrow(A)
  n_b = ncol(A)
  lowest_rows = .lowest_lev(A)
  n_u_low = length(lowest_rows)  # number of lowest upper
  n_u_upp = n_u - n_u_low        # number of "upper upper" 
  
  # Get mean and covariance matrix of the MVN upper base forecasts
  mu_u    = fc_upper$mu
  Sigma_u = as.matrix(fc_upper$Sigma)
  
  ### Get upper samples
  if (n_u == n_u_low) {     
    # If all the upper are lowest-upper, just sample from the base distribution
    U = .MVN_sample(num_samples, mu_u, Sigma_u)   # (dim: num_samples x n_u_low)
    U = round(U)                 # round to integer
    U_js = asplit(U, MARGIN = 2) # split into list of column vectors
    
  } else {
    # Else, analytically reconcile the upper and then sample from the lowest-uppers
    
    # Get the aggregation matrix A_u for the upper sub-hierarchy
    A_u = .get_Au(A, lowest_rows)
    
    # Analytically reconcile the upper
    # The entries of mu_u must be in the correct order, i.e. rows of A_u (upper), columns of A_u (bottom)
    mu_u_ord = c(mu_u[-lowest_rows],mu_u[lowest_rows])  
    # Same for Sigma_u
    Sigma_u_ord = matrix(nrow=n_u, ncol=n_u)
    Sigma_u_ord[1:n_u_upp, 1:n_u_upp]             = Sigma_u[-lowest_rows,-lowest_rows]
    Sigma_u_ord[1:n_u_upp, (n_u_upp+1):n_u]       = Sigma_u[-lowest_rows,lowest_rows]
    Sigma_u_ord[(n_u_upp+1):n_u, 1:n_u_upp]       = Sigma_u[lowest_rows,-lowest_rows]
    Sigma_u_ord[(n_u_upp+1):n_u, (n_u_upp+1):n_u] = Sigma_u[lowest_rows,lowest_rows]
    rec_gauss_u = reconc_gaussian(A_u, mu_u_ord, Sigma_u_ord)
    
    # Sample from reconciled MVN on the lowest level of the upper (dim: num_samples x n_u_low)
    U = .MVN_sample(n_samples = num_samples,
                    mu    = rec_gauss_u$bottom_reconciled_mean, 
                    Sigma = rec_gauss_u$bottom_reconciled_covariance)  
    U = round(U)                 # round to integer
    U_js = asplit(U, MARGIN = 2) # split into list of column vectors
  }
  
  # TODO: check that samples are in the support?
  num_samples_ok = num_samples
  
  # Probabilistic top-down using empirical pmf 
  B = matrix(nrow = n_b, ncol = num_samples_ok)
  for (j in 1:n_u_low) {
    mask_j = as.logical(A[lowest_rows[j], ])  # mask for the position of the bottom referring to lowest upper j
    B[mask_j, ] = .TD_emp(U_js[[j]], bottom_train[mask_j,])
  }
  U = A %*% B              # dim: n_upper x num_samples
  
  # Prepare output: include the marginal pmfs and/or the samples (depending on "return" inputs)
  out = list(bottom_reconciled=list(), upper_reconciled=list())
  if (return_type %in% c('pmf', 'all')) {
    upper_pmf  = lapply(1:n_u, function(i) PMF.from_samples(U[i,]))
    bottom_pmf = lapply(1:n_b, function(i) PMF.from_samples(B[i,]))
    out$bottom_reconciled$pmf = bottom_pmf
    out$upper_reconciled$pmf = upper_pmf
  }
  if (return_type %in% c('samples','all')) {
    out$bottom_reconciled$samples = B
    out$upper_reconciled$samples = U
  } 
  
  return(out)
}






