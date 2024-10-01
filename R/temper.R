# TD conditioning using tempering
# Upper distr: pi_U * pi_bu^(1/T)
# temp = 1    --> conditioning
# temp = Inf  --> TD-cond
reconc_MixTemp = function(fc_bottom, fc_upper, temp,
                         bottom_in_type = "pmf", distr = NULL,
                         num_samples = 2e4, return_type = "pmf", 
                         suppress_warnings = FALSE, seed = NULL,
                         S = NULL, A = NULL, lowest_rows = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # TODO: modify .check_input_TD (input: A instead of S)
  # Check inputs 
  # .check_input_TD(S, fc_bottom, fc_upper, 
  #                 bottom_in_type, distr,
  #                 return_type)
  
  # Get aggr. matrix A and find the "lowest upper" 
  if (is.null(A)) A = .get_A_from_S(S)$A
  if (is.null(lowest_rows)) lowest_rows = .lowest_lev(A)

  n_u = nrow(A)
  n_b = ncol(A)
  n_u_low = length(lowest_rows)  # number of lowest upper
  
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
  
  # Get mean and covariance matrix of the MVN upper base forecasts
  mu_u    = fc_upper$mu
  Sigma_u = as.matrix(fc_upper$Sigma)
  
  # Get parameters of the MVN on the lowest upper level
  if (n_u == n_u_low) {  
    # If all the upper are lowest-upper, no reconciliation is needed  
    MVN_lower_params = list(
      mu    = mu_u,
      Sigma = Sigma_u
    )   
  } else {  
    # Else, analytically reconcile the upper
    
    # Get the aggregation matrix A_u and the summing matrix S_u for the upper sub-hierarchy
    A_u = .get_Au(A, lowest_rows)
    S_u = matrix(nrow = n_u, ncol = n_u_low)
    S_u[-lowest_rows,] = A_u
    S_u[lowest_rows,] = diag(n_u_low)
    
    # Analytically reconcile the upper
    rec_gauss_u = reconc_gaussian(S_u, mu_u, Sigma_u)
    MVN_lower_params = list(
      mu    = rec_gauss_u$bottom_reconciled_mean,
      Sigma = rec_gauss_u$bottom_reconciled_covariance
    )
  }
  
  if (is.infinite(temp)) {
    # If temp = Inf --> TD-cond 
    # Sample from MVN on the lowest level of the upper (dim: num_samples x n_u_low)
    U = .MVN_sample(n_samples = num_samples,
                    mu    = MVN_lower_params$mu, 
                    Sigma = MVN_lower_params$Sigma)  
    U = round(U)  # round to integer
  } else {
    # sample from MVN * pi_bu^(1/temp)
    bu_pmf = lapply(L_pmf_js, PMF.bottom_up)
    bu_pmf = lapply(bu_pmf, PMF.tempering, temp = temp)
    U = temper_sample(MVN_lower_params, bu_pmf, num_samples)
  }
  
  U_js = asplit(U, MARGIN = 2)  # split into list of column vectors
  
  # Check that each multiv. sample of U is contained in the supp of the bottom-up distr
  samp_ok = mapply(PMF.check_support, U_js, L_pmf_js)
  samp_ok = rowSums(samp_ok) == n_u_low
  # Only keep the "good" upper samples, and throw a warning if some samples are discarded:
  U_js = lapply(U_js, "[", samp_ok) 
  if (sum(samp_ok) != num_samples & !suppress_warnings) {
    # We round down to the nearest decimal
    warning(paste0("Only ", floor(sum(samp_ok)/num_samples*1000)/10, "% of the upper samples ",
                   "are in the support of the bottom-up distribution; ",
                   "the others are discarded."))
  }
  
  # Get bottom samples via the prob top-down
  B = list()
  for (j in 1:n_u_low) {
    B[[j]] = .TD_sampling(U_js[[j]], L_pmf_js[[j]])
  }
  B = do.call("rbind", B)  # dim: n_bottom x num_samples
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






















