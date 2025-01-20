
.compute_upper_fc = function(upper_train, 
                             upper_fc_model,
                             upper_fc_args,
                             parallel_run = F, n_cpu = NULL,
                             H = 1) {
  
  if (upper_fc_model == "ets") {
    model_func = function(x) forecast::ets(x, additive.only = T)
  } else if (upper_fc_model == "arima") {
    model_func = forecast::auto.arima
  } else {
    stop("The specified model for the upper forecasts is not implemented")
  }
  
  L = ncol(upper_train)
  n_u = nrow(upper_train)
  
  mus = matrix(nrow = n_u, ncol = H)
  sds = matrix(nrow = n_u, ncol = H)
  res = matrix(nrow = n_u, ncol = L)
  
  if (!parallel_run) {
    for (j in 1:n_u) {
      u = ts(upper_train[j,], frequency = upper_fc_args$freq)
      model = model_func(u)
      fc = forecast::forecast(model, h=H, level=90)
      mus[j,] = fc$mean
      sds[j,] = c(fc$upper - fc$mean) / qnorm(0.95)
      res[j,] = model$residuals
    }
    
  } else {
    if (is.null(n_cpu)) {
      n_cpu = detectCores() - 2  
    }
    
    plan(multisession, workers = n_cpu)
    
    l = foreach(j = 1:n_u, .options.future = list(seed = TRUE)) %dofuture% {
      u = ts(upper_train[j,], frequency = upper_fc_args$freq)
      model = model_func(u)
      fc = forecast::forecast(model, h=H, level=90)
      list(
        mu  = fc$mean,
        sd  = c(fc$upper - fc$mean) / qnorm(0.95),
        res = model$residuals
      )
    }
    
    mus = do.call(rbind, lapply(l, "[[", "mu"))
    sds = do.call(rbind, lapply(l, "[[", "sd"))
    res = do.call(rbind, lapply(l, "[[", "res"))
  }

  fc_u = list()
  corr = cov2cor(schaferStrimmer_cov(t(res))$shrink_cov)
  for (h in 1:H) {
    fc_u[[h]] = list(
      mu = mus[,h],
      Sigma = corr * outer(sds[,h], sds[,h])
    )
  }
  
  return(fc_u)
}

.get_upper_sample = function(mu_u, Sigma_u, 
                             lowest_rows, A_u, n_u, n_u_low, num_samples) {
  
  if (n_u == n_u_low) {     
    # If all the upper are lowest-upper, just use the base distribution parameters
    mu    = mu_u
    Sigma = Sigma_u
    
  } else {
    # Else, analytically reconcile the upper and then sample from the lowest-uppers
    
    n_u_upp = n_u - n_u_low
    
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
    
    mu    = rec_gauss_u$bottom_reconciled_mean
    Sigma = rec_gauss_u$bottom_reconciled_covariance
  }
  
  # Sample from reconciled MVN on the lowest level of the upper (dim: num_samples x n_u_low)
  U = .MVN_sample(n_samples = num_samples,
                  mu    = mu, 
                  Sigma = Sigma) 
  U = round(U)          # round to integer
  mode(U) = "integer"   # change type for memory reasons
  
  # TODO: fix upper threshold (e.g. qnorm(1-1e-5, ...)) and remove samples 
  # that fall above; also remove negative samples
  # TODO: return correct number of samples
  U[U<0] = 0  # temporary
  
  return(U)
}

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
  
  if (any(u<0)) {
    stop("There are negative samples!")
  }
  
  # If one of the dimensions of the bpmf is 1, no sampling is needed
  if (s1 == 0) {
    return(list(rep(0,length(u)), u))
  } else if (s2 == 0) {
    return(list(u, rep(0,length(u))))
  }
  
  b1 = rep(NA, length(u))  # initialize empty vector
  
  Sums = row(bpmf)-1 + col(bpmf)-1  # matrix with the values of the support of the sum
  for (u_uniq in unique(u)) {  # loop over different values of u
    
    if (u_uniq > s1+s2) {
      # If the upper sample is outside the support of the sum of the bottom distributions,
      # split it proportionally to the lengths of the supports
      # E.g. supp1 = {0,1,2,3}; supp2 = {0,1,2}; u = 10
      # --> b1 = 6, b2 = 4
      b1_prop = u_uniq * s1 / (s1 + s2)
      b1[u == u_uniq] = round(jitter(b1_prop, amount = 1e-9)) # add jitter so that .5 is split randomly
      
    } else if (u_uniq == s1+s2) {
      b1[u == u_uniq] = s1
      
    } else if (u_uniq == 0) {
      b1[u == u_uniq] = 0
      
    } else {
      
      mask = (Sums == u_uniq)
      p = bpmf[mask]
      
      if (any(p<0)) {
        warning("Some probabilities are lower than zero")
        p[p<0] = 0   # set neg to zero
      }
      if (sum(p) == 0) {
        l = length(p)
        p = rep(1/l, l)   # if p is zero, set uniform probabilities
        warning("Probabilities are all zero")
      }
      
      p = p / sum(p)
      
      # supp1 = row(bpmf)[mask] - 1
      # more efficient implementation:
      supp1 = min(u_uniq, s1) : max(u_uniq - s2, 0)
      
      u_posit = (u == u_uniq)
      b1[u_posit] = sample(supp1, size = sum(u_posit), replace = TRUE, prob = p)
      
    }
    
  }
  return(list(b1, u-b1))
}

# (...) TODO: add explanation
# TODO: test the function
.TD_emp = function(u, b_train, copula_family, L_max_copula = NULL, max_frac_NA = 0.05) {
  
  n_b = nrow(b_train)
  if (n_b == 1) {return(u)}
  
  l_l_idx = rev(.BU_idxs(n_b))
  
  b_old = matrix(u, nrow = 1)
  for (l_idx in l_l_idx) {
    L = length(l_idx)
    b_new = matrix(ncol = length(u), nrow = L)
    for (j in 1:(L%/%2)) {
      # fictitious aggregated series (+ imputation):
      b1 = .aggr_ts(b_train[l_idx[[2*j-1]],,drop=F], max_frac_NA = max_frac_NA)  
      b2 = .aggr_ts(b_train[l_idx[[2*j]],,drop=F],   max_frac_NA = max_frac_NA)
      u_  = b_old[j,]
      # Ensure that all the upper samples are covered by the distribution:
      min_supp1 = ceiling(max(u_) * max(b1, na.rm=T) / (max(b1, na.rm=T)+max(b2, na.rm=T)))
      min_supp2 = ceiling(max(u_) * max(b2, na.rm=T) / (max(b1, na.rm=T)+max(b2, na.rm=T)))
      bpmf = bPMF.from_samples(b1, b2,
                               copula_family = copula_family,
                               # TODO: add parameters
                               # min_supp1 = min_supp1, min_supp2 = min_supp2,  # no need to specify the min_supp anymore!
                               L_max_copula = L_max_copula
                               )
      b = .cond_joint_biv_sampling(u_, bpmf)
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
                   copula_family,
                   L_max_copula = NULL, max_frac_NA = 0.05,
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
    
    # TODO: fix upper threshold (e.g. qnorm(1-1e-5, ...)) and remove samples 
    # that fall above; also remove negative samples
    # TODO: return correct number of samples
    U[U<0] = 0  # temporary
    
    U_js = asplit(U, MARGIN = 2) # split into list of column vectors
  }

  # Probabilistic top-down 
  B = matrix(nrow = n_b, ncol = num_samples)
  for (j in 1:n_u_low) {
    print(paste0("Lowest upper n.", j))
    mask_j = as.logical(A[lowest_rows[j], ])  # mask for the position of the bottom referring to lowest upper j
    B[mask_j, ] = .TD_emp(U_js[[j]], bottom_train[mask_j,],
                          copula_family = copula_family, L_max_copula = L_max_copula)
    # TODO: pass parameters to .TD_emp
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



# Function similar to hier_TD, but for multiple steps ahead
# Much faster than running H times the function hier_TD
# There might be memory issues; set a low num_samples 
# fc_upper is a list of length H, each entry is a list with $mu and $Sigma 
hier_TD_Hstep = function(A, fc_upper, bottom_train, 
                         copula_family,
                         L_max_copula = NULL, max_frac_NA = 0.05,
                         num_samples = 1e3, return_type = "pmf", 
                         parallel_run = FALSE, n_cpu = NULL,
                         suppress_warnings = FALSE, seed = NULL) {
                   
  if (!is.null(seed)) set.seed(seed)
  
  ### Check input ###
  # TODO: check_input
  # The forecast horizon is the length of the input upper fc list
  H = length(fc_upper)
  
  
  ### Pre-process aggregation matrix ###
  # Find the "lowest upper" 
  n_u = nrow(A)
  n_b = ncol(A)
  lowest_rows = .lowest_lev(A)
  n_u_low = length(lowest_rows)  # number of lowest upper
  n_u_upp = n_u - n_u_low        # number of "upper upper" 
  if (n_u_upp > 0) {
    # Get the aggregation matrix A_u for the upper sub-hierarchy
    A_u = .get_Au(A, lowest_rows)
  }
  
  
  ### Get upper samples ###
  U = matrix(nrow = num_samples*H, ncol = n_u_low)
  for (h in 1:H) {
    pos_h = (1+num_samples*(h-1)):(num_samples*h)
    U[pos_h,] = .get_upper_sample(fc_upper[[h]]$mu, as.matrix(fc_upper[[h]]$Sigma), 
                                  lowest_rows, A_u, n_u, n_u_low, num_samples)
  }
  # Create a list of column vectors
  U_js = lapply(seq_len(n_u_low), function(i) U[,i])
  rm(U)   # for memory reasons
  
  
  ### Probabilistic top-down ###
  if (!parallel_run) {
    B = matrix(nrow = n_b, ncol = num_samples*H)
    for (j in 1:n_u_low) {
      # print(paste0("Lowest upper n.", j))
      mask_j = as.logical(A[lowest_rows[j], ])  # mask for the position of the bottom referring to lowest upper j
      B[mask_j, ] = .TD_emp(U_js[[j]], bottom_train[mask_j,],
                            copula_family = copula_family, L_max_copula = L_max_copula)
      # TODO: pass parameters to .TD_emp
    }
    
  } else {
    if (is.null(n_cpu)) {
      n_cpu = detectCores() - 2  
    }
    
    plan(multisession, workers = n_cpu)
    
    par_data = lapply(1:n_u_low, function(j) list(U = U_js[[j]],
                                                  mask = as.logical(A[lowest_rows[j], ]),
                                                  bottom_train = bottom_train[as.logical(A[lowest_rows[j], ]), ]))
    B = foreach(data_j = par_data, .combine = rbind, .options.future = list(seed = TRUE)) %dofuture% {
      .TD_emp(data_j$U, data_j$bottom_train,
              copula_family = copula_family, L_max_copula = L_max_copula)
    }
    
    # Recover correct order of rows in B
    order_bottom = unlist(sapply(lowest_rows, function(lr) which(as.logical(A[lr,]))))
    B[order_bottom,] = B
  }
  U = A %*% B              # dim: n_upper x num_samples

  # Output: a list of length H, each entry is a list with the results
  # Include the marginal pmfs and/or the samples (depending on "return" inputs)
  out = rep(list(list(bottom_reconciled=list(), upper_reconciled=list())), H)
  for (h in 1:H) {
    pos_h = (1+num_samples*(h-1)):(num_samples*h)
    if (return_type %in% c('pmf', 'all')) {
      upper_pmf  = lapply(1:n_u, function(i) PMF.from_samples(U[i,pos_h]))
      bottom_pmf = lapply(1:n_b, function(i) PMF.from_samples(B[i,pos_h]))
      out[[h]]$bottom_reconciled$pmf = bottom_pmf
      out[[h]]$upper_reconciled$pmf = upper_pmf
    }
    if (return_type %in% c('samples','all')) {
      out[[h]]$bottom_reconciled$samples = B[,pos_h]
      out[[h]]$upper_reconciled$samples = U[,pos_h]
    }
  }

  return(out)
}


# ...
hier_TD_e2e = function(A, 
                       bottom_train, 
                       upper_fc_model,
                       upper_fc_args,
                       H = 1,
                       copula_family = "gauss",
                       L_max_copula = NULL, 
                       max_frac_NA = 0.1,
                       num_samples = 1e3, 
                       return_type = "pmf", 
                       parallel_run = FALSE, 
                       n_cpu = NULL,
                       suppress_warnings = FALSE, 
                       seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  tic.clear()
  tic.clearlog()
  
  ### Check input ###
  # TODO: check_input
  n_u = nrow(A)
  n_b = ncol(A)
  
  ### Compute upper forecasts ###
  # First, compute aggregated series
  print("Computing aggregated series...")
  tic()
  # find first time for which at least (1-max_frac_NA) of the bottom are not NA
  first_ind = min(which(colSums(is.na(bottom_train)) < max_frac_NA*n_b))
  # TODO: da sistemare? Se ci sono NA in mezzo non li vede...
  # impute bottom and compute upper series
  upper_train = A %*% t(apply(bottom_train[,first_ind:ncol(bottom_train)], 1,
                              .impute_ts))
  toc(log = T, quiet = T)
  
  # Compute upper forecasts
  print("Computing upper forecasts...")
  tic()
  fc_upper = .compute_upper_fc(upper_train, 
                               upper_fc_model,
                               upper_fc_args,
                               parallel_run = parallel_run, 
                               n_cpu = n_cpu,
                               H) 
  toc(log = T, quiet = T)
  
  
  ### Pre-process aggregation matrix ###
  # Find the "lowest upper" 
  print("Preprocessing of the A matrix...")
  tic()
  lowest_rows = .lowest_lev(A)
  n_u_low = length(lowest_rows)  # number of lowest upper
  n_u_upp = n_u - n_u_low        # number of "upper upper" 
  if (n_u_upp > 0) {
    # Get the aggregation matrix A_u for the upper sub-hierarchy
    A_u = .get_Au(A, lowest_rows)
  }
  toc(log = T, quiet = T)
  
  
  ### Get upper samples ###
  print("Drawing reconciled upper samples...")
  tic()
  U_low = matrix(nrow = n_u_low, ncol = num_samples*H) # matrix of samples from the lowest upper
  for (h in 1:H) {
    pos_h = (1+num_samples*(h-1)):(num_samples*h)
    U_low[,pos_h] = .get_upper_sample(fc_upper[[h]]$mu, as.matrix(fc_upper[[h]]$Sigma), 
                                      lowest_rows, A_u, n_u, n_u_low, num_samples)
  }
  toc(log = T, quiet = T)
  
  ### Probabilistic top-down ###
  if (!parallel_run) {
    print("Top-down algorithm, sequentially")
    tic()
    B = matrix(nrow = n_b, ncol = num_samples*H)
    for (j in 1:n_u_low) {
      # print(paste0("Lowest upper n.", j))
      mask_j = as.logical(A[lowest_rows[j], ])  # mask for the position of the bottom referring to lowest upper j
      B[mask_j, ] = .TD_emp(U_low[j,], bottom_train[mask_j,],
                            copula_family = copula_family, L_max_copula = L_max_copula)
      # TODO: pass parameters to .TD_emp
    }
    toc(log = T, quiet = T)
    
  } else {
    if (is.null(n_cpu)) {
      n_cpu = detectCores() - 2  
    }
    print(paste0("Top-down algorithm, in parallel, using ", n_cpu, " cores..."))
    tic()
    
    plan(multisession, workers = n_cpu)
    
    par_data = lapply(1:n_u_low, function(j) list(U = U_low[j,],
                                                  mask = as.logical(A[lowest_rows[j], ]),
                                                  bottom_train = bottom_train[as.logical(A[lowest_rows[j], ]), ]))
    B = foreach(data_j = par_data, .combine = rbind, .options.future = list(seed = TRUE)) %dofuture% {
      .TD_emp(data_j$U, data_j$bottom_train,
              copula_family = copula_family, L_max_copula = L_max_copula)
    }
    
    # Recover correct order of rows in B
    order_bottom = unlist(sapply(lowest_rows, function(lr) which(as.logical(A[lr,]))))
    B[order_bottom,] = B
    # TODO: optimize reordering of the rows of B
    # 2 possibilities (1 is probably more efficient but also more complicated):
    # 1) in-place reordering of rows (see chatGPT)
    # 2) modify foreach to return list of rows; then use a for loop to assign one
    #    row at a time, canceling then from the list
    toc(log = T, quiet = T)
  }
  print("Computing upper samples")
  tic()
  U = rbind(A_u %*% U_low, 
            U_low)
  toc(log = T, quiet = T)
  
  # Output: a list of length H, each entry is a list with the results
  # Include the marginal pmfs and/or the samples (depending on "return" inputs)
  print("Preparing output")
  tic()
  out = rep(list(list(bottom_reconciled=list(), upper_reconciled=list())), H)
  for (h in 1:H) {
    pos_h = (1+num_samples*(h-1)):(num_samples*h)
    if (return_type %in% c('pmf', 'all')) {
      upper_pmf  = lapply(1:n_u, function(i) PMF.from_samples(U[i,pos_h]))
      bottom_pmf = lapply(1:n_b, function(i) PMF.from_samples(B[i,pos_h]))
      out[[h]]$bottom_reconciled$pmf = bottom_pmf
      out[[h]]$upper_reconciled$pmf = upper_pmf
    }
    if (return_type %in% c('samples','all')) {
      out[[h]]$bottom_reconciled$samples = B[,pos_h]
      out[[h]]$upper_reconciled$samples = U[,pos_h]
    }
  }
  toc(log = T, quiet = T)
  
  times = sapply(tic.log(format=FALSE), function(x) x$toc - x$tic)
  print(paste0("Time for computing aggregated series: ",      times[[1]]))
  print(paste0("Time for computing upper forecasts: ",        times[[2]]))
  print(paste0("Time for preprocessing the A matrix: ",       times[[3]]))
  print(paste0("Time for drawing reconciled upper samples: ", times[[4]]))
  print(paste0("Time for top-down algorithm: ",               times[[5]]))
  print(paste0("Time for computing upper samples: ",          times[[6]]))
  print(paste0("Time for preparing output: ",                 times[[7]]))
  
  return(out)
}