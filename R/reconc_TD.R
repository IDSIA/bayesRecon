### TODO:
# optimize cond_biv_sampling
# fix TD_sampling
# fix reconc_TD

# Sample from the distribution p(B_1, B_2 | B_1 + B_2 = u)
# B_1 and B_2 are distributed as pmf1 and pmf2
# u can be a vector!
cond_biv_sampling = function(u, pmf1, pmf2) {
  
  # In this way then we iterate over the one with shorter support:
  sw = FALSE
  if (length(pmf1) > length(pmf2)) {
    pmf_ = pmf1
    pmf1 = pmf2
    pmf2 = pmf_
    sw = TRUE
  }
  
  len_u = length(u)
  len_supp1 = length(pmf1)
  supp1 = 0:(len_supp1-1)
  
  W1 = t(matrix(pmf1, ncol=len_u, nrow=len_supp1))
  supp2 = outer(u, supp1, `-`)
  supp2[supp2<0] = Inf  # trick to get NA when we access pmf2 outside the support
  W2 = pmf2[supp2+1]    # add 1 because support starts from 0
  W2[is.na(W2)] = 0     # set NA to zero
  W2 = matrix(W2, nrow = len_u)  # back to matrix shape
  W = W1 * W2
  W = W / rowSums(W)  # normalize
  cumW = W %*% upper.tri(diag(len_supp1), diag = TRUE)  # cumulative sum on each row
  
  unif = runif(len_u)
  idxs = rowSums(unif > cumW) + 1
  b1 = supp1[idxs]
  if (sw) b1 = u - b1   # if we have switched, switch back
  
  return(list(b1, u-b1))
}

# Given a vector u of the upper values
# and a list of pmf (tables) of the bottom distr,
# returns samples (dim: n_bottom x length(u))
# from the conditional distr of the bottom given the upper values
TD_sampling = function(u, L_tab, Ltoll=1e-16, Rtoll=1e-7, 
                       smoothing=FALSE, al_smooth=NULL, lap_smooth=FALSE) {
  
  l_l_tab = rev(bottom_up_tab(L_tab, toll = toll, Rtoll = Rtoll, return_all = TRUE, 
                              smoothing=smoothing, al_smooth=al_smooth, lap_smooth=lap_smooth))
  
  b_old = matrix(u, nrow = 1)
  for (l_tab in l_l_tab[2:length(l_l_tab)]) {
    L = length(l_tab)
    b_new = matrix(ncol = length(u), nrow = L)
    for (j in 1:(L%/%2)) {
      b = cond_biv_sampling(b_old[j,], l_tab[[2*j-1]], l_tab[[2*j]])
      b_new[2*j-1,] = b[[1]]
      b_new[2*j,]   = b[[2]]
    }
    if (L%%2 == 1) b_new[L,] = b_old[L%/%2 + 1,]
    b_old = b_new
  }
  
  return(b_new)
}

# Reconciliation top-down using (...)
reconc_TD <- function(A, fc_bottom, mu_u, Sigma_u, names_u, 
                      N_samples = 1e3, Ltoll=1e-16, Rtoll=1e-7,
                      smooth_bottom=TRUE, al_smooth=NULL, lap_smooth=FALSE, 
                      return_samples = FALSE, seed = NULL) {
  
  set.seed(seed)
  
  n_b = ncol(A)
  n_u = nrow(A)
  
  lowest_rows = .lowest_lev(A)
  A_u = get_Au(A, lowest_rows)
  
  # Reconcile the upper
  rec_gauss_u = rec_gauss(A_u, 
                          mu_u=mu_u[-lowest_rows], Sigma_u=Sigma_u[-lowest_rows, -lowest_rows], 
                          mu_b=mu_u[lowest_rows],  Sigma_b=Sigma_u[lowest_rows, lowest_rows], 
                          Sigma_ub=Sigma_u[-lowest_rows, lowest_rows])
  
  # Sample from reconciled MVN on the lowest level of the upper
  U = mvtnorm::rmvnorm(N_samples, mean=rec_gauss_u$mu_b_tilde, sigma=rec_gauss_u$Sigma_b_tilde)  # n_samp x n_lowest_upp
  U = round(U)  # round to integer
  U_js = asplit(U, MARGIN = 2) # split into list of column vectors
  
  # Prepare lists of upper samples and bottom pmf
  n_low_u = length(lowest_rows)
  L_tab = lapply(fc_bottom, "[[", "samples")  # list of bottom empirical pmf
  L_tab_js = list()   # list of the lists of bottom pmf relative to the low upp
  for (j in lowest_rows) {
    Aj = A[j,]
    L_tab_js = c(L_tab_js, list(L_tab[as.logical(Aj)]))
  }
  
  # Check that each multivariate sample of U is contained in the support of the bottom-up distr
  samp_ok = mapply(check_support, U_js, L_tab_js, 
                   MoreArgs = list(Ltoll=Ltoll, Rtoll=Rtoll,
                                   smoothing=smooth_bottom, al_smooth=al_smooth, lap_smooth=lap_smooth))
  samp_ok = rowSums(samp_ok) == n_low_u
  U_js = lapply(U_js, "[", samp_ok) # Only keep the "good" upper samples
  print(paste0("Only ", round(sum(samp_ok)/N_samples, 3)*100, "% of the samples ",
               "are in the support of the bottom-up distribution"))
  
  # Get bottom samples via the prob top-down
  B = list()
  for (j in 1:n_low_u) {
    print(j)
    B[[j]] = TD_sampling(U_js[[j]], L_tab_js[[j]], 
                         Ltoll=Ltoll, Rtoll=Rtoll, 
                         smoothing=smooth_bottom, al_smooth=al_smooth, lap_smooth=lap_smooth)
  }
  B = do.call("rbind", B)  # dim: n_bottom x n_samples
  
  U = A %*% B
  
  rec_bottom = list()
  for (i in 1:n_b) {
    name = names(fc_bottom)[[i]]
    rec_bottom[[name]] = table(B[i,])
  }
  rec_upper = list()
  for (j in 1:n_u) {
    name = names_u[[j]]
    rec_upper[[name]] = table(U[j,])
  }
  
  out = list(
    bottom = rec_bottom,
    upper = rec_upper,
    n_samp = sum(samp_ok)
  )
  if (return_samples) {
    out$B_samples = B
    out$U_samples = U
  }
  
  return(out)
}



































