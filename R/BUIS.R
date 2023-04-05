### Dependencies
# - stats

###############################################################################

### Utils
shape <- function(m) {
  print(paste0("(",nrow(m),",",ncol(m),")"))
}

resample <- function(S, weights, num_samples=NA) {
  if (is.na(num_samples)) {
    num_samples = length(weights)
  }
  return( S[sample(x=1:num_samples, num_samples, replace=TRUE, prob=weights),] )
}

emp_pmf <- function(l, density_samples) {
  empirical_pmf = sapply(0:max(density_samples), function(i) sum(density_samples == i)/length(density_samples))
  w = sapply(l, function(i) empirical_pmf[i+1])
  return( w )
}
###############################################################################

H = t(matrix(data=c(1,1,0,0,0,0,0,0,
                    0,0,1,1,0,0,0,0,
                    0,0,0,0,1,1,0,0,
                    0,0,0,0,0,0,1,1,
                    1,1,1,1,0,0,0,0,
                    0,0,0,0,1,1,1,1,
                    1,1,1,1,1,1,1,1), nrow=8, ncol=7))
G = NULL

### Reconciliation using BUIS
params_b = list(
  matrix(data=c(5.634849165190505, 9.833589192410015, 6.3023800293289165, 
                9.486182621822866, 6.883748580948357, 6.681108721672265, 
                7.256882352376998, 9.201275416306906))
)
params_u = list(
  matrix(data=c(20.108969864880674, 20.525131446497316, 17.63431449340681, 
                21.395605099289075, 40.63410131137799, 39.02991959269589, 
                79.66402090407388))
)
params_hier = list(
  matrix(data=c(20.108969864880674, 20.525131446497316, 17.63431449340681, 
                21.395605099289075, 40.63410131137799, 39.02991959269589, 
                79.66402090407388))
)
params_addconstr = NULL

recBUIS <- function(
    distr,
    params_b,
    params_hier,
    H,
    params_addconstr = NULL,
    G = NULL,
    num_samples = 1e5) {
  
  if ( !(distr %in% c("Poisson", "Gaussian")) ) {
    stop("Please specify the type of distribution: 'Gaussian' or 'Poisson'")
  }
  
  n_bottom = ncol(H)
  n_hier = nrow(H)
  
  #Sample from the base distribution of b
  B = matrix(data=NA, nrow=num_samples, ncol=n_bottom)
  for (i in 1:n_bottom) {
    if (distr == "Poisson") {
      B[,i] = stats::rpois(n=num_samples, lambda=params_b[[1]][i])
    }
    if (distr == "Gaussian") {
      # TODO B[,i] = ...
    }
  }
  
  #Bottom-Up IS on the hierarchical part
  for (i in 1:n_hier) {   # It's important that the rows of H are in the 
                          # correct order (i.e. bottom-up)!
    if (distr == "Poisson") {
      w_i = stats::dpois(x=(B %*% H[i,]), lambda = params_u[[1]][i])
    }
    if (distr == "Gaussian") {
      # TODO w_i = ...
    }
    bottom_mask = (H[i,] != 0)
    B[,bottom_mask] = resample(B[,bottom_mask], w_i)
  }
  
  if (!is.null(G)) {
    # Plain IS on the additional constraints
    # TODO
  }
  
  A = if (is.null(G)) H else rbind(H, G)
  u = B %*% t(A)
  y = cbind(u, B)
  
  out = list(B_samples=B, U_samples=u, Y_samples=y)
  return(out)
}

system.time({
  recBUIS.res = recBUIS(
    distr = "Poisson",
    params_b = params_b,
    params_hier = params_hier,
    H = H,
    params_addconstr = params_addconstr,
    G = G,
  )
})


### Reconciliation using BUIS with samples
path = "../LorenzoRepoBUIS/"
samples_b = reticulate::py_load_object(filename = paste0(path,"samples_b.pkl"))
samples_hier = reticulate::py_load_object(filename = paste0(path,"samples_hier.pkl"))

recBUISSamples <- function(
    samples_b, 
    samples_hier, 
    H, 
    samples_addconstr = NULL, 
    G = NULL,
    distr = NULL) {
  
  if ( !(distr %in% c("discrete", "continuous")) ) {
    stop("Please specify the type of distribution: 'discrete' or 'continuous'")
  }
  
  n_bottom = ncol(H)
  n_hier = nrow(H)
  
  B = samples_b
  
  # Bottom-Up IS on the hierarchical part 
  for (i in 1:n_hier) {   # It's important that the rows of H are in the 
                          # correct order (i.e. bottom-up)!
    if (distr == "discrete") {
      w_i = emp_pmf( (B %*% H[i,]), samples_hier[,i] )
      bottom_mask = H[i,] != 0
      B[,bottom_mask] = resample(B[,bottom_mask], w_i)
    }
    if (distr == "continuous") {
      # TODO
    }
  }
  
  if (!is.null(G)) {
    # Plain IS on the additional constraints
    # TODO
  }
  
  A = if (is.null(G)) H else rbind(H, G)
  u = B %*% t(A)
  y = cbind(u, B)
  
  out = list(B_samples=B, U_samples=u, Y_samples=y)
  return(out)
}

system.time({
  recBUISSamples.res = recBUISSamples(
    samples_b = samples_b,
    samples_hier = samples_hier,
    H = H,
    samples_addconstr = NULL,
    G = NULL,
    distr = "discrete"
  )
})


colMeans(recBUIS.res$B_samples)
colMeans(recBUISSamples.res$B_samples)