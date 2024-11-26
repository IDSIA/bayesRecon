# A bivariate pmf is represented as a normalized numeric matrix M: 
# for each (i,j) in {0,...,Mi} x {0,...,Mj}, the probability of (i,j) is the value M[i+1,j+1]

# .KERNELS = c("gauss")
# .KDE_SIGMAS = c("I")

library(ks)

###

# bPMF.get_kde_matrix = function(bpmf, met) {
#   if (!(met %in% .KDE_SIGMAS)) {
#     stop("The KDE covariance is not valid")
#   }
#   
#   if (met == "I") {
#     return(diag(2))  # identity matrix
#   }
# }


bPMF.smoothing = function(bpmf, kde_type, params_kernel = NULL) {
  
  if (kde_type == "gauss") {
    # Create matrix of all points of the grid, 
    # dim: n_points x 2, where n_points = nrow(bpmf)*ncol(bpmf)
    grid = as.matrix(expand.grid(x = 0:(nrow(bpmf)-1), 
                                 y = 0:(ncol(bpmf)-1)))
    
    if (is.null(params_kernel)) {
      H = 1e-5 * diag(2)
    }
    
    # TODO: optimize. compute the cholesky outside the for loop.
    # Is there a way to vectorize everything?
    new_bpmf = rep(0, nrow(grid))
    for (i in 1:nrow(grid)) {
      new_bpmf = new_bpmf + .MVN_density(grid, grid[i,], H)
    }
  }
  
  # # Use ks package
  # if (kde_type == "ks") {
    # grid = expand.grid(x = 0:(nrow(bpmf)-1),
    #                    y = 0:(ncol(bpmf)-1))
  #   samples = bPMF.sample(bpmf, 1e5)
    # kde_ = ks::kde(x = samples,
    #                eval.points = grid,
    #                binned=FALSE)
  #   new_bpmf = matrix(kde_$estimate, nrow = nrow(bpmf))
  # }
  
  return(new_bpmf)
  
}


# Funzione di prova, usa pacchetto ks
bPMF.from_samples_kde = function(v1, v2) {
  
  # Check that samples are discrete and non-negative
  .check_discrete_samples(v1)
  .check_discrete_samples(v2)
  if (any(v1[!is.na(v1)]<0) | any(v2[!is.na(v2)]<0)) {
    stop("Samples must be non-negative")
  }
  
  # reduced vectors: I only keep indices for which both vectors have values 
  v1_red = v1[!is.na(v1) & !is.na(v2)]
  v2_red = v2[!is.na(v1) & !is.na(v2)]
  
  # KDE
  grid = as.matrix(expand.grid(x = 0:max(v1_red), 
                               y = 0:max(v2_red)))
  kde_ = ks::kde(x = cbind(v1_red, v2_red),
                 eval.points = grid,
                 binned=FALSE)
  bpmf_kde = matrix(kde_$estimate, nrow = max(v1_red)+1)
  bpmf_kde = bpmf_kde / sum(bpmf_kde)

  
  # TODO: where do i do KDE?
  # If here: avoid holes in the bpmf that cannot be filled with the marginals correction
  # however then the marginal are re-scaled by something that has not been smoothened...
  
  # Correct the marginals using values that have been ignored before
  # If there are values in v1 for which there is NA on v2:
  if (sum(!is.na(v1)) > length(v1_red)) {  
    new_marginal1 = ks::kde(v1[!is.na(v1)], eval.points = 0:max(v1, na.rm=T))
    new_marginal1 = new_marginal1 / sum(new_marginal1)
    bpmf = bpmf / rowSums(bpmf) * new_marginal1  # divide and multiply each row
    bpmf[is.na(bpmf)] = 0  # because NA appear where rowsum is zero 
  }
  # Same for v2
  if (sum(!is.na(v2)) > length(v2_red)) {  
    new_marginal2 = ks::kde(v2[!is.na(v2)], eval.points = 0:max(v2, na.rm=T))
    new_marginal2 = new_marginal2 / sum(new_marginal2)
    bpmf =  bpmf %*% diag(new_marginal2/colSums(bpmf))  # divide and multiply each column
    bpmf[is.na(bpmf)] = 0  # because NA appear where colsum is zero 
  }

  return(bpmf)
}

# Compute the empirical bivariate pmf from a pair of vectors of samples
# v1 and v2 have to be of the same length. 
# They can, however, have NA, corresponding to values that are not observed
# First, compute the bivariate empirical distribution using only pairs where
# both values are not NA
# Then, correct each marginal by looking at all the available data for that vector
# Returns a bpmf (where v1 is on the rows, v2 on the columns)
bPMF.from_samples = function(v1, v2, smoothing = TRUE, kde_type = "gauss") {
  
  # Check that samples are discrete and non-negative
  .check_discrete_samples(v1)
  .check_discrete_samples(v2)
  if (any(v1[!is.na(v1)]<0) | any(v2[!is.na(v2)]<0)) {
    stop("Samples must be non-negative")
  }
  
  # reduced vectors: I only keep indices for which both vectors have values 
  v1_red = v1[!is.na(v1) & !is.na(v2)]
  v2_red = v2[!is.na(v1) & !is.na(v2)]
  
  tab = table(factor(v1_red, levels = 0:max(v1[!is.na(v1)])),
              factor(v2_red, levels = 0:max(v2[!is.na(v2)])))
  bpmf = as.matrix(tab)
  bpmf = bpmf / sum(bpmf)
  rownames(bpmf) = NULL
  colnames(bpmf) = NULL
  
  # TODO: where do i do KDE?
  # If here: avoid holes in the bpmf that cannot be filled with the marginals correction
  # however then the marginal are re-scaled by something that has not been smoothened...
  
  # Correct the marginals using values that have been ignored before
  # If there are values in v1 for which there is NA on v2:
  if (sum(!is.na(v1)) > length(v1_red)) {  
    new_marginal1 = PMF.from_samples(v1[!is.na(v1)])
    bpmf = bpmf / rowSums(bpmf) * new_marginal1  # divide and multiply each row
    bpmf[is.na(bpmf)] = 0  # because NA appear where rowsum is zero 
  }
  # Same for v2
  if (sum(!is.na(v2)) > length(v2_red)) {  
    new_marginal2 = PMF.from_samples(v2[!is.na(v2)])
    bpmf =  bpmf %*% diag(new_marginal2/colSums(bpmf))  # divide and multiply each column
    bpmf[is.na(bpmf)] = 0  # because NA appear where colsum is zero 
  }
  
  # Kernel density estimation
  if (smoothing) {
    bpmf = bPMF.smoothing(bpmf, kde_type)
    bpmf = bpmf / sum(bpmf)
  }
  
  return(bpmf)
}


#' @title Sample from the bivariate distribution given as a bPMF object. TODO!!
#'
#' @description
#' 
#' Samples (with replacement) from the bivariate probability distribution specified by `bpmf`.
#' 
#' @param bpmf the bPMF object.
#' @param N_samples number of samples. 
#' 
#' @return Matrix N_samples x 2 of bivariate samples drawn from the distribution 
#' specified by `bpmf`.  
#' 
#' @export
bPMF.sample = function(bpmf, N_samples) {
  Mi = nrow(bpmf)
  Mj = ncol(bpmf)
  s = sample(1:(Mi * Mj), 
             prob = as.vector(bpmf), 
             replace = TRUE, size = N_samples)
  si = (s-1) %% Mi
  sj = (s-1) %/% Mi
  S = matrix(c(si,sj), ncol=2)
  return(S)
}

