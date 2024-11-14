# A bivariate pmf is represented as a normalized numeric matrix M: 
# for each (i,j) in {0,...,Mi} x {0,...,Mj}, the probability of (i,j) is the value M[i+1,j+1]

###

# Compute the empirical bivariate pmf from a pair of vectors of samples
bPMF.from_samples = function(v1,v2) {
  .check_discrete_samples(v1)
  .check_discrete_samples(v2)
  if (any(v1<0) | any(v2<0)) {
    stop("Samples must be non-negative")
  }
  
  tab = table(factor(v1, levels = 0:max(v1)),
              factor(v2, levels = 0:max(v2)))
  bpmf = as.matrix(tab)
  rownames(bpmf) = NULL
  colnames(bpmf) = NULL
  
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

