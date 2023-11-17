#' @title Schäfer Strimmer covariance shrinkage 
#'
#' @description
#'
#' Computes the Schäfer Strimmer shrinkage estimator for a covariance matrix
#' from a matrix of samples. 
#'
#' @details
#'
#' This function computes the shrinkage to a diagonal covariance with unequal variances. 
#' Note that here we use the estimators \eqn{S = X X^T/n} and \eqn{T = diag(S)} and we internally
#' use the correlation matrix in place of the covariance to compute the optimal shrinkage factor.
#'
#' @param x matrix of samples with dimensions nxp (n samples, p dimensions).
#'
#' @return A list containing the shrinkage estimator and the optimal lambda. The list has the following named elements:
#'
#' * `shrink_cov`: the shrinked covariance matrix (`p` x `p`);
#' * `lambda_star`: the optimal lambda for the shrinkage;
#'
#' @examples
#'
#' # Generate some multivariate normal samples
#' # Parameters
#' nSamples <- 200
#' pTrue <- 2
#' 
#' # True moments
#' trueSigma <- matrix(c(3,2,2,2), nrow=2)
#' chol_trueSigma <- chol(trueSigma)
#' trueMean <- c(0,0) 
#' 
#' # Generate samples
#' set.seed(42)
#' x <- replicate(nSamples, trueMean) +  t(chol_trueSigma)%*%matrix(rnorm(pTrue*nSamples), 
#'                                                                  nrow=pTrue,ncol=nSamples)
#' x <- t(x) 
#' res_shrinkage <- schaferStrimmer_cov(x)
#' res_shrinkage$lambda_star # should be 0.01287923
#'
#' @references
#' Schäfer, Juliane, and Korbinian Strimmer. (2005). *A Shrinkage Approach to Large-Scale Covariance Matrix Estimation and Implications for Functional Genomics.* Statistical Applications in Genetics and Molecular Biology 4: Article32. \doi{10.2202/1544-6115.1175}.
#' 
#'
#'
#' @export
schaferStrimmer_cov <- function(x){
  n <- nrow(x)
  p <- ncol(x)
  
  # Compute the two estimators
  # Here we used the biased estimator
  sMat <- crossprod(x)/(n)
  
  # Alternative: unbiased estimator
  # sMat <- crossprod(x)/(n-1)
  
  tMat <- diag(diag(sMat))
  
  
  
  # Scale x by the standard deviation
  xscale <- scale(x,center=FALSE,scale=sqrt(diag(sMat)))
  # Compute the correlation equivalents
  rSmat <- stats::cov2cor(sMat)
  rTmat <- stats::cov2cor(tMat)
  
  
  
  # Compute numerator
  # we divide by (1 / (n*(n - 1))) because of biased estimator
  varSij <- (1 / (n*(n - 1))) * (crossprod(xscale ^ 2) - 1 / n * (crossprod(xscale)) ^ 2)
  
  # Alternative for unbiased estimator: divide by (n / ((n - 1)^3))
  # varSij <- (n / ((n - 1)^3)) * (crossprod(xscale ^ 2) - 1 / n * (crossprod(xscale)) ^ 2)
  
  
  diag(varSij) <- 0
  
  # Compute denominator
  sqSij <- (rSmat - rTmat) ^ 2
  
  # Compute Lambda
  lambda_star <- sum(varSij) / sum(sqSij)
  lambda_star <- max(min(lambda_star, 1), 0)
  
  # Compute shrinkage
  shrink_cov <- lambda_star * tMat +  (1 - lambda_star) * sMat
  
  return(list(shrink_cov=shrink_cov, lambda_star=lambda_star))
}
