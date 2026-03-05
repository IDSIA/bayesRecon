# Schäfer Strimmer covariance shrinkage

Computes the Schäfer Strimmer shrinkage estimator for a covariance
matrix from a matrix of samples.

## Usage

``` r
schaferStrimmer_cov(x)
```

## Arguments

- x:

  matrix of samples with dimensions nxp (n samples, p dimensions).

## Value

A list containing the shrinkage estimator and the optimal lambda. The
list has the following named elements:

- `shrink_cov`: the shrinked covariance matrix (`p` x `p`);

- `lambda_star`: the optimal lambda for the shrinkage;

## Details

This function computes the shrinkage to a diagonal covariance with
unequal variances. Note that here we use the estimators \\S = X X^T/n\\
and \\T = diag(S)\\ and we internally use the correlation matrix in
place of the covariance to compute the optimal shrinkage factor.

## References

Schäfer, Juliane, and Korbinian Strimmer. (2005). *A Shrinkage Approach
to Large-Scale Covariance Matrix Estimation and Implications for
Functional Genomics.* Statistical Applications in Genetics and Molecular
Biology 4: Article32.
[doi:10.2202/1544-6115.1175](https://doi.org/10.2202/1544-6115.1175) .

## Examples

``` r
# Generate some multivariate normal samples
# Parameters
nSamples <- 200
pTrue <- 2

# True moments
true_cov <- matrix(c(3, 2, 2, 2), nrow = 2)
chol_true_cov <- chol(true_cov)
true_mean <- c(0, 0)

# Generate samples
set.seed(42)
x <- replicate(nSamples, true_mean) +
  t(chol_true_cov) %*% matrix(stats::rnorm(pTrue * nSamples),
    nrow = pTrue, ncol = nSamples
  )
x <- t(x)
res_shrinkage <- schaferStrimmer_cov(x)
res_shrinkage$lambda_star # should be 0.01287923
#> [1] 0.01287923
```
