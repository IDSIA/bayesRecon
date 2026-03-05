# Optimize Degrees of Freedom (nu) via LOO Cross-Validation

Optimize Degrees of Freedom (nu) via LOO Cross-Validation

## Usage

``` r
multi_log_score_optimization(res, prior_mean, trim = 0.1)
```

## Arguments

- res:

  Matrix of residuals (n_obs x n_var).

- prior_mean:

  The prior mean covariance matrix (n_var x n_var).

- trim:

  Fraction of observations to trim (0 to 1). Default is 0.1 (10%).

## Value

A list containing the optimization results:

- `optimal_nu`: The optimal degrees of freedom found.

- `min_neg_log_score`: The minimum negative log score achieved.

- `convergence`: The convergence status of the optimizer.

- `time_elapsed`: The time taken for the optimization.

## Details

TODO: check these details

**Leave-One-Out (LOO) Cross-Validation:** This function estimates the
optimal degrees of freedom \\\nu\\ by maximizing the out-of-sample
predictive performance. This is achieved by computing the log-density of
each held-out observation \\\mathbf{r}\_i\\ given the remaining data
\\\mathbf{R}\_{-i}\\. The total objective function is the sum of these
predictive log-densities: \$\$\mathcal{L}(\nu) = \sum\_{i=1}^T \log
f(\mathbf{r}\_i \| \mathbf{R}\_{-i}, \nu)\$\$

**The Log-Density Function:** For each LOO step, the residuals are
assumed to follow a Multivariate Student-t distribution. The density is
expressed directly as a function of the posterior sum-of-squares matrix
\\\Psi\\, where \\\Psi\\ scales implicitly with \\\nu\\:
\$\$f(\mathbf{r}\_i \| \Psi, \nu) = \frac{\Gamma(\frac{\nu +
T}{2})}{\Gamma(\frac{\nu + T - p}{2}) \pi^{p/2}} \|\Psi\|^{-1/2} \left(
1 + \mathbf{r}\_i^\top \Psi^{-1} \mathbf{r}\_i \right)^{-\frac{\nu +
T}{2}}\$\$ In the code, \\\Psi\\ is constructed as: \$\$\Psi = (\nu -
p - 1)\bar{\Sigma}\_{prior} + \mathbf{R}^\top\mathbf{R}\$\$ By using
this formulation, the standard scaling factors \\1/\nu\\ and
\\\nu^{-p/2}\\ are absorbed into the matrix inverse and determinant,
respectively.

**Efficient Computation via Sherman-Morrison:** Rather than recomputing
\\\Psi\_{-i}\\ and its inverse \\T\\ times, the function uses the
full-sample matrix \\\Psi\\ and adjusts it using the leverage \\h_i =
\mathbf{r}\_i^\top \Psi^{-1} \mathbf{r}\_i\\.

Through the Matrix Determinant Lemma and the Sherman-Morrison formula,
the internal term \\(1 + \mathbf{r}\_i^\top \Psi\_{-i}^{-1}
\mathbf{r}\_i)\\ simplifies to \\(1 - h_i)^{-1}\\. The final log-density
contribution used in the code is: \$\$\log f_i \propto \textrm{const} -
\frac{1}{2}\log\|\Psi\| + \frac{\nu + T - 1}{2} \log(1 - h_i)\$\$
