# Core reconciliation via multivariate t-distribution.

Internal function that performs the core reconciliation logic for the
t-distribution based reconciliation method. This function assumes an
uncertain covariance matrix with an Inverse-Wishart prior.

## Usage

``` r
.core_reconc_t(
  A,
  base_fc_mean,
  Psi_post,
  nu_post,
  return_upper = FALSE,
  return_parameters = FALSE
)
```

## Arguments

- A:

  Matrix (n_upper x n_bottom) defining the hierarchy where upper = A
  %\*% bottom.

- base_fc_mean:

  Vector of length (n_upper + n_bottom) containing the base forecast
  means for both upper and bottom levels (upper first, then bottom).

- Psi_post:

  Scale matrix (n_upper + n_bottom x n_upper + n_bottom) of the
  posterior multivariate t-distribution.

- nu_post:

  Degrees of freedom of the posterior multivariate t-distribution.

- return_upper:

  Logical, whether to return the reconciled parameters for the upper
  variables (default is FALSE).

## Value

A list containing:

- `bottom_rec_mean`: reconciled bottom-level mean.

- `bottom_rec_scale_matrix`: reconciled bottom-level scale matrix.

- `bottom_rec_df`: reconciled degrees of freedom.

If `return_upper` is TRUE, also returns:

- `upper_rec_mean`: reconciled upper-level mean.

- `upper_rec_scale_matrix`: reconciled upper-level scale matrix.

- `upper_rec_df`: reconciled upper-level degrees of freedom.
