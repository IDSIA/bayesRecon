# Core Reconciliation via Multivariate Student-t Distribution.

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
  return_parameters = FALSE,
  suppress_warnings = FALSE
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
  posterior Student-t distribution.

- nu_post:

  Degrees of freedom of the posterior Student-t distribution.

- return_upper:

  Logical, whether to return the reconciled parameters for the upper
  variables (default is FALSE).

- return_parameters:

  Logical. If TRUE, returns internal parameters (posterior nu, posterior
  Psi, C) for debugging or advanced use. Default is FALSE.

- suppress_warnings:

  Logical. If TRUE, suppresses warnings about numerical issues. Default
  is FALSE.

## Value

A list containing:

- `bottom_rec_mean`: Reconciled mean vector for bottom level.

- `bottom_rec_scale_matrix`: Reconciled scale matrix for bottom level.

- `bottom_rec_df`: Reconciled degrees of freedom for bottom level.

- `upper_rec_mean`: (only if `return_upper=TRUE`) Reconciled mean vector
  for upper level.

- `upper_rec_scale_matrix`: (only if `return_upper=TRUE`) Reconciled
  scale matrix for upper level.

- `upper_rec_df`: (only if `return_upper=TRUE`) Reconciled degrees of
  freedom for upper level.

- `posterior_nu`: (only if `return_parameters=TRUE`) Posterior degrees
  of freedom.

- `posterior_Psi`: (only if `return_parameters=TRUE`) Posterior scale
  matrix.

- `C`: (only if `return_parameters=TRUE`) Scaling factor used in
  reconciliation.
