# Core Reconciliation via Bayesian Universality Information Sharing

Internal function that performs the core reconciliation logic for the
BUIS method, which reconciles forecasts using importance sampling based
on both hierarchical and equality constraints through a Bayesian
framework.

## Usage

``` r
.core_reconc_BUIS(
  A,
  H,
  G,
  B,
  upper_base_fc_H,
  in_typeH,
  distr_H,
  upper_base_fc_G,
  in_typeG,
  distr_G,
  .comp_w = .compute_weights,
  suppress_warnings = FALSE,
  return_upper = TRUE
)
```

## Arguments

- A:

  Matrix defining the overall hierarchy.

- H:

  Matrix defining hierarchical constraints.

- G:

  Matrix defining general linear constraints.

- B:

  Matrix of bottom level base forecast samples.

- upper_base_fc_H:

  List of upper base forecasts for hierarchical constraints.

- in_typeH:

  Character string specifying input type for H forecasts ('pmf',
  'samples', or 'params').

- distr_H:

  Character string specifying distribution type for H forecasts
  ('poisson' or 'nbinom').

- upper_base_fc_G:

  List of upper base forecasts for general constraints.

- in_typeG:

  Character string specifying input type for G forecasts ('pmf',
  'samples', or 'params').

- distr_G:

  Character string specifying distribution type for G forecasts
  ('poisson' or 'nbinom').

- .comp_w:

  Function to compute weights for importance sampling. Default is
  `.compute_weights`.

- suppress_warnings:

  Logical. If TRUE, suppresses warnings about sample quality. Default is
  FALSE.

## Value

A list containing:

- `bottom_rec`: List with reconciled bottom forecasts (pmf and/or
  samples).

- `upper_rec_H`: List with reconciled upper forecasts for H constraints.

- `upper_rec_G`: List with reconciled upper forecasts for G constraints.
