To Do file for package
================
2023-03-30

## General

- Create a test for the case where .fix_weights returns a warning.

- Decide on notation upper vs aggregate

- Add start in temporal aggregation

- Remove bullet points from all argument descriptions

- Vignettes: comparison between MCMC and BUIS on different hierarchies

- tscounts

## Docs

- Docs MCMC: add link/reference for MCMC convergence checks

- update README with badges, installation

- Once we make github public: add getting help and contribute sections
  to README

## Long term to do

- add `usethis::use_github_action()`

- Change documentation of reconc_buis and reconc_mcmc: distr can be a
  list of strings (also in_type?)

- Check on list of distributions (if a bottom is continuous, then all
  the upper must be continuous)

- Control on effective sample size to check for weights collapse

- Extension to batch-BUIS?

- check use of parallel computation on local cores for batch BUIS.

- Add zero-inflation

- Add “fixed forecast” (i.e. Dirac delta)
