To Do file for package
================
2023-03-30

## MCMC

- line 51 reconc_MH: use function to get distr_b, params_b

- line 31 reconc_MH: use function to check input

- line 94 reconc_MH: add u_samples

- uniform parameters with reconc_IS (create list of distr from a string,
  check rest)

## General

- decide format for samples output (n_samples x n) vs (n x n_samples)

- Create a test for the case where .fix_weights returns a warning.

## Docs

- Finish docs MCMC

- update main package docs to link to MCMC function

- Add basic examples for all main functions

## Pre-CRAN submission

- Add references to all papers in documentation

## Long term to do

- Control on effective sample size to check for weights collapse

- Extension to batch-BUIS?

- check use of parallel computation on local cores for batch IS.
