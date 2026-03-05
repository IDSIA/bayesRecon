# bayesRecon: Probabilistic Reconciliation via Conditioning

Provides methods for probabilistic reconciliation of hierarchical
forecasts of time series. The available methods include analytical
Gaussian reconciliation (Corani et al., 2021)
[doi:10.1007/978-3-030-67664-3_13](https://doi.org/10.1007/978-3-030-67664-3_13)
, MCMC reconciliation of count time series (Corani et al., 2024)
[doi:10.1016/j.ijforecast.2023.04.003](https://doi.org/10.1016/j.ijforecast.2023.04.003)
, Bottom-Up Importance Sampling (Zambon et al., 2024)
[doi:10.1007/s11222-023-10343-y](https://doi.org/10.1007/s11222-023-10343-y)
, methods for the reconciliation of mixed hierarchies (Mix-Cond and
TD-cond) (Zambon et al., 2024)
<https://proceedings.mlr.press/v244/zambon24a.html>, analytical
reconciliation with Bayesian treatment of the covariance matrix (Carrara
et al., 2025) [doi:
10.48550/arXiv.2506.19554](https://doi.org/%2010.48550/arXiv.2506.19554)
.

## Learn more

To learn more about `bayesRecon`, start with the vignettes:
`browseVignettes(package = "bayesRecon")`

## Reconciliation functions

The package implements reconciliation via conditioning for probabilistic
forecasts of hierarchical time series. The main functions are:

- [`reconc_gaussian()`](https://idsia.github.io/bayesRecon/reference/reconc_gaussian.md):
  reconciliation via conditioning, assumes multivariate Gaussian base
  forecasts; this is done analytically;

- [`reconc_BUIS()`](https://idsia.github.io/bayesRecon/reference/reconc_BUIS.md):
  reconciliation via conditioning of any probabilistic forecast via
  importance sampling; this is the recommended option for non-Gaussian
  base forecasts;

- [`reconc_MCMC()`](https://idsia.github.io/bayesRecon/reference/reconc_MCMC.md):
  reconciliation via conditioning of discrete probabilistic forecasts
  via Markov Chain Monte Carlo;

- [`reconc_MixCond()`](https://idsia.github.io/bayesRecon/reference/Mixed_reconciliation.md):
  reconciliation via conditioning of mixed hierarchies, where the upper
  forecasts are multivariate Gaussian and the bottom forecasts are
  discrete distributions;

- [`reconc_TDcond()`](https://idsia.github.io/bayesRecon/reference/Mixed_reconciliation.md):
  reconciliation via top-down conditioning of mixed hierarchies, where
  the upper forecasts are multivariate Gaussian and the bottom forecasts
  are discrete distributions;

- [`reconc_t()`](https://idsia.github.io/bayesRecon/reference/reconc_t.md):
  reconciliation via conditioning with uncertain covariance matrix; the
  reconciled forecasts are multivariate Student-t; this is done
  analytically.

## Utility functions

- [`temporal_aggregation()`](https://idsia.github.io/bayesRecon/reference/temporal_aggregation.md):
  temporal aggregation of a given time series object of class
  [ts](https://rdrr.io/r/stats/ts.html);

- [`get_reconc_matrices()`](https://idsia.github.io/bayesRecon/reference/get_reconc_matrices.md):
  aggregation and summing matrices for a temporal hierarchy of time
  series from user-selected list of aggregation levels;

- [`schaferStrimmer_cov()`](https://idsia.github.io/bayesRecon/reference/schaferStrimmer_cov.md):
  computes the Schäfer-Strimmer shrinkage estimator for the covariance
  matrix;

- [`multi_log_score_optimization()`](https://idsia.github.io/bayesRecon/reference/multi_log_score_optimization.md):
  estimates the optimal degrees of freedom for
  [`reconc_t()`](https://idsia.github.io/bayesRecon/reference/reconc_t.md)
  by maximizing the leave-one-out (LOO) multivariate log density;

- [PMF](https://idsia.github.io/bayesRecon/reference/PMF.md): functions
  for handling PMF objects (sampling, computing statistics like mean,
  variance, quantiles, and summaries).

## References

Corani, G., Azzimonti, D., Augusto, J.P.S.C., Zaffalon, M. (2021).
*Probabilistic Reconciliation of Hierarchical Forecast via Bayes' Rule*.
ECML PKDD 2020. Lecture Notes in Computer Science, vol 12459.
[doi:10.1007/978-3-030-67664-3_13](https://doi.org/10.1007/978-3-030-67664-3_13)
.

Corani, G., Azzimonti, D., Rubattu, N. (2024). *Probabilistic
reconciliation of count time series*. International Journal of
Forecasting 40 (2), 457-469.
[doi:10.1016/j.ijforecast.2023.04.003](https://doi.org/10.1016/j.ijforecast.2023.04.003)
.

Zambon, L., Azzimonti, D. & Corani, G. (2024). *Efficient probabilistic
reconciliation of forecasts for real-valued and count time series*.
Statistics and Computing 34 (1), 21.
[doi:10.1007/s11222-023-10343-y](https://doi.org/10.1007/s11222-023-10343-y)
.

Zambon, L., Agosto, A., Giudici, P., Corani, G. (2024). *Properties of
the reconciled distributions for Gaussian and count forecasts*.
International Journal of Forecasting (in press).
[doi:10.1016/j.ijforecast.2023.12.004](https://doi.org/10.1016/j.ijforecast.2023.12.004)
.

Zambon, L., Azzimonti, D., Rubattu, N., Corani, G. (2024).
*Probabilistic reconciliation of mixed-type hierarchical time series*.
Proceedings of the Fortieth Conference on Uncertainty in Artificial
Intelligence, PMLR 244:4078-4095.
<https://proceedings.mlr.press/v244/zambon24a.html>.

Carrara, C., Corani, G., Azzimonti, D., & Zambon, L. (2025). *Modeling
the uncertainty on the covariance matrix for probabilistic forecast
reconciliation*. arXiv preprint arXiv:2506.19554.
<https://arxiv.org/abs/2506.19554>.

## See also

Useful links:

- <https://github.com/IDSIA/bayesRecon>

- <https://idsia.github.io/bayesRecon/>

- Report bugs at <https://github.com/IDSIA/bayesRecon/issues>

## Author

**Maintainer**: Dario Azzimonti <dario.azzimonti@gmail.com>
([ORCID](https://orcid.org/0000-0001-5080-3061))

Authors:

- Lorenzo Zambon <lorenzo.zambon@idsia.ch>
  ([ORCID](https://orcid.org/0000-0002-8939-993X))

- Stefano Damato <stefano.damato@idsia.ch>
  ([ORCID](https://orcid.org/0009-0001-0225-3736))

- Nicolò Rubattu <nicolo.rubattu@idsia.ch>
  ([ORCID](https://orcid.org/0000-0002-2703-1005))

- Giorgio Corani <giorgio.corani@idsia.ch>
  ([ORCID](https://orcid.org/0000-0002-1541-8384))
