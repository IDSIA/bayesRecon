# bayesRecon 0.3.2

* Fixed issues in tests that would run into compatibility issues with the new version of waldo.

* Added new tests for some limit usage of mixCond, TDcond and BUIS. 

* Added GitHub URL and BugReports URL to DESCRIPTION. 

# bayesRecon 0.3.1

* IMPORTANT CHANGE IN THE API OF THE `reconc_*` functions: they now require the aggregation matrix A and not the summing matrix S.

* The examples section of the `reconc_TDcond` now contains an example showing how to handle the case of an unbalanced hierarchy. 

# bayesRecon 0.3.0

* Added `reconc_MixCond`, the implementation of Mixed conditioning for the reconciliation of mixed-type hierarchical forecasts.

* Added `reconc_TDcond`, the implementation of Top-down conditioning for the reconciliation of mixed-type hierarchical forecasts.

* Vignette "Reconciliation of M5 hierarchy with mixed-type forecasts".

* Added functions `PMF.get_mean`, `PMF.get_quantile`, `PMF.get_var`, `PMF.sample`, `PMF.summary` that return point estimates and samples from a probability mass function object.

* Added the dataset `M5_CA1_basefc` which contains the base forecasts for the store "CA1" of the hierarchical time series in the M5 competition. 

# bayesRecon 0.2.0

* Vignette "Properties of the reconciled distribution via conditioning".

* Added option in `reconc_BUIS` to pass some base forecast as parameters and some as samples.

* Added option in `reconc_BUIS` to input a list of distributions instead of a string.

* Fixed bugs and closed github issues.

# bayesRecon 0.1.2

* Vignette "Probabilistic Reconciliation via Conditioning with bayesRecon".

* Added the `schaferStrimmer_cov` function.

* Fixed bugs.


# bayesRecon 0.1.1

* Added a `NEWS.md` file to track changes to the package.

* Fixed issue that could cause test fails on some machines.

# bayesRecon 0.1.0

* Initial release of the package.
