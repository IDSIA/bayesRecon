# Swiss tourism: monthly Swiss overnight stay data

A dataset of 27 time series with monthly observations of the number of
overnight stays in Switzerland, from Jan 2005 to Jan 2025.

## Usage

``` r
swiss_tourism
```

## Format

A list that contains:

- `ts`: a multivariate time series of class
  [ts](https://rdrr.io/r/stats/ts.html);

- `n_bottom`: integer with the number of bottom time series (26);

- `n_upper`: integer with the number of upper time series (1);

- `agg_mat`: the aggregation matrix for the hierarchy. (1x26)

## Source

Federal Statistical Office of Switzerland (2025). *Hotel sector:
arrivals and overnight stays of open establishments by year, month,
canton and visitors' country of residence*.
<https://www.bfs.admin.ch/asset/en/px-x-1003020000_102>

## Details

This is a hierarchy of time series that contains:

- 26 bottom-level series corresponding to the 26 Swiss cantons

- 1 top-level series corresponding to the Swiss total. Each time series
  has 241 observations.

## References

Carrara, C., Azzimonti, D., Corani, G., Zambon, L. (2025). *Modeling the
uncertainty on the covariance matrix for probabilistic forecast
reconciliation*. arXiv:2506.19554
[doi:10.48550/arXiv.2506.19554](https://doi.org/10.48550/arXiv.2506.19554)
.
