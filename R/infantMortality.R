#' Infant Mortality grouped time series dataset
#'
#' A yearly grouped time series dataset, from 1901 to 2003, of infant mortality counts (deaths) in Australia; 
#' disaggregated by state (see below), and sex (male and female).
#' 
#' States: New South Wales (NSW), Victoria (VIC), Queensland (QLD), South Australia (SA), Western Australia 
#' (WA), Northern Territory (NT), Australian Capital Territory (ACT), and Tasmania (TAS).
#'
#' @format
#' List of time series of class \link[stats]{ts}.
#' 
#' @source 
#' hts package [CRAN](https://cran.r-project.org/package=hts)
#' 
#' @references 
#' R. J. Hyndman, R. A. Ahmed, G. Athanasopoulos and H.L. Shang (2011) Optimal combination forecasts for hierarchical time series. Computational Statistics and Data Analysis, 55(9), 2579-2589.
"infantMortality"