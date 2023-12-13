#' Extreme market events dataset
#'
#' Count time series of extreme market events in five economic sectors.
#' The data refer to the trading days between 2004/12/31 and 2018/12/19 (3508 trading days in total).
#' 
#' The counts are computed by considering 29 companies included in the Euro Stoxx 
#' 50 index and observing if the value of the CDS spread on a given day exceeds 
#' the 90-th percentile of its distribution in the last trading year.
#' The companies are divided in the following  sectors: Financial (FIN), Information 
#' and Communication Technology (ICT), Manufacturing (MFG), Energy (ENG), and Trade (TRD). 
#' 
#' There are 6 time series: 
#' - 5 bottom time series, corresponding to the daily counts for each sector
#' - 1 upper time series, which is the sum of all the bottom (ALL)
#'
#' @format
#' A multivariate time series of class \link[stats]{ts}.
#' 
#' @source 
#' Zambon, L., Agosto, A., Giudici, P., Corani, G. (2023). *Properties of the reconciled distributions for Gaussian and count forecasts*. \doi{10.48550/arXiv.2303.15135}.
#' 
#' @references 
#' Zambon, L., Agosto, A., Giudici, P., Corani, G. (2023). *Properties of the reconciled distributions for Gaussian and count forecasts*. \doi{10.48550/arXiv.2303.15135}.
#'
#' Agosto, A. (2022). *Multivariate Score-Driven Models for Count Time Series to Assess Financial Contagion*. \doi{10.2139/ssrn.4119895}
"extr_mkt_events"
