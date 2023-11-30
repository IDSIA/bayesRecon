library(readxl)

dir_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_path)

data        <- read_excel("NB_score_driven_forecasts.xlsx", sheet = "Data")
predictions <- read_excel("NB_score_driven_forecasts.xlsx", sheet = "Predictions")
alpha       <- read_excel("NB_score_driven_forecasts.xlsx", sheet = "alpha", n_max = 1)

base_fc <- list(
  mu = predictions,
  size = 1 / alpha
)

extreme_market_events <- list(actuals = data,
                              base_fc = base_fc)

usethis::use_data(extreme_market_events, overwrite = TRUE)
