library(readxl)

dir_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_path)

predictions <- read_excel("NB_score_driven_forecasts.xlsx", sheet = "Predictions")
alpha       <- read_excel("NB_score_driven_forecasts.xlsx", sheet = "alpha", n_max = 1)

extr_mkt_events_basefc <- list(mu   = data.frame(predictions),
                               size = data.frame(1 / alpha))

usethis::use_data(extr_mkt_events_basefc, overwrite = TRUE)
