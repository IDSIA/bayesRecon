library(readxl)

dir_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_path)

data        <- read_excel("NB_score_driven_forecasts.xlsx", sheet = "Data")
extr_mkt_events <- ts(data)

usethis::use_data(extr_mkt_events, overwrite = TRUE)
