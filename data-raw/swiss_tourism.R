# Swiss Tourism data generator
rm(list=ls())

dir_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_path)


#Function for reading and formatting the Swiss Tourism dataset

#INPUT
#file_name: the "SwissTourism.csv" file

#OUTPUT
# the list to be saved as data in the package:
# ts: the data as time series object
# n_bottom: integer with # of bottoms
# n_upper: integer with # of uppers 
# agg_mat: aggregation matrix for the hierarchy
format_SwissDataset = function(data_, end_month=1, end_year=2025){
  current_year = format(Sys.Date(), "%Y")
  current_month = format(Sys.Date(), "%m")
  
  years = unique(data_$Anno)
  month_names = c("Gennaio",   "Febbraio",  "Marzo",     "Aprile",    "Maggio",
                  "Giugno",    "Luglio",    "Agosto",   "Settembre", "Ottobre",
                  "Novembre",  "Dicembre" )
  months = match(data_$Mese, month_names)
  start_year = as.numeric(years[1])
  start_month = as.numeric(months[1])
  
  bottom = unique(data_$Cantone)[2:length(unique(data_$Cantone))]
  upper = unique(data_$Cantone)[1]
  
  mdata = matrix(t(as.numeric(data_$Pernottamenti[!is.na(as.numeric(data_$Pernottamenti))])), 
                 ncol = 27, byrow = T)
  colnames(mdata) = c(upper, bottom)
  
  dfData = tibble::as_tibble(mdata)
  
  data_ts = ts(dfData,start=c(as.numeric(start_year), as.numeric(start_month)),
               frequency=12)
  
  if(!is.null(end_year) || !is.null(end_month)){
    data_ts <- window(data_ts, end=c(as.numeric(end_year), as.numeric(end_month)))
  }
  
  n_bottom = 26
  n_upper = 1
  ts_length = nrow(data_ts)
  freq = frequency(data_ts)

    
  A = matrix(1, nrow = n_upper, ncol = n_bottom, dimnames = list(upper,
                                                                 bottom))
  
  out = list(ts = data_ts, n_bottom = n_bottom, n_upper = n_upper, agg_mat = A)
  
  
  return(out)
}




library(pxweb)
library(tidyverse)

# The working directory is data-raw (see setwd() above), so the file is read
# from and written to that folder.
file_name <- "SwissTourism.csv"

# Set to TRUE to download a fresh copy of the data from the BFS API.
# NOTE: the BFS revises the published figures, so a new download does NOT
# reproduce the dataset currently shipped in the package: it adds the months
# published since the last download and updates the provisional values.
download_data <- FALSE

if (download_data) {
  # URL of the BFS API
  url <- "https://www.pxweb.bfs.admin.ch/api/v1/it/px-x-1003020000_102/px-x-1003020000_102.px"

  # Building the query
  query <- list(
    "Jahr" = c("*"),
    "Monat" = c("*"),
    "Kanton" = c("*"),
    "Herkunftsland" = c("00"),       # Tutte le provenienze
    "Indikator" = c("2")            # Pernottamenti
  )

  # Download the data
  data <- pxweb_get(url, query = query, verbose = TRUE)

  # Geneate the dataframe
  df <- as.data.frame(data)

  # make the tibble
  tb_swiss <- df |> as_tibble()

  # filter the tibble, remove useless columns and rename the column with the data
  tb_swiss <- tb_swiss |> filter(Mese != "Totale dell'anno") |>
    select(-c("Indicatore", "Paese di provenienza")) |>
    rename("Pernottamenti" = `Settore alberghiero: arrivi e pernottamenti degli stabilimenti aperti`)

  # Write the result
  write_csv(tb_swiss, file_name)
}

data_ = read_csv(file_name)
swiss_tourism <- format_SwissDataset(data_)

usethis::use_data(swiss_tourism, overwrite = TRUE)
