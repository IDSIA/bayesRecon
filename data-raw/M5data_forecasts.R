rm(list = ls())
library(m5)
library(smooth)


source("./data-raw/M5_utils.R")

##################################

seed=42

set.seed(seed)
STORE="CA_1"
n_samples_b = 2e4
round_up = FALSE

path_m5_data <- "./data-raw/M5_data"
m5_forecast_path <- "./data-raw/M5_for_data/"

m5::m5_download(path_m5_data)

dset = m5::m5_get_raw_evaluation(path = path_m5_data)
names(dset) = c("sales_train_evaluation", "sales_test_evaluation", 'sell_prices',
                "calendar","weights_evaluation")
dset.store.train = dset$sales_train_evaluation
dset.store.train = dset.store.train[dset.store.train$store_id == STORE]
dset.store.test = dset$sales_test_evaluation
dset.store.test = dset.store.test[dset.store.test$store_id == STORE]

CAT = c("HOBBIES", "HOUSEHOLD", "FOODS")
DEPT = c("HOBBIES_1", "HOBBIES_2", "HOUSEHOLD_1", "HOUSEHOLD_2", 
         "FOODS_1", "FOODS_2", "FOODS_3")

store_path = paste0(m5_forecast_path, STORE)
len = 1941

###############
### Base forecasts upper time series
print("Computing base forecasts of upper time series...")

h=1

if (!dir.exists(store_path)) dir.create(store_path,recursive = TRUE)

str_base_fc = paste0(store_path,"/CA_1_h_",h,"_base_fc.rds")


M5_CA1_basefc = list()

train_u <- list()

# Store
ts.agg = get_upp_ts(dset.store.train, dset.store.test, "store_id", 
                    STORE, h = h, len = len)
train.agg = ts.agg$train
test.agg = ts.agg$test

model = model_upper(train.agg)       
fc.model = forecast(model, h = 1) 
M5_CA1_basefc$upper[[STORE]] = list(
  mu = fc.model$mean,
  sigma = model$scale,
  actual = test.agg,
  residuals = model$residuals
)
train_u[[STORE]] <- train.agg

# Category
for (cat.id in CAT) {
  ts.agg = get_upp_ts(dset.store.train, dset.store.test, "cat_id", 
                      cat.id, h = h, len = len)
  train.agg = ts.agg$train
  test.agg = ts.agg$test
  
  model = model_upper(train.agg)       
  fc.model = forecast(model, h = 1) 
  M5_CA1_basefc$upper[[cat.id]] = list(
    mu = fc.model$mean,
    sigma = model$scale,
    actual = test.agg,
    residuals = model$residuals
  )
  train_u[[cat.id]] <- train.agg
}

# Department
for (dept.id in DEPT) {
  ts.agg = get_upp_ts(dset.store.train, dset.store.test, "dept_id", 
                      dept.id, h = h, len = len)
  train.agg = ts.agg$train
  test.agg = ts.agg$test
  
  model = model_upper(train.agg)       
  fc.model = forecast(model, h = 1) 
  M5_CA1_basefc$upper[[dept.id]] = list(
    mu = fc.model$mean,
    sigma = model$scale,
    actual = test.agg,
    residuals = model$residuals
  )
  train_u[[dept.id]] <- train.agg
}


saveRDS(M5_CA1_basefc, str_base_fc)


train_b <- list()

for (item.id in unique(dset.store.train$item_id)) {
  
  bts = get_bott_ts(dset.store.train, dset.store.test, 
                    item_id, h, len = 1941)
  train = bts$train
  test = bts$test
  
  model = model_bottom(train, model_str = "MNN",
                       occ_str = "auto", #"odds-ratio", 
                       distr = "dgamma")
  fc.model = forecast(model, h = 1, interval="simulated",
                      scenarios=TRUE, nsim = n_samples_b)
  
  # round to integer (up or to the closest integer, depending on round_up)
  samples = if (round_up) ceiling(fc.model$scenarios[1,]) else round(fc.model$scenarios[1,])
  samples[samples<0] = 0    # set negative to zero
  samples <- as.integer(samples)
  pmf = PMF.from_samples(samples)  # empirical pmf
  
  M5_CA1_basefc$bottom[[item.id]] = list(
    pmf = pmf,
    actual = test
    #residuals = model$residuals
  )
  train_b[[item.id]] <- train 
}


hier = get_hier_M5(save_ = FALSE)

M5_CA1_basefc$A <- hier$A
M5_CA1_basefc$S <- hier$S


Q_u = unlist(lapply(train_u, function(x) mean(abs(x[-1] - x[-length(x)]))))
Q_b = unlist(lapply(train_b, function(x) mean(abs(x[-1] - x[-length(x)]))))

M5_CA1_basefc$Q_u <- Q_u
M5_CA1_basefc$Q_b <- Q_b


usethis::use_data(M5_CA1_basefc, overwrite = TRUE)

unlink("./data-raw/M5_for_data/",recursive=TRUE)
unlink("./data-raw/M5_data/",recursive=TRUE)
