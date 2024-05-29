get_hier_M5 = function(save_ = FALSE, save_path = "../../../results/") {
  
  n_b = 3049
  n_u = 11
  n = n_u + n_b
  
  A = matrix(0, nrow = n_u, ncol = n_b)
  A[1,] = 1           # store
  A[2,1:565] = 1      # categories
  A[3,566:1612] = 1
  A[4,1613:3049] = 1
  A[5,1:416] = 1      # departments
  A[6,417:565] = 1
  A[7,566:1097] = 1
  A[8,1098:1612] = 1
  A[9,1613:1828] = 1
  A[10,1829:2226] = 1
  A[11,2227:3049] = 1
  
  S = rbind(A, diag(rep(1, n_b)))
  
  if (save_) {
    saveRDS(S, paste0(save_path, "S.rds"))
    saveRDS(A, paste0(save_path, "A.rds"))
  }
  
  return(list(A = A, S = S))
}

get_upp_ts = function(dset.store.train, dset.store.test, lev, 
                      name, h, len = 1941) {
  
  if (!(lev %in% c("store_id", "cat_id", "dept_id"))) stop("Wrong lev name")
  
  df1 = dset.store.train[dset.store.train[[lev]] == name]
  df2 = dset.store.test[dset.store.test[[lev]] == name]
  st = which(colnames(df1)=="d_1")
  
  serie = cbind(
    df1[,st:(st+len-1)],
    df2[,st:(st+h-1)] )
  data.train = serie[,1:(len+h-1)]
  data.test = serie[,(len+h):(len+h)]
  train.agg = colSums(data.train)
  test.agg = colSums(data.test)
  return(list(train = train.agg, test = test.agg))
}

get_bott_ts = function(dset.store.train, dset.store.test, 
                       item_id, h, len = 1941) {
  
  df1 = dset.store.train[dset.store.train$item_id == item.id]
  df2 = dset.store.test[dset.store.test$item_id == item.id]
  st = which(colnames(df1)=="d_1")
  
  serie = c(df1[,st:(st+len-1)], df2[,st:(st+h-1)] )
  train = as.numeric(serie)[1:(len+h-1)]
  test = as.numeric(serie)[[len+h]]
  return(list(train = train, test = test))
}

model_upper = function(data, mod = "AXX") {
  model = auto.adam(data, mod, lags = c(7), distribution = "dnorm",
                    occurrence = "none")
  # uses: orders = list(ar = c(3, 3), i = c(2, 1), ma = c(3, 3), select = TRUE)
  return(model)
}

model_bottom = function(data, model_str = "MNN", occ_str = "auto", distr = c("dgamma")) {
  # distr = c("dnorm","dlaplace","ds","dgnorm", "dlnorm", "dinvgauss", "dgamma")
  model = adam(data, model_str, lags = c(7), 
               occurrence = occ_str, distribution = distr)
  return(model)
}