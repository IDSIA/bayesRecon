source("reconc.R")

S = t(matrix(data=c(
  1,0,0,0,0,0,0,0,
  0,1,0,0,0,0,0,0,
  0,0,1,0,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,0,1,0,0,0,
  0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,1,
  1,1,0,0,0,0,0,0,
  0,0,1,1,0,0,0,0,
  0,0,0,0,1,1,0,0,
  0,0,0,0,0,0,1,1,
  1,1,1,1,0,0,0,0,
  0,0,0,0,1,1,1,1,
  1,1,1,1,1,1,1,1), nrow=8, ncol=15))

base_forecasts = list(
  list(5.634849165190505),
  list(9.833589192410015),
  list(6.3023800293289165),
  list(9.486182621822866),
  list(6.883748580948357),
  list(6.681108721672265),
  list(7.256882352376998),
  list(9.201275416306906),
  list(20.108969864880674),
  list(20.525131446497316),
  list(17.63431449340681),
  list(21.395605099289075),
  list(40.63410131137799),
  list(39.02991959269589),
  list(79.66402090407388))

### Reconciliation using BUIS with parameters
system.time({
  recBUIS.res = reconc(
    S=S,
    base_forecasts = base_forecasts,
    in_type = "params",
    distr = "poisson"
  )
})

### Reconciliation using BUIS with samples
path = "../../LorenzoRepoBUIS/"
samples_b = reticulate::py_load_object(filename = paste0(path,"samples_b.pkl"))
samples_hier = reticulate::py_load_object(filename = paste0(path,"samples_hier.pkl"))
base_forecasts_samples = apply(cbind(samples_b, samples_hier), 2, as.list)

system.time({
  recBUISSamples.res = reconc(
    S=S,
    base_forecasts = base_forecasts_samples,
    in_type = "samples",
    distr = "poisson"
  )
})

colMeans(recBUIS.res$bottom_reconciled_samples)
colMeans(recBUISSamples.res$bottom_reconciled_samples)
