# This script runs an example for the 'bayesReco' package.
# A monthly count time series is aggregated and forecasted at different time
# scales. Reconciliation is performed given discrete samples.
# The script uses the following external packages:
# - [X] thief: to build temporal hierarchy
# - tscount: to model count data with Poisson
# - scoringRules: to compute CRPS
# - ggplot2: to plot

library(ggplot2)
library(ggfortify)
library(tscount)
library(bayesReco)
library(scoringRules)

set.seed(42)

### Load data
# Time series from McCleary and Hay (1980). The data are the number of purses
# snatched in the Hyde Park neighborhood of Chicago from Jan 1969 to Nov 1974
HydePark = data.frame(
  Snatched.Purses = c(10,15,10,10,12,10,7,7,10,14,8,17,14,18,3,9,11,10,6,12,14,
                      10,25,29,33,33,12,19,16,19,19,12,34,15,36,29,26,21,17,19,
                      13,20,24,12,6,14,6,12,9,11,17,12,8,14,14,12,5,8,10,3,16,
                      8,8,7,12,6,10,8,10,5,7))
HydePark.ts = ts(data = HydePark, frequency = 12, start = c(1969, 1))
ggplot2::autoplot(HydePark.ts)

### Here we use a simple train/test set
train = window(HydePark.ts, end = c(1973, 12))
test = window(HydePark.ts, start = c(1974))

### Build hierarchy
aggf = c(1,2,3,4,6,12)
train.aggregate = bayesReco::temporal_aggregation(train, aggf)
levels <- attributes(train.aggregate)$names

### Hierarchical forecasting
# n.b. here we use 'tscount' package
fc.samples = list()
fc.count = 1
for (l in seq_along(train.aggregate)) {
  f.level <- frequency(train.aggregate[[l]])
  print(paste("Forecasting at ", levels[l], "..", sep = ""))
  model <- tscount::tsglm(
    ts = train.aggregate[[l]],
    model = list(past_obs = 1, past_mean = NULL),
    distr = "poisson",
    link = "log"
  )
  H <- if (f.level == 1) 2 else f.level
  preds_obj <- predict(model,
                       n.ahead = H,
                       level = 0.9,
                       B = 20000)
  preds_samples <- preds_obj$samples
  for (h in seq(1:f.level)) {
    fc.samples[[fc.count]] <- preds_samples[h,]
    fc.count = fc.count + 1
  }
}

### Reconciliation
tmp = bayesReco::reconc_matrices(aggf, bottom.f = frequency(train), bottom.H = 12)
S = tmp$S
A = tmp$A
reconc.res = bayesReco::reconc_IS(S, base_forecasts = fc.samples, in_type = "samples", distr = "discrete")

### Evaluate bottom time series
# - MAPE for point forecasts
# - CRPS for distributions
ape.fc = list()
ape.reconc = list()
crps.fc = list()
crps.reconc = list()

y.hat = list()
y.hat.lo90 = list()
y.hat.hi90 = list()
y.reconc = list()
y.reconc.lo90 = list()
y.reconc.hi90 = list()

for (h in 1:length(test)) {
  y.hat_ = median(fc.samples[[nrow(A)+h]])
  y.reconc_ = median(reconc.res$bottom_reconciled_samples[, h])
  # Compute Absolute Percentage Errors (APE)
  ape.fc[[h]] = abs(test[h] - y.hat_) / abs(test[h])
  ape.reconc[[h]] = abs(test[h] - y.reconc_) / abs(test[h])
  # Compute Continuous Ranked Probability Score (CRPS)
  crps.fc[[h]] = crps_sample(y = test[h], dat = fc.samples[[nrow(A)+h]])
  crps.reconc[[h]] = crps_sample(y = test[h], dat = reconc.res$bottom_reconciled_samples[, h])

  y.hat[[h]] = y.hat_
  y.hat.lo90[[h]] = quantile(fc.samples[[nrow(A)+h]], (1 - 0.9) / 2)[[1]]
  y.hat.hi90[[h]] = quantile(fc.samples[[nrow(A)+h]], 1 - (1 - 0.9) / 2)[[1]]
  y.reconc[[h]] = y.reconc_
  y.reconc.lo90[[h]] = quantile(reconc.res$bottom_reconciled_samples[, h], (1 -
                                                                              0.9) / 2)[[1]]
  y.reconc.hi90[[h]] = quantile(reconc.res$bottom_reconciled_samples[, h], 1 -
                                  (1 - 0.9) / 2)[[1]]
}

mape.fc = mean(unlist(ape.fc))
mape.reconc = mean(unlist(ape.reconc))
crps.fc = mean(unlist(crps.fc))
crps.reconc = mean(unlist(crps.reconc))
metrics = data.frame(row.names = c("MAPE", "CRPS"), base = c(mape.fc, crps.fc), reconc = c(mape.reconc, crps.reconc)
)
metrics

xtrain = seq(as.Date("1969-01-01"), as.Date("1973-12-01"), by = "months")
xtest = seq(as.Date("1974-01-01"), as.Date("1974-11-01"), by = "months")

df.test = data.frame(time = append(xtrain, xtest), y = as.numeric(HydePark.ts), lo = NA, hi = NA, t = "actual")
df.fc.pf = data.frame(time = xtest, y = unlist(y.hat), lo = unlist(y.hat.lo90), hi = unlist(y.hat.hi90), t = "base")
df.reconc.pf = data.frame(time = xtest, y = unlist(y.reconc), lo = unlist(y.reconc.lo90), hi = unlist(y.reconc.hi90), t = "reconc")

ggplot(df.test, aes(x = time, y = y)) + geom_line() +
  geom_line(aes(x = time, y = y), data = df.fc.pf, colour = "#ff8200") +
  geom_ribbon(aes(ymin = lo, ymax = hi), data = df.fc.pf, fill = alpha("#ff8200", 0.3)) +
  geom_line(aes(x = time, y = y), data = df.reconc.pf, colour = "#0082ff") +
  geom_ribbon(aes(ymin = lo, ymax = hi), data = df.reconc.pf, fill = alpha("#0082ff", 0.3)
)
