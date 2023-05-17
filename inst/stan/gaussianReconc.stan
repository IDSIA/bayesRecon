data {
  //number of bottom and upper time series
  int<lower=0> n_b;
  int<lower=0> n_a;

  //base: bottom predictives
  vector[n_b] bottom_means_base;
  vector<lower=0>[n_b] bottom_sigma_base;

  //base: aggr predictives
  vector[n_a] aggr_means_base;
  vector<lower=0>[n_a] aggr_sigma_base;

  int sum_a; //sum(A)
  int flat_array[sum_a];  // indexes of bottom to be summed
  int start[n_a];  // starting indexes of bottom to be summed
  int stop[n_a];  // stopping indexes of bottom to be summed
}


parameters {
  vector[n_b] bottom;
}


transformed parameters {
  vector[n_a] aggr;

  for (i in 1:n_a) {
    aggr[i] = sum(bottom[flat_array[start[i]:stop[i]]]);
  }
}

model {
  //gaussian case
  //predittive bottom, vettorializzate
  bottom  ~ normal(bottom_means_base, bottom_sigma_base);
  //predittive aggr, vettorializzate
  target += normal_lpdf(aggr | aggr_means_base, aggr_sigma_base);
}



