data {
  int<lower=0> N;
  int<lower=0> y[N];              // count outcomes
  vector<lower=0>[N] E;           // exposure
}
transformed data {
  vector[N] log_E = log(E);
}
parameters {
  real beta0;                // intercept
}
model {
  y ~ poisson_log(log_E + beta0);  // intercept only, no covariates
  beta0 ~ normal(0.0, 2.5);
}
generated quantities {
  vector[N] eta = log_E + beta0;
  vector[N] mu = exp(eta);
}
