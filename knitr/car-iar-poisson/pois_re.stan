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
  vector[N] theta;           // heterogeneous random effects
  real<lower=0> sigma;       // non-centered re variance 
}
model {
  y ~ poisson_log(log_E + beta0 + theta * sigma);
  beta0 ~ normal(0.0, 2.5);
  theta ~ normal(0, 1);
  sigma ~ normal(0, 5);
}
generated quantities {
  vector[N] eta = log_E + beta0 + theta * sigma;
  vector[N] mu = exp(eta);
}
