data {
  int<lower=0> N;
  int<lower=0> y[N];
  vector<lower=0>[N] x;
}
parameters {
  real beta_1;
  real beta_2;
}
model {
  y ~ poisson_log(beta_1 + beta_2 * x);
  beta_1 ~ normal(0,10);
  beta_2 ~ normal(0,2.5);
}
generated quantities {
  vector[N] mu = exp(beta_1 + beta_2 * x);
}
