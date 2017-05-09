data {
  int<lower=0> N;
  int<lower=0> y[N];
  vector[N] x;
}
transformed data {
  vector[N] x_nonzero = x + 0.001;
}
parameters {
  real beta_1;
  real beta_2;
  vector[N] re;  // individual-level random effect
  real<lower=0> sigma;   // scale of random effect
}
model {
  y ~ poisson_log(beta_1 + beta_2 * x_nonzero + re);
  beta_1 ~ normal(0, 10);
  beta_2 ~ normal(0, 2.5);
  re ~ normal(0, sigma);
  sigma ~ normal(0, 5);
}
generated quantities {
  vector[N] mu = exp(beta_1 + beta_2 * x_nonzero + re);
}

