data {
  int<lower=0> N; // items
  array[N] int<lower=0> K; // initial trials
  array[N] int<lower=0> y; // initial successes
}
parameters {
  real mu; // population mean of success log-odds
  real<lower=0> sigma; // population sd of success log-odds
  vector[N] alpha; // success log-odds
}
model {
  mu ~ normal(-1, 1); // hyperprior
  sigma ~ normal(0, 1); // hyperprior
  alpha ~ normal(mu, sigma); // prior (hierarchical)
  y ~ binomial_logit(K, alpha); // likelihood
}
generated quantities {
  vector[N] theta = inv_logit(alpha);
}
