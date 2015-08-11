## ---- irt-1pl-hier-stan ----
data {
  int<lower=0> I;
  int<lower=0> J;
  int<lower=0,upper=1> y[I,J];
}
parameters {
  vector[I] b;
  vector[J] theta;
  real mu_b;
  real<lower=0> sigma_b;
  real<lower=0> sigma_theta;
}
model {
  // hyperpriors
  mu_b ~ normal(0, 5);
  sigma_b ~ cauchy(0, 2);
  sigma_theta ~ cauchy(0, 2);

  // priors
  b ~ normal(mu_b, sigma_b);
  theta ~ normal(0, sigma_theta);

  // likelihood
  for (i in 1:I)
    y[i] ~ bernoulli_logit(theta - b[i]);
}
