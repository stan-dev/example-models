## ---- irt-1pl-vague-stan ----
data {
  int<lower=0> I;
  int<lower=0> J;
  int<lower=0,upper=1> y[I,J];
  real mu_theta;
  real<lower=0> sigma_theta;
  real mu_b;
  real<lower=0> sigma_b;
}
parameters {
  vector[I] b;
  vector[J] theta;
}
model {
  theta ~ normal(mu_theta, sigma_theta);
  b ~ normal(mu_b, sigma_b);
  for (i in 1:I)
    y[i] ~ bernoulli_logit(theta - b[i]);
}
