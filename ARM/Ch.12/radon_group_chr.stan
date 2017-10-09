data {
  int<lower=1> N;
  int<lower=1> J; # number of counties
  int<lower=1,upper=J> county[N];
  vector[N] u;
  vector[N] x;
  vector[N] y;
}
parameters {
  vector[2] beta;
  vector[J] eta;
  real mu_b;
  real<lower=0> sigma;
  real<lower=0> sigma_b;
}
transformed parameters {
  vector[J] b;
  vector[N] y_hat;

  b = mu_b + sigma_b * eta;

  for (i in 1:N)
    y_hat[i] = b[county[i]] + x[i] * beta[1] + u[i] * beta[2];
}
model {
  mu_b ~ normal(0, 1);
  eta ~ normal(0, 1);
  beta ~ normal(0, 100);
  sigma ~ cauchy(0, 2.5);
  sigma_b ~ cauchy(0, 2.5);

  y ~ normal(y_hat, sigma);
}
