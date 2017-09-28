data {
  int<lower=1> N;
  int<lower=1> J; # number of counties
  int<lower=1,upper=J> county[N];
  vector[N] u;
  vector[N] x;
  vector[N] y;
}
parameters {
  vector[J] alpha;
  vector[2] beta;
  real mu_alpha;
  real mu_beta;
  real<lower=0> sigma;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
}
model {
  vector[N] y_hat;
  for (i in 1:N)
    y_hat[i] = alpha[county[i]] + x[i] * beta[1] + u[i] * beta[2];

  alpha ~ normal(mu_alpha, sigma_alpha);
  beta ~ normal(mu_beta, sigma_beta);
  sigma ~ cauchy(0, 2.5);
  mu_alpha ~ normal(0, 1);
  sigma_alpha ~ cauchy(0, 2.5);
  mu_beta ~ normal(0, 1);
  sigma_beta ~ cauchy(0, 2.5);

  y ~ normal(y_hat, sigma);
}
