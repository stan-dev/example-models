data {
  int<lower=1> N;
  int<lower=1> J; # number of counties
  int<lower=1,upper=J> county[N];
  vector[N] x;
  vector[N] y;
}
parameters {
  vector[J] a;
  real beta;
  real<lower=0> sigma_a;
  real<lower=0> sigma_y;
  real mu_a;
}
model {
  vector[N] y_hat;
  for (i in 1:N)
    y_hat[i] = beta * x[i] + a[county[i]];

  beta ~ normal(0, 1);
  mu_a ~ normal(0, 1);
  sigma_a ~ cauchy(0, 2.5);
  sigma_y ~ cauchy(0, 2.5);

  a ~ normal(mu_a, sigma_a);
  y ~ normal(y_hat, sigma_y);
}
