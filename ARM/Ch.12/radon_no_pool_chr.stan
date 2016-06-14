data {
  int<lower=0> N;
  int<lower=1,upper=85> county[N];
  vector[N] x;
  vector[N] y;
}
parameters {
  real beta;
  vector[85] eta;
  real mu_a;
  real<lower=0> sigma_a;
  real<lower=0> sigma_y;
}
transformed parameters {
  vector[85] a;
  vector[N] y_hat;

  a <- mu_a + sigma_a * eta;

  for (i in 1:N)
    y_hat[i] <- beta * x[i] + a[county[i]];
}
model {
  beta ~ normal(0, 1);
  mu_a ~ normal(0, 1);
  eta ~ normal(0, 1);
  sigma_a ~ cauchy(0, 2.5);
  sigma_y ~ cauchy(0, 2.5);

  y ~ normal(y_hat, sigma_y);
}
