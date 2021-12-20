data {
  int<lower=0> N;
  int<lower=0> J;
  vector[N] y;
  array[N] int<lower=0, upper=1> x;
  array[N] int county;
}
parameters {
  array[J] real a;
  real b;
  real mu_a;
  real<lower=0> sigma_y;
  real<lower=0> sigma_a;
}
model {
  a ~ normal(mu_a, sigma_a);
  for (n in 1 : N) {
    y[n] ~ normal(a[county[n]] + b * x[n], sigma_y);
  }
}
