data {
  int<lower=1> N;  // observations
  int<lower=1> J;  // counties
  array[N] int<lower=1, upper=J> county;
  vector[N] y;     // radon
  vector[N] x;     // floor
  vector[N] u;     // uranium
}
parameters {
  vector[J] alpha;
  vector[2] beta;
  real<lower=0> sigma;
}
model {
  y ~ normal(alpha[county] + beta[1] * x + beta[2] * u, sigma);
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  sigma ~ normal(0, 10);
}
generated quantities {
  array[N] real y_rep = normal_rng(alpha[county] + beta[1] * x + beta[2] * u, sigma);
}
