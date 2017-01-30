data {
  int<lower=1> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  vector[2] beta;
  real<lower=0> sigma;
}
model {
  sigma ~ cauchy(0, 2.5);
  y ~ normal(beta[1] + beta[2] * x, sigma);
}
