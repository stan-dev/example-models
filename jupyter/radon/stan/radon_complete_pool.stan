data {
  int<lower=1> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  sigma ~ normal(0, 5);
  y ~ normal(alpha + beta * x, sigma);
}

