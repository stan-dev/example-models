data {
  int<lower=1> N;  // observations
  int<lower=1> J;  // counties
  array[N] int<lower=1, upper=J> county;
  vector[N] x;     // floor
  vector[N] y;     // radon
}
parameters {
  vector[J] alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  y ~ normal(alpha[county] + beta * x, sigma);  
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  sigma ~ normal(0, 10);
}
generated quantities {
  array[N] real y_rep = normal_rng(alpha[county] + beta * x, sigma);
}
