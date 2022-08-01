data {
  int<lower=1> N;  // observations
  int<lower=1> J;  // counties
  array[N] int<lower=1, upper=J> county;
  vector[N] x;
  vector[N] y;
}
parameters {
  real beta;
  vector[J] alpha;
  real<lower=0> sigma;
}
model {
  vector[N] eta;
  for (i in 1:N)
    eta[i] = alpha[county[i]] + beta * x[i];
  sigma ~ normal(0, 5);
  y ~ normal(eta, sigma);  
}
