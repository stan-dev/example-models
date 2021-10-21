// BUGS pump example (Vol 1, Example 2)
data {
  int<lower=0> N;
  array[N] int<lower=0> x;
  vector[N] t;
}
parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  vector<lower=0>[N] theta;
}
model {
  alpha ~ exponential(1.0);
  beta ~ gamma(0.1, 1.0);
  theta ~ gamma(alpha, beta);
  x ~ poisson(theta .* t);
}
