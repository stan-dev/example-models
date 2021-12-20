data {
  int<lower=1> K;
  int<lower=1> J;
  int<lower=0> N;
  array[N] vector[J] x;
  array[N] vector[K] y;
}
parameters {
  matrix[K, J] beta;
  cov_matrix[K] Sigma;
}
model {
  array[N] vector[K] mu;
  for (n in 1 : N) {
    mu[n] = beta * x[n];
  }
  y ~ multi_normal(mu, Sigma);
}
