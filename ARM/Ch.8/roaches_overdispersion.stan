data {
  int<lower=0> N;
  vector[N] exposure2;
  vector[N] roach1;
  vector[N] senior;
  vector[N] treatment;
  array[N] int y;
}
transformed data {
  vector[N] log_expo = log(exposure2);
  matrix[N, 3] x = [roach1', senior', treatment']';
}
parameters {
  real alpha;
  vector[3] beta;
  real<lower=0> sigma;
  vector<multiplier=sigma>[N] lambda;
}
model {
  lambda ~ normal(0, sigma);
  y ~ poisson_log_glm(x, alpha + lambda + log_expo, beta);
}
