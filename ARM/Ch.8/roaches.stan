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
  matrix[N, 3] x = [roach1', treatment', senior']';
}
parameters {
  real alpha;
  vector[3] beta;
}
model {
  y ~ poisson_log_glm(x, alpha + log_expo, beta);
}
