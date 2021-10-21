data {
  int<lower=0> N;
  matrix[N, 1] x;
  vector[N] y;
}
parameters {
  real alpha;
  vector[1] beta;
  real<lower=0> sigma;
}
model {
  y ~ normal_id_glm(x, alpha, beta, sigma);
}
