data {
  int<lower=0> N;
  vector[N] watched_hat;
  vector[N] y;
}
transformed data {
  matrix[N, 1] x = [watched_hat']';
}
parameters {
  real alpha;
  vector[1] beta;
  real<lower=0> sigma;
}
model {
  y ~ normal_id_glm(x, alpha, beta, sigma);
}
