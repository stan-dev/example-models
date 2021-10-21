data {
  int<lower=0> N;
  vector[N] y;
  vector[N] y_lag;
}
transformed data {
  matrix[N, 1] x = [y_lag']';
}
parameters {
  real alpha;
  vector[1] beta;
  real<lower=0> sigma;
}
model {
  y ~ normal_id_glm(x, alpha, beta, sigma);
}
