data {
  int<lower=0> N;
  vector[N] earn;
  vector[N] height;
}
transformed data {
  matrix[N, 1] x = [height']';
}
parameters {
  real alpha;
  vector[1] beta;
  real<lower=0> sigma;
}
model {
  earn ~ normal_id_glm(x, alpha, beta, sigma);
}
